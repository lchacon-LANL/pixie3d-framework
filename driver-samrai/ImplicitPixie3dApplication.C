#include "ImplicitPixie3dApplication.h"
#include "source/AMRUtilities.h"

#include "SAMRAI/tbox/IEEE.h"
#include "SAMRAI/solv/PETSc_SAMRAIVectorReal.h"
#include <algorithm>

namespace SAMRAI{

ImplicitPixie3dApplication::ImplicitPixie3dApplication():
  pixie3dApplication()
{
  /*
   * Set default values for data members.
   */
  
  d_current_time                = tbox::IEEE::getSignalingNaN();
  d_current_dt                  = tbox::IEEE::getSignalingNaN();
  d_old_dt                      = tbox::IEEE::getSignalingNaN();
  d_new_time                    = tbox::IEEE::getSignalingNaN();
  d_initial_dt                  = tbox::IEEE::getSignalingNaN();

  d_debug_print_info_level=0;

  d_newSolutionVector.setNull();
  d_currentSolutionVector.setNull();
  d_previousSolutionVector.setNull();
  d_scratchVector.setNull();
}
  
ImplicitPixie3dApplication::ImplicitPixie3dApplication(ImplicitPixie3dApplicationParameters *parameters):
  pixie3dApplication(parameters)
{  
  /*
   * Set default values for data members.
   */
  
  d_current_time                = tbox::IEEE::getSignalingNaN();
  d_current_dt                  = tbox::IEEE::getSignalingNaN();
  d_old_dt                      = tbox::IEEE::getSignalingNaN();
  d_new_time                    = tbox::IEEE::getSignalingNaN();
  d_initial_dt                  = tbox::IEEE::getSignalingNaN();
  
  d_debug_print_info_level=0;

  d_newSolutionVector.setNull();
  d_currentSolutionVector.setNull();
  d_previousSolutionVector.setNull();
  d_scratchVector.setNull();

  initialize(parameters);
}
  
void
ImplicitPixie3dApplication::initialize(ImplicitPixie3dApplicationParameters *parameters)
{

  tbox::pout << "Begin: Initializing implicit run" << std::endl;
  
  pixie3dApplication::initialize(parameters);

  tbox::Pointer<tbox::Database> db = parameters->d_db;

#ifdef DEBUG_CHECK_ASSERTIONS
  assert(!db.isNull());
#endif
  
  if (db->keyExists("initial_timestep")) {
    d_initial_dt = db->getDouble("initial_timestep");
  } else {
    TBOX_ERROR(d_object_name << " -- Key data `initial_timestep'"
	       << " missing in input.");
  }
  
  d_current_dt      = d_initial_dt;
  d_old_dt          = 0.0;
  
  d_first_step      = true;
  
  // the initial condition vector is assumed to have all the necessary information
  // about components
  tbox::Pointer< solv::SAMRAIVectorReal<double> > ic_vector = this->get_ic();
  
#ifdef DEBUG_CHECK_ASSERTIONS
  assert(!ic_vector.isNull());
#endif
  
  // Initialize them to have the same values as the initial condition vector
  d_currentSolutionVector->copyVector(ic_vector);
  d_previousSolutionVector->copyVector(ic_vector);

  // this sets the weight_id to be the volume of each cell in cells that do not lie
  // under finer cells and to zero for cells that lie under finer cells. The weight id
  // is important for SAMRAIVectors and is used in computing vector norms correctly
  AMRUtilities::setVectorWeights(d_hierarchy, d_weight_id);

  tbox::pout << "End: Initializing implicit run" << std::endl;

}

void
ImplicitPixie3dApplication::resetHierarchyConfiguration( const tbox::Pointer<hier::PatchHierarchy> hierarchy,
							 const int coarsest_level,
							 const int finest_level )
{
  // Check if the application has been initialized
  if ( d_hierarchy.isNull() )
        return;
  tbox::pout << "Begin: implicit resetHierarchyConfiguration" << std::endl;
  // call the base class function
  pixie3dApplication::resetHierarchyConfiguration(hierarchy, coarsest_level, finest_level);
  
  // the initial condition vector is assumed to have all the necessary information
  // about components
  tbox::Pointer< solv::SAMRAIVectorReal<double> > ic_vector = this->get_ic();
  
#ifdef DEBUG_CHECK_ASSERTIONS
  assert(!ic_vector.isNull());
#endif
  

  if(!d_vectorsCloned)
    {
      // allocate components and space for the current and previous solution vectors
      d_currentSolutionVector = ic_vector->cloneVector("CurrentSolutionVector");
      d_currentSolutionVector->allocateVectorData(0.0);

      d_previousSolutionVector = ic_vector->cloneVector("PreviousSolutionVector");
      d_previousSolutionVector->allocateVectorData(0.0);

      d_scratchVector = ic_vector->cloneVector("ScratchVector");
      d_scratchVector->allocateVectorData(0.0);
      
      d_vectorsCloned = true;
    }
#ifdef DEBUG_CHECK_ASSERTIONS
  assert(!d_currentSolutionVector.isNull());
  assert(!d_previousSolutionVector.isNull());
  assert(!d_scratchVector.isNull());
#endif

  tbox::pout << "Middle: implicit resetHierarchyConfiguration" << std::endl;
  if(d_newSolutionVector)
    {
      d_newSolutionVector->resetLevels(0, hierarchy->getNumberOfLevels()-1);
    }
  
  // reset the level fields in the vectors after regridding
  d_currentSolutionVector->resetLevels(0, hierarchy->getNumberOfLevels()-1);
  tbox::pout << "Middle: implicit resetHierarchyConfiguration" << std::endl;
  d_previousSolutionVector->resetLevels(0, hierarchy->getNumberOfLevels()-1);
  d_scratchVector->resetLevels(0, hierarchy->getNumberOfLevels()-1);
  
  // this sets the weight_id to be the volume of each cell in cells that do not lie
  // under finer cells and to zero for cells that lie under finer cells. The weight id
  // is important for SAMRAIVectors and is used in computing vector norms correctly
  AMRUtilities::setVectorWeights(d_hierarchy, d_weight_id);
  tbox::pout << "End: implicit resetHierarchyConfiguration" << std::endl;
  
}

ImplicitPixie3dApplication::~ImplicitPixie3dApplication()
{
}

/*
*************************************************************************
*                                                                       *
* Intermediate call between the PETSc abstract inteface and SAMRAI.     *
*                                                                       *
*************************************************************************
*/
int ImplicitPixie3dApplication::evaluateNonlinearFunction(Vec xcur, Vec fcur)
{
  tbox::Pointer< solv::SAMRAIVectorReal<double> > x =
    solv::PETSc_SAMRAIVectorReal<double>::getSAMRAIVector(xcur);
  tbox::Pointer< solv::SAMRAIVectorReal<double> > f =
    solv::PETSc_SAMRAIVectorReal<double>::getSAMRAIVector(fcur);
  
  return ( evaluateNonlinearFunction(x, f) );
  
}

/*
*************************************************************************
*                                                                       *
* Intermediate call between the PETSc abstract interface and SAMRAI.    *
*                                                                       *
*************************************************************************
*/
int
ImplicitPixie3dApplication::applyPreconditioner(Vec r, Vec z)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(r != NULL);
   assert(z != NULL);
#endif

   tbox::Pointer< solv::SAMRAIVectorReal<double> > rhs  = solv::PETSc_SAMRAIVectorReal<double>::getSAMRAIVector(r);
   tbox::Pointer< solv::SAMRAIVectorReal<double> > soln = solv::PETSc_SAMRAIVectorReal<double>::getSAMRAIVector(z);

  this->applyPreconditioner(rhs, soln);

   return 0;
}

/*
*************************************************************************
*                                                                       *
* Intermediate call between the PETSc abstract interface and SAMRAI.    *
*                                                                       *
*************************************************************************
*/
int
ImplicitPixie3dApplication::setupPreconditioner(Vec xcur)
{
  tbox::Pointer< solv::SAMRAIVectorReal<double> > soln = solv::PETSc_SAMRAIVectorReal<double>::getSAMRAIVector(xcur);  
  
  return(setupPreconditioner(soln));
  
}
 
int
ImplicitPixie3dApplication::evaluateNonlinearFunction(tbox::Pointer< solv::SAMRAIVectorReal<double> > x,
						      tbox::Pointer< solv::SAMRAIVectorReal<double> > f)
{
  tbox::pout << "Begin: evaluateNonlinearFunction" << std::endl;
  tbox::Pointer< solv::SAMRAIVectorReal<double> > nullVector;

  tbox::pout << "Begin: Calling f(u)" << std::endl;
  // compute f(u^{n+1})
  this->apply(nullVector, x, f, 1.0, 0.0);
  tbox::pout << "End: Calling f(u)" << std::endl;

  if (d_first_step) 
    {
      tbox::pout << "Begin: calculating BE residual" << std::endl;
      // if it's the first step use backward Euler instead of BDF2
      d_scratchVector->subtract(x, d_currentSolutionVector);
      f->axpy(-d_current_dt, f, d_scratchVector);
      tbox::pout << "End: calculating BE residual" << std::endl;
    }
  else
    {
      
      const double factor1 = pow(d_current_dt+d_old_dt,2)/(d_old_dt*(2.0*d_current_dt+d_old_dt));
      const double factor2 = pow(d_current_dt,2)/(d_old_dt*(2.0*d_current_dt+d_old_dt));
      const double factor3 = d_current_dt*(d_current_dt+d_old_dt)/(2.0*d_current_dt+d_old_dt);
      
      tbox::pout << "Begin: calculating BDF2 residual" << std::endl;
      d_scratchVector->linearSum(-factor1, d_currentSolutionVector, factor2, d_previousSolutionVector);
      d_scratchVector->add(x, d_scratchVector);

      f->axpy(-factor3, f, d_scratchVector);
      tbox::pout << "End: calculating BDF2 residual" << std::endl;

    }

  tbox::pout << "End: evaluateNonlinearFunction" << std::endl;

  d_currentSolutionVector->print(tbox::pout);
  x->print(tbox::pout);
  f->print(tbox::pout);
  
  return 0;
}

int
ImplicitPixie3dApplication::setupPreconditioner(tbox::Pointer< solv::SAMRAIVectorReal<double> > x)
{
  abort();
}

int
ImplicitPixie3dApplication::applyPreconditioner(tbox::Pointer< solv::SAMRAIVectorReal<double> > r,
						tbox::Pointer< solv::SAMRAIVectorReal<double> > z)
{
  abort();
}

void
ImplicitPixie3dApplication::setupSolutionVector(tbox::Pointer< solv::SAMRAIVectorReal<double> > solution)
{
#ifdef DEBUG_CHECK_ASSERTIONS
  assert(!solution.isNull());
#endif
  d_newSolutionVector = solution;
  
  const tbox::Dimension &dimObj = d_hierarchy->getDim();
  
  hier::IntVector zeroGhosts(dimObj, 0);

  /* 
   * Define variable required contexts.
   */
  hier::VariableDatabase* variable_db = hier::VariableDatabase::getDatabase();
  
  tbox::Pointer< hier::VariableContext > solverContext      = variable_db->getContext("SolverContext");
  
  // the initial condition vector is assumed to have all the necessary information
  // about components
  tbox::Pointer< solv::SAMRAIVectorReal<double> > ic_vector = this->get_ic();
  
#ifdef DEBUG_CHECK_ASSERTIONS
  assert(!ic_vector.isNull());
#endif

  for(int component=0; component<ic_vector->getNumberOfComponents(); component++)
    {
      // get the variable associated with the component
      tbox::Pointer< hier::Variable > var = ic_vector->getComponentVariable(component);
      
#ifdef DEBUG_CHECK_ASSERTIONS
      assert(!var.isNull());
#endif
      // get the control volume id for the component
      const int controlVolId = ic_vector->getControlVolumeIndex(component);

      // register a new component id using the solver context
      const int componentId = variable_db->registerVariableAndContext(var,
								      solverContext,
								      zeroGhosts);
      // add a new component to the vector
      solution->addComponent(var,
			     componentId,
			     controlVolId);
      
    }

  d_newSolutionVector->allocateVectorData(0.0);
}

double
ImplicitPixie3dApplication::getInitialDt(void)
{
    double returnVal = 0.0;

#if 0    
    if((d_restart_data!=NULL)&& d_restart_data->readTimeStepFromRestart())
        returnVal = d_restart_data->getNextDt();
    else
        returnVal = d_initial_dt;
#else
    returnVal = d_initial_dt;
#endif

    return returnVal;
}

double
ImplicitPixie3dApplication::getNextDt(const bool good_solution, const int solver_retcode)
{
  // setting d_old_dt all the time, BP: 03/23/04
  d_old_dt = d_current_dt;

  if (good_solution) 
    {
    }
  else
    {
      d_current_dt = 0.0;
    }
  
  d_current_dt=std::min(d_current_dt, d_max_timestep);
  
  return(d_current_dt);
}

void
ImplicitPixie3dApplication::setInitialGuess(const bool first_step,
					    const double current_time,
					    const double current_dt,
					    const double old_dt)
{
  d_first_step = first_step;
  
  // the initial condition vector is assumed to have all the necessary information
  // about components
  tbox::Pointer< solv::SAMRAIVectorReal<double> > ic_vector = this->get_ic();
  
#ifdef DEBUG_CHECK_ASSERTIONS
  assert(!ic_vector.isNull());
#endif

  if(d_first_step)
    {
      // Initialize them to have the same values as the initial condition vector
      d_currentSolutionVector->copyVector(ic_vector);
      d_previousSolutionVector->copyVector(ic_vector);
    }

  // for now just copy over from the previous timestep
  d_newSolutionVector->copyVector(d_currentSolutionVector);
  
  tbox::pout << "In the initial guess routine, not doing anything for now" << std::endl;
}
  
bool
ImplicitPixie3dApplication::checkNewSolution(const int solver_retcode)
{
  if (solver_retcode == 1)
    {
      return(true);
    }
  else
    {
      return(false);
    }
}

void
ImplicitPixie3dApplication::updateSolution(const double new_time)
{
  static int counter;
  counter++;
  
  d_new_time = d_current_time = new_time;

  // swap vectors between current and previous so a copy is not incurred
  d_previousSolutionVector->swapVectors(d_currentSolutionVector);

  // copy the new solution to the current solution. A swap is not
  // done so as not to cause side effects in this case
  d_currentSolutionVector->copyVector(d_newSolutionVector);

}
  
void
ImplicitPixie3dApplication::putToDatabase(tbox::Pointer<tbox::Database> db)
{
  abort();
}
  
}
