#include "ImplicitPixie3dApplication.h"

#include "SAMRAI/tbox/IEEE.h"
#include <algorithm>

namespace SAMRAI{

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
  
  d_current_dt      = d_initial_dt;
  d_old_dt          = 0.0;
  
  d_first_step      = true;
  
  // the initial condition vector is assumed to have all the necessary information
  // about components
  tbox::Pointer< solv::SAMRAIVectorReal<double> > ic_vector = this->get_ic();
  
#ifdef DEBUG_CHECK_ASSERTIONS
  assert(!ic_vector.isNull());
#endif
  
  // allocate components and space for the current and previous solution vectors
  d_currentSolutionVector = ic_vector->cloneVector("CurrentSolutionVector");
  d_previousSolutionVector = ic_vector->cloneVector("PreviousSolutionVector");
  d_scratchVector = ic_vector->cloneVector("ScratchVector");
  
#ifdef DEBUG_CHECK_ASSERTIONS
  assert(!d_currentSolutionVector.isNull());
  assert(!d_previousSolutionVector.isNull());
  assert(!d_scratchVector.isNull());
#endif
  
  // Initialize them to have the same values as the initial condition vector
  d_currentSolutionVector->copyVector(ic_vector);
  d_previousSolutionVector->copyVector(ic_vector);
  
}

ImplicitPixie3dApplication::~ImplicitPixie3dApplication()
{
}

int
ImplicitPixie3dApplication::evaluateNonlinearFunction(tbox::Pointer< solv::SAMRAIVectorReal<double> > x,
						      tbox::Pointer< solv::SAMRAIVectorReal<double> > f)
{
  tbox::Pointer< solv::SAMRAIVectorReal<double> > nullVector;

  // compute f(u^{n+1})
  this->apply(nullVector, x, f, 1.0, 0.0);

  if (d_first_step) 
    {
      // if it's the first step use backward Euler instead of BDF2
      d_scratchVector->subtract(x, d_currentSolutionVector);
      f->axpy(-d_current_dt, f, d_scratchVector);
    }
  else
    {
      
      const double factor1 = pow(d_current_dt+d_old_dt,2)/(d_old_dt*(2.0*d_current_dt+d_old_dt));
      const double factor2 = pow(d_current_dt,2)/(d_old_dt*(2.0*d_current_dt+d_old_dt));
      const double factor3 = d_current_dt*(d_current_dt+d_old_dt)/(2.0*d_current_dt+d_old_dt);
      
      d_scratchVector->linearSum(-factor1, d_currentSolutionVector, factor2, d_previousSolutionVector);
      d_scratchVector->add(x, d_scratchVector);

      f->axpy(-factor3, f, d_scratchVector);

    }

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
  abort();
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
