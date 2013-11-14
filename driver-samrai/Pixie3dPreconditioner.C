#include "SAMRAI/geom/CartesianGridGeometry.h"
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/pdat/FaceData.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/tbox/Timer.h"
#include "SAMRAI/tbox/TimerManager.h"

#include "source/AMRUtilities.h"

#include "operators/MultilevelOperatorParameters.h"

#include"PCComponentFACSolver.h"
#include "Pixie3dPreconditioner.h"

namespace SAMRAI{
namespace Pixie3d{

Pixie3dPreconditioner::Pixie3dPreconditioner(Pixie3dPreconditionerParameters *parameters)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(parameters!=NULL);
#endif

   d_preconditioner_print_flag   = false;

   d_hierarchy                   = parameters->d_hierarchy;

   d_cf_interpolant              = parameters->d_cf_interpolant;

   d_dt                          = 0.0;

   getFromInput(parameters->d_db);

   initializeOperators(parameters->d_db);

   initializeSolvers(parameters->d_db);
}

Pixie3dPreconditioner::~Pixie3dPreconditioner()
{

}

void
Pixie3dPreconditioner::getFromInput( boost::shared_ptr<tbox::Database> &db,
				     bool is_from_restart)
{
   if (!is_from_restart)
   {
     
   }
   else
   {
     // currently there is nothing to be read in
     if(db==NULL)
       {
	 tbox::pout << "WARNING::currently a NULL database is being passed to the preconditioner" << std::endl;
       }
   }
}


void
Pixie3dPreconditioner::initializeOperators(boost::shared_ptr<tbox::Database> &db)
{
  
  // first read in the names of the M Operators and construct the M Operators
  if (db->keyExists("MOperators"))
    {
      tbox::Array<std::string> MOperatorNames = db->getStringArray("MOperators");

      d_MOperators.resizeArray(MOperatorNames.getSize());
      
      for(int i=0; i<MOperatorNames.getSize(); i++)
	{
	  
	  if (db->keyExists(MOperatorNames[i]))
	    {
	      boost::shared_ptr<tbox::Database> MOperatorDB = db->getDatabase(MOperatorNames[i]);

#ifdef DEBUG_CHECK_ASSERTIONS
	      assert(MOperatorDB!=NULL);
#endif

	      boost::shared_ptr<SAMRSolvers::MultilevelOperatorParameters> mlParams ( new SAMRSolvers::MultilevelOperatorParameters(MOperatorDB) );

	      mlParams->d_hierarchy      = d_hierarchy;
	      mlParams->d_cf_interpolant = d_cf_interpolant;
	      mlParams->d_set_boundary_ghosts.reset();

	      d_MOperators[i].reset( new PCDiagonalMultilevelOperator(mlParams) );
	      
#ifdef DEBUG_CHECK_ASSERTIONS
	      assert(d_MOperators[i]!=NULL);
#endif

	    }
	  else
	    {
	      TBOX_ERROR( "Pixie3dPreconditioner"
			  << " -- Required key "  << MOperatorNames[i]
			  << " missing in input.");
	    }
	}
    }
  else
    {
      TBOX_ERROR( "Pixie3dPreconditioner"
		  << " -- Required key `MOperators'"
		  << " missing in input.");
    }

  // next construct the U Operator
  if (db->keyExists("UOperator"))
    {
      boost::shared_ptr<tbox::Database> UOperatorDB = db->getDatabase("UOperator");
      
#ifdef DEBUG_CHECK_ASSERTIONS
      assert(UOperatorDB!=NULL);
#endif
      
      boost::shared_ptr<SAMRSolvers::MultilevelOperatorParameters> mlParams ( new SAMRSolvers::MultilevelOperatorParameters(UOperatorDB) );
      
      mlParams->d_hierarchy      = d_hierarchy;
      mlParams->d_cf_interpolant = d_cf_interpolant;
      mlParams->d_set_boundary_ghosts.reset();
      
      d_UOperator.reset( new PCDiagonalMultilevelOperator(mlParams) );
      
#ifdef DEBUG_CHECK_ASSERTIONS
      assert(d_UOperator!=NULL);
#endif      
      
    }  
  else
    {
      TBOX_ERROR( "Pixie3dPreconditioner"
		  << " -- Required key `UOperator'"
		  << " missing in input.");
    }
  
  // next construct the L Operator
  if (db->keyExists("LOperator"))
    {
      boost::shared_ptr<tbox::Database> LOperatorDB = db->getDatabase("LOperator");
      
#ifdef DEBUG_CHECK_ASSERTIONS
      assert(LOperatorDB!=NULL);
#endif
      
      boost::shared_ptr<SAMRSolvers::MultilevelOperatorParameters> mlParams ( new SAMRSolvers::MultilevelOperatorParameters(LOperatorDB) );
      
      mlParams->d_hierarchy      = d_hierarchy;
      mlParams->d_cf_interpolant = d_cf_interpolant;
      mlParams->d_set_boundary_ghosts.reset();
      
      d_LOperator.reset( new PCDiagonalMultilevelOperator(mlParams) );
      
#ifdef DEBUG_CHECK_ASSERTIONS
      assert(d_LOperator!=NULL);
#endif      
      
    }  
  else
    {
      TBOX_ERROR( "Pixie3dPreconditioner"
		  << " -- Required key `LOperator'"
		  << " missing in input.");
    }
  

  // next construct the PSchur Operator
  if (db->keyExists("PSchurOperator"))
    {
      boost::shared_ptr<tbox::Database> PSchurOperatorDB = db->getDatabase("PSchurOperator");
      
#ifdef DEBUG_CHECK_ASSERTIONS
      assert(PSchurOperatorDB!=NULL);
#endif
      
      boost::shared_ptr<SAMRSolvers::MultilevelOperatorParameters> mlParams ( new SAMRSolvers::MultilevelOperatorParameters(PSchurOperatorDB) );
      
      mlParams->d_hierarchy      = d_hierarchy;
      mlParams->d_cf_interpolant = d_cf_interpolant;
      mlParams->d_set_boundary_ghosts.reset();
      
      d_PSchurOperator.reset( new PCDiagonalMultilevelOperator(mlParams) );
      
#ifdef DEBUG_CHECK_ASSERTIONS
      assert(d_PSchurOperator!=NULL);
#endif      
      
    }  
  else
    {
      TBOX_ERROR( "Pixie3dPreconditioner"
		  << " -- Required key `PSchurOperator'"
		  << " missing in input.");
    }

}

void
Pixie3dPreconditioner::initializeSolvers(boost::shared_ptr<tbox::Database> &db)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(d_hierarchy!=NULL);
#endif
  // first read in the names of the M Operators and construct the M Operators
  if (db->keyExists("MSolvers"))
    {
      tbox::Array<std::string> MSolverNames = db->getStringArray("MSolvers");

      d_MSolvers.resizeArray(MSolverNames.getSize());
      
      for(int i=0; i<MSolverNames.getSize(); i++)
	{
	  
	  if (db->keyExists(MSolverNames[i]))
	    {
	      boost::shared_ptr<tbox::Database> MSolverDB = db->getDatabase(MSolverNames[i]);

#ifdef DEBUG_CHECK_ASSERTIONS
	      assert(MSolverDB!=NULL);
#endif

	      SAMRSolvers::MultilevelSolverParameters *mlParams = new SAMRSolvers::MultilevelSolverParameters(MSolverDB) ;

	      mlParams->d_hierarchy                = d_hierarchy;

	      d_MSolvers[i].reset( new PCComponentFACSolver(mlParams) );

	      delete mlParams;
	      
#ifdef DEBUG_CHECK_ASSERTIONS
	      assert(d_MSolvers[i]!=NULL);
#endif

	    }
	  else
	    {
	      TBOX_ERROR( "Pixie3dPreconditioner"
			  << " -- Required key "  << MSolverNames[i]
			  << " missing in input.");
	    }
	}
    }
  else
    {
      TBOX_ERROR( "Pixie3dPreconditioner"
		  << " -- Required key `MSolvers'"
		  << " missing in input.");
    }

  // next construct the PSchur Solver
  if (db->keyExists("PSchurSolver"))
    {
      boost::shared_ptr<tbox::Database> PSchurSolverDB = db->getDatabase("PSchurSolver");
      
#ifdef DEBUG_CHECK_ASSERTIONS
      assert(PSchurSolverDB!=NULL);
#endif
      
      SAMRSolvers::MultilevelSolverParameters *mlParams = new SAMRSolvers::MultilevelSolverParameters(PSchurSolverDB) ;
      
      mlParams->d_hierarchy                = d_hierarchy;
      
      d_PSchurSolver.reset( new PCComponentFACSolver(mlParams) );
      
#ifdef DEBUG_CHECK_ASSERTIONS
      assert(d_PSchurSolver!=NULL);
#endif      

      delete mlParams;
      
    }  
  else
    {
      TBOX_ERROR( "Pixie3dPreconditioner"
		  << " -- Required key `PSchurSolver'"
		  << " missing in input.");
    }

}

void
Pixie3dPreconditioner::setRefinementBoundaryInterpolant(boost::shared_ptr<RefinementBoundaryInterpolation> cf_interpolant)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(cf_interpolant!=NULL);

#endif

   d_cf_interpolant = cf_interpolant;

}

/*
*************************************************************************
*                                                                       *
* Apply the FAC Preconditioner to the right-hand side r, placing the    *
* result in the vector z.                                               *
*                                                                       *
*************************************************************************
*/
int
Pixie3dPreconditioner::applyPreconditioner(
   boost::shared_ptr< solv::SAMRAIVectorReal<double> > r,
   boost::shared_ptr< solv::SAMRAIVectorReal<double> > z)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(r!=NULL);
   assert(z!=NULL);
#endif

   // First try out the identity preconditioner
   z->copyVector(r);

   return(0);
}

int
Pixie3dPreconditioner::setupPreconditioner( SAMRSolvers::PreconditionerParameters *parameters )
{
   static boost::shared_ptr<tbox::Timer> t_setup_pc = tbox::TimerManager::getManager()->getTimer("rd2t::Pixie3dPreconditioner::setupPreconditioner");
   t_setup_pc->start();

   Pixie3dPreconditionerParameters *ptr = dynamic_cast<Pixie3dPreconditionerParameters *>(parameters);

#ifdef DEBUG_CHECK_ASSERTIONS
   assert(ptr!=NULL);
#endif

   t_setup_pc->stop();

   return(0);
}

void
Pixie3dPreconditioner::coarsenVariable( const int var_id,
                                     std::string coarsen_op_str)
{
  xfer::CoarsenAlgorithm coarsen_var(d_hierarchy->getDim());
   boost::shared_ptr<hier::CoarsenOperator > soln_coarsen_op;
   boost::shared_ptr< hier::Variable > cvar;
   boost::shared_ptr<geom::CartesianGridGeometry> grid_geometry = 
      boost::dynamic_pointer_cast<geom::CartesianGridGeometry>(d_hierarchy->getGridGeometry());

   hier::VariableDatabase* variable_db = hier::VariableDatabase::getDatabase();

   variable_db->mapIndexToVariable(var_id, cvar);
   soln_coarsen_op = grid_geometry->lookupCoarsenOperator(cvar, coarsen_op_str);

   coarsen_var.registerCoarsen(var_id, var_id, soln_coarsen_op);

   // coarsen data to be consistent across all levels
   for ( int ln = d_hierarchy->getFinestLevelNumber();
         ln >= 0;
         ln-- )
   {

      boost::shared_ptr<hier::PatchLevel > level = d_hierarchy->getPatchLevel(ln);

      if(ln<d_hierarchy->getFinestLevelNumber())
      {
         boost::shared_ptr< xfer::CoarsenSchedule > coarsen_schedule = getCoarsenSchedule(ln, coarsen_var);
         coarsen_schedule->coarsenData();
      }
   }
}

boost::shared_ptr<xfer::CoarsenSchedule >
Pixie3dPreconditioner::getCoarsenSchedule(int ln, xfer::CoarsenAlgorithm &crs_fill_alg)
{
   boost::shared_ptr<xfer::CoarsenSchedule > crs_fill_schedule;

   boost::shared_ptr<hier::PatchLevel > level = d_hierarchy->getPatchLevel(ln);
   boost::shared_ptr<hier::PatchLevel > flevel = d_hierarchy->getPatchLevel(ln+1);

   crsList* coarsen_fill_schedules = getCoarsenSchedules(ln);

   bool found_compatible_schedule = false;

   for (crsList::iterator crs=coarsen_fill_schedules->begin(); crs!=coarsen_fill_schedules->end(); crs++)
   {
      if (crs_fill_alg.checkConsistency(*crs))
      {
         crs_fill_alg.resetSchedule(*crs);
         crs_fill_schedule = *crs;
         found_compatible_schedule = true;
      }

      if (found_compatible_schedule) break;
   }

   if (!found_compatible_schedule)
   {
      crs_fill_schedule = crs_fill_alg.createSchedule(level, flevel);
      coarsen_fill_schedules->push_back(crs_fill_schedule);
   }

   return crs_fill_schedule;

}

}
}
