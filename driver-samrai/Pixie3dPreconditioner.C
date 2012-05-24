#include "SAMRAI/geom/CartesianGridGeometry.h"
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/pdat/FaceData.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/tbox/Timer.h"
#include "SAMRAI/tbox/TimerManager.h"

#include "source/AMRUtilities.h"

#include "operators/MultilevelOperatorParameters.h"

#include "Pixie3dPreconditioner.h"

namespace SAMRAI{
namespace SAMRSolvers{

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
Pixie3dPreconditioner::getFromInput( tbox::Pointer<tbox::Database> &db,
				     bool is_from_restart)
{
#ifdef DEBUG_CHECK_ASSERTIONS
  //   assert(!db.isNull());
#endif
   if (!is_from_restart)
   {
     
   }
   else
   {
     // currently there is nothing to be read in
     if(db.isNull())
       {
	 tbox::pout << "WARNING::currently a NULL database is being passed to the preconditioner" << std::endl;
       }
   }
}


void
Pixie3dPreconditioner::initializeOperators(tbox::Pointer<tbox::Database> &db)
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
	      tbox::Pointer<tbox::Database> MOperatorDB = db->getDatabase(MOperatorNames[i]);

#ifdef DEBUG_CHECK_ASSERTIONS
	      assert(!MOperatorDB.isNull());
#endif

	      tbox::Pointer<SAMRSolvers::MultilevelOperatorParameters> mlParams ( new SAMRSolvers::MultilevelOperatorParameters(MOperatorDB) );

	      d_MOperators[i].reset( new PCDiagonalMultilevelOperator(mlParams) );
	      
#ifdef DEBUG_CHECK_ASSERTIONS
	      assert(!d_MOperators[i].isNull());
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
      tbox::Pointer<tbox::Database> UOperatorDB = db->getDatabase("UOperator");
      
#ifdef DEBUG_CHECK_ASSERTIONS
      assert(!UOperatorDB.isNull());
#endif
      
      tbox::Pointer<SAMRSolvers::MultilevelOperatorParameters> mlParams ( new SAMRSolvers::MultilevelOperatorParameters(UOperatorDB) );
      
      d_UOperator.reset( new PCDiagonalMultilevelOperator(mlParams) );
      
#ifdef DEBUG_CHECK_ASSERTIONS
      assert(!d_UOperator.isNull());
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
      tbox::Pointer<tbox::Database> LOperatorDB = db->getDatabase("LOperator");
      
#ifdef DEBUG_CHECK_ASSERTIONS
      assert(!LOperatorDB.isNull());
#endif
      
      tbox::Pointer<SAMRSolvers::MultilevelOperatorParameters> mlParams ( new SAMRSolvers::MultilevelOperatorParameters(LOperatorDB) );
      
      d_LOperator.reset( new PCDiagonalMultilevelOperator(mlParams) );
      
#ifdef DEBUG_CHECK_ASSERTIONS
      assert(!d_LOperator.isNull());
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
      tbox::Pointer<tbox::Database> PSchurOperatorDB = db->getDatabase("PSchurOperator");
      
#ifdef DEBUG_CHECK_ASSERTIONS
      assert(!PSchurOperatorDB.isNull());
#endif
      
      tbox::Pointer<SAMRSolvers::MultilevelOperatorParameters> mlParams ( new SAMRSolvers::MultilevelOperatorParameters(PSchurOperatorDB) );
      
      d_PSchurOperator.reset( new PCDiagonalMultilevelOperator(mlParams) );
      
#ifdef DEBUG_CHECK_ASSERTIONS
      assert(!d_PSchurOperator.isNull());
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
Pixie3dPreconditioner::initializeSolvers(tbox::Pointer<tbox::Database> &db)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!d_hierarchy.isNull());
#endif

}

void
Pixie3dPreconditioner::setRefinementBoundaryInterpolant(RefinementBoundaryInterpolation *cf_interpolant)
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
   tbox::Pointer< solv::SAMRAIVectorReal<double> > r,
   tbox::Pointer< solv::SAMRAIVectorReal<double> > z)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!r.isNull());
   assert(!z.isNull());
#endif

   // First try out the identity preconditioner
   z->copyVector(r);

   return(0);
}

int
Pixie3dPreconditioner::setupPreconditioner( PreconditionerParameters *parameters )
{
   static tbox::Pointer<tbox::Timer> t_setup_pc = tbox::TimerManager::getManager()->getTimer("rd2t::Pixie3dPreconditioner::setupPreconditioner");
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
   tbox::Pointer<hier::CoarsenOperator > soln_coarsen_op;
   tbox::Pointer< hier::Variable > cvar;
   tbox::Pointer<hier::GridGeometry > grid_geometry = d_hierarchy->getGridGeometry();

   hier::VariableDatabase* variable_db = hier::VariableDatabase::getDatabase();

   variable_db->mapIndexToVariable(var_id, cvar);
   soln_coarsen_op = grid_geometry->lookupCoarsenOperator(cvar, coarsen_op_str);

   coarsen_var.registerCoarsen(var_id, var_id, soln_coarsen_op);

   // coarsen data to be consistent across all levels
   for ( int ln = d_hierarchy->getFinestLevelNumber();
         ln >= 0;
         ln-- )
   {

      tbox::Pointer<hier::PatchLevel > level = d_hierarchy->getPatchLevel(ln);

      if(ln<d_hierarchy->getFinestLevelNumber())
      {
         tbox::Pointer< xfer::CoarsenSchedule > coarsen_schedule = getCoarsenSchedule(ln, coarsen_var);
         coarsen_schedule->coarsenData();
      }
   }
}

tbox::Pointer<xfer::CoarsenSchedule >
Pixie3dPreconditioner::getCoarsenSchedule(int ln, xfer::CoarsenAlgorithm &crs_fill_alg)
{
   tbox::Pointer<xfer::CoarsenSchedule > crs_fill_schedule;

   tbox::Pointer<hier::PatchLevel > level = d_hierarchy->getPatchLevel(ln);
   tbox::Pointer<hier::PatchLevel > flevel = d_hierarchy->getPatchLevel(ln+1);

   crsList* coarsen_fill_schedules = getCoarsenSchedules(ln);

   bool found_compatible_schedule = false;

   for (crsList::Iterator crs(*coarsen_fill_schedules); crs; crs++)
   {
      if (crs_fill_alg.checkConsistency(crs()))
      {
         crs_fill_alg.resetSchedule(crs());
         crs_fill_schedule = crs();
         found_compatible_schedule = true;
      }

      if (found_compatible_schedule) break;
   }

   if (!found_compatible_schedule)
   {
      crs_fill_schedule = crs_fill_alg.createSchedule(level, flevel);
      coarsen_fill_schedules->appendItem(crs_fill_schedule);
   }

   return crs_fill_schedule;

}

}
}
