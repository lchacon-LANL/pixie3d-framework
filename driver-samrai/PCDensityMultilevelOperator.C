
#include "PCDensityMultilevelOperator.h"


#include "SAMRAI/xfer/CoarsenAlgorithm.h"
#include "SAMRAI/geom/CartesianGridGeometry.h"

// SAMRAIUTILS headers
#include "source/AMRUtilities.h"

namespace SAMRAI {
namespace SAMRSolvers {

// default constructor
PCDensityMultilevelOperator::PCDensityMultilevelOperator()
{
   d_flux_id                      = -1;
   d_coarsen_diffusive_fluxes     = true;
   d_schedules_initialized        = false;
   d_adjust_cf_coefficients       = false; 
   d_variable_order_interpolation = false;
   d_reset_ghost_values           = true;
   d_face_coarsen_op_str          = "CONSERVATIVE_COARSEN";
   d_cell_coarsen_op_str          = "CONSERVATIVE_COARSEN";
   d_cell_refine_op_str           = "CONSTANT_REFINE";
   d_face_refine_op_str           = "CONSTANT_REFINE";
   d_flux.setNull();

   const int hierarchy_size       = d_hierarchy->getNumberOfLevels();

   d_level_operators.resizeArray(hierarchy_size);
   d_flux_coarsen_schedule.resizeArray(hierarchy_size-1);
   d_src_coarsen_schedule.resizeArray(hierarchy_size-1);
   d_var_refine_schedule.resizeArray(hierarchy_size);
   d_interpolate_schedule.resizeArray(hierarchy_size);
}

PCDensityMultilevelOperator::PCDensityMultilevelOperator(MultilevelOperatorParameters *parameters):MultilevelOperator(parameters)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(parameters!=NULL);
#endif
   d_flux_id                      = -1;
   d_coarsen_diffusive_fluxes     = true;
   d_schedules_initialized        = false;
   d_adjust_cf_coefficients       = false; 
   d_variable_order_interpolation = false;
   d_reset_ghost_values           = true;
   d_face_coarsen_op_str          = "CONSERVATIVE_COARSEN";
   d_cell_coarsen_op_str          = "CONSERVATIVE_COARSEN";
   d_cell_refine_op_str           = "CONSTANT_REFINE";
   d_face_refine_op_str           = "CONSTANT_REFINE";
   d_flux.setNull();

   const int hierarchy_size       = d_hierarchy->getNumberOfLevels();

   d_level_operators.resizeArray(hierarchy_size);
   d_flux_coarsen_schedule.resizeArray(hierarchy_size-1);
   d_src_coarsen_schedule.resizeArray(hierarchy_size-1);
   d_var_refine_schedule.resizeArray(hierarchy_size);
   d_interpolate_schedule.resizeArray(hierarchy_size);

   // read in parameters from database
   getFromInput(parameters->d_db);
     
   initializeInternalVariableData();

   // make sure the flux id is non zero before initializing the level operators
   // else the level operators will all create individual spaces for the fluxes
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(d_flux_id>=0);
#endif

   initializeLevelOperators();
     
}

PCDensityMultilevelOperator::~PCDensityMultilevelOperator()
{

   for(int ln=0; ln<=d_hierarchy->getFinestLevelNumber(); ln++)
   {
     tbox::Pointer<hier::PatchLevel > level = d_hierarchy->getPatchLevel(ln);

     if(level->checkAllocated(d_flux_id))
       {
	 level->deallocatePatchData(d_flux_id);
       }
   }

   // the next two lines are kept at the bottom of the dtor so that
   // if the remove happens in the levels this does not prevent
   // deallocation of memory on the levels potentially
   hier::VariableDatabase* variable_db = hier::VariableDatabase::getDatabase();
   variable_db->removePatchDataIndex(d_flux_id);
}

void
PCDensityMultilevelOperator::initializeInternalVariableData()
{
   hier::VariableDatabase* variable_db = hier::VariableDatabase::getDatabase();

   const tbox::Dimension &dim = d_hierarchy->getDim();
   const hier::IntVector zero_ghosts = hier::IntVector::getZero(dim);

   const tbox::Pointer< hier::VariableContext > scratch_cxt = variable_db->getContext("SCRATCH");

   std::ostringstream ibuffer;
   ibuffer<<(long)d_object_id;
   std::string object_str=ibuffer.str();

   std::string cellFlux("PCDensityMultilevelOperator_InternalFlux");
   cellFlux+=object_str;

   d_flux = variable_db->getVariable(cellFlux);

   if (!d_flux) 
   {
     d_flux = new pdat::FaceVariable<double>(dim,cellFlux,1);
   }

   d_flux_id = variable_db->registerVariableAndContext(d_flux,
                                                       scratch_cxt,
                                                       zero_ghosts);

   for(int ln=0; ln<d_hierarchy->getNumberOfLevels(); ln++)
   {
      tbox::Pointer<hier::PatchLevel > level = d_hierarchy->getPatchLevel(ln);
      // allocate storage space
      if(!level->checkAllocated(d_flux_id))
      {
         level->allocatePatchData(d_flux_id);
      }

#ifdef DEBUG_CHECK_ASSERTIONS
      assert(level->checkAllocated(d_flux_id));
#endif

   }
}

void 
PCDensityMultilevelOperator::getFromInput(tbox::Pointer<tbox::Database> db)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!db.isNull());
#endif

   if (db->keyExists("tangent_interp_scheme")) 
     {
       d_tangent_interp_scheme_str = db->getString("tangent_interp_scheme");
       d_tangent_interp_scheme = SAMRAI::RefinementBoundaryInterpolation::lookupInterpolationScheme(d_tangent_interp_scheme_str);
   } 
   else 
   {
      TBOX_ERROR( "PCDensityMultilevelOperator" 
                 << " -- Required key `tangent_interp_scheme'"
                 << " missing in input.");
   }

   if (db->keyExists("normal_interp_scheme")) 
   {
     d_normal_interp_scheme_str = db->getString("normal_interp_scheme");
     d_normal_interp_scheme = SAMRAI::RefinementBoundaryInterpolation::lookupInterpolationScheme(d_normal_interp_scheme_str);
   } 
   else 
   {
      TBOX_ERROR( "PCDensityMultilevelOperator" 
                 << " -- Required key `normal_interp_scheme'"
                 << " missing in input.");
   }

   if (db->keyExists("coarsen_diffusive_fluxes")) 
   {
      d_coarsen_diffusive_fluxes = db->getBool("coarsen_diffusive_fluxes");
   } 
   else 
   {
      TBOX_ERROR("PCDensityMultilevelOperator" 
                 << " -- Required key `coarsen_diffusive_fluxes'"
                 << " missing in input.");
   }
   
   if(db->keyExists("face_coarsen_op"))
   {
      d_face_coarsen_op_str = db->getString("face_coarsen_op");
   }
   else
   {
      TBOX_ERROR("PCDensityMultilevelOperator" 
                 << " -- Required key `face_coarsen_op'"
                 << " missing in input.");
   }

   if(db->keyExists("face_refine_op"))
   {
      d_face_refine_op_str = db->getString("face_refine_op");
   }
   else
   {
      TBOX_ERROR("PCDensityMultilevelOperator" 
                 << " -- Required key `face_refine_op'"
                 << " missing in input.");
   }

   if(db->keyExists("cell_coarsen_op"))
   {
      d_cell_coarsen_op_str = db->getString("cell_coarsen_op");
   }
   else
   {
      TBOX_ERROR("PCDensityMultilevelOperator" 
                 << " -- Required key `cell_coarsen_op'"
                 << " missing in input.");
   }

   if (db->keyExists("adjust_cf_coefficients")) 
   {
      d_adjust_cf_coefficients = db->getBool("adjust_cf_coefficients");      
   } 
   else 
   {
      TBOX_ERROR( "PCDensityMultilevelOperator" 
                 << " -- Required key `adjust_cf_coefficients'"
                 << " missing in input.");
   }

   if(db->keyExists("use_cf_interpolant"))
   {
      d_use_cf_interpolant = db->getBool("use_cf_interpolant");
      
      if(d_use_cf_interpolant)
      {
         d_cell_refine_op_str = "CONSTANT_REFINE";

         if (db->keyExists("variable_order_interpolation")) 
         {
            d_variable_order_interpolation = db->getBool("variable_order_interpolation");
         } 
         else 
         {
            TBOX_ERROR("PCDensityMultilevelOperator"
                       << " -- Required key `variable_order_interpolation'"
                       << " missing in input.");
         }
      }
      else
      {
         if(db->keyExists("cell_refine_op"))
         {
            d_cell_refine_op_str = db->getString("cell_refine_op");
         }
         else
         {
            TBOX_ERROR( "PCDensityMultilevelOperator"
                       << " -- Required key `cell_refine_op'"
                       << " missing in input.");
         }
      }
   }
   else
   {
      TBOX_ERROR("PCDensityMultilevelOperator"
                 << " -- Required key `use_cf_interpolant'"
                 << " missing in input.");
   }

   if (db->keyExists("boundary_conditions")) 
     {
       const int dim = d_hierarchy->getDim().getValue();
      // get the database object for boundary conditions
      db->getIntegerArray("boundary_conditions", d_bdry_types, 2*dim);
      
      tbox::Pointer<geom::CartesianGridGeometry > grid_geometry = 
         d_hierarchy->getGridGeometry();
      hier::IntVector shift = grid_geometry->getPeriodicShift(hier::IntVector::getOne(d_hierarchy->getDim()));
      
      for (int d=0; d<dim; d++) 
	{
	  //	  tbox::pout << "Boundary conditions in direction " << d << " are " << d_bdry_types[2*d] << ", " << d_bdry_types[2*d+1] << std::endl; 	  
	  if (shift[d]!=0)
	    {
	      d_bdry_types[2*d+0] = PERIODIC;
	      d_bdry_types[2*d+1] = PERIODIC;
         }
	}

   } 
   else 
   {
      TBOX_ERROR("PCDensityMultilevelOperator" 
                 << " -- Required key `boundary_conditions'"
                 << " missing in input.");
   }   

}


void 
PCDensityMultilevelOperator::apply(const int ln,
				   const int *f_id,
				   const int *u_id, 
				   const int *r_id,
				   const int *f_idx,
				   const int *u_idx,
				   const int *r_idx,
				   const double a,
				   const double b)
{
#ifdef DEBUG_CHECK_ASSERTIONS
  assert(f_id!=NULL);
  assert(u_id!=NULL);
  assert(r_id!=NULL);
  assert(!d_level_operators[ln].isNull());
#endif
  
  d_level_operators[ln]->apply(f_id , u_id , r_id ,
			       f_idx, u_idx, r_idx,
			       a, b);
  
}

LevelOperator *
PCDensityMultilevelOperator::getLevelOperator(const int ln)
{
#ifdef DEBUG_CHECK_ASSERTIONS
  assert(!d_level_operators[ln].isNull());
#endif
  
  return d_level_operators[ln].getPointer();
}
  

void
PCDensityMultilevelOperator::apply(const int coarse_ln,
				   const int fine_ln,
				   const int *f_id,
				   const int *u_id, 
				   const int *r_id,
				   const int *f_idx,
				   const int *u_idx,
				   const int *r_idx,
				   const double a,
				   const double b)
{
#ifdef DEBUG_CHECK_ASSERTIONS
  assert(f_id!=NULL);
  assert(u_id!=NULL);
  assert(r_id!=NULL);
#endif
  
  setupTransferSchedules();
  
  // the order of coomunications required
  // is different depending on whether fluxes are
  // coarsened or not.
  if(d_coarsen_diffusive_fluxes)
    {
      for(int ln=fine_ln; ln>=coarse_ln; ln--)
	{
	  PCDensityLevelOperator *pcDensityLevelOp = dynamic_cast<PCDensityLevelOperator *>(this->getLevelOperator(ln));
	  
#ifdef DEBUG_CHECK_ASSERTIONS
	  assert(pcDensityLevelOp!=NULL);
#endif
	  
	  pcDensityLevelOp->setFlux(d_flux_id, u_id, u_idx);
	}
      
      for(int ln=fine_ln-1; ln>=coarse_ln; ln--)
	{
	  bool coarsen_rhs = (b!=0.0)?true:false;         
	  coarsenSolutionAndSourceTerm(ln, u_id[0], f_id[0], coarsen_rhs);      
	  
#ifdef DEBUG_CHECK_ASSERTIONS
	  assert(d_schedules_initialized==true);
         assert(!d_flux_coarsen_schedule[ln].isNull());
#endif
         d_flux_coarsen_schedule[ln]->coarsenData();
      }
      
      for(int ln=fine_ln; ln>=coarse_ln; ln--)
	{
	  PCDensityLevelOperator *pcDensityLevelOp = dynamic_cast<PCDensityLevelOperator *>(this->getLevelOperator(ln));
	  
#ifdef DEBUG_CHECK_ASSERTIONS
	  assert(pcDensityLevelOp!=NULL);
#endif
	  
	  pcDensityLevelOp->apply(d_flux_id,
				  f_id , u_id , r_id, 
				  f_idx, u_idx, r_idx,
				  a, b);
	}
    }
  else
    {
      for(int ln=fine_ln; ln>=coarse_ln; ln--)
	{
	  PCDensityLevelOperator *pcDensityLevelOp = dynamic_cast<PCDensityLevelOperator *>(this->getLevelOperator(ln));
	  
#ifdef DEBUG_CHECK_ASSERTIONS
	  assert(pcDensityLevelOp!=NULL);
#endif
	  
	  if(ln<fine_ln)
	    {
	      bool coarsen_rhs = (b!=0.0)?true:false;         
	      coarsenSolutionAndSourceTerm(ln, u_id[0], f_id[0], coarsen_rhs);      
	    }
	  
	  pcDensityLevelOp->setFlux(d_flux_id, u_id, u_idx);
	  
	  if((ln<fine_ln)&&(d_coarsen_diffusive_fluxes))
	    {
#ifdef DEBUG_CHECK_ASSERTIONS
	      assert(d_schedules_initialized==true);
	      assert(!d_flux_coarsen_schedule[ln].isNull());
#endif
	      d_flux_coarsen_schedule[ln]->coarsenData();
	    } 
	  
	  pcDensityLevelOp->apply(d_flux_id,
				  f_id , u_id , r_id, 
				  f_idx, u_idx, r_idx,
				  a, b);
	  
	  if(ln<fine_ln && (!d_coarsen_diffusive_fluxes))
	    {
	      tbox::pout << "ERROR:: currently diffusive fluxes must be coarsened" << std::endl;
	      abort();
	    }
	}
    }
}
  
  
void
PCDensityMultilevelOperator::setupTransferSchedules(void)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!d_hierarchy.isNull());
   assert(d_flux_id>=0);
#endif

   int ln;
   
   if(!d_schedules_initialized)
   {
      tbox::Pointer<hier::PatchLevel > flevel = d_hierarchy->getPatchLevel(0);
      tbox::Pointer<hier::PatchLevel > clevel;

      tbox::Pointer<hier::GridGeometry > geometry = flevel->getGridGeometry();

      xfer::CoarsenAlgorithm flux_coarsen_alg(d_hierarchy->getDim());

      flux_coarsen_alg.registerCoarsen(d_flux_id, d_flux_id,
                                       geometry->lookupCoarsenOperator(d_flux,d_face_coarsen_op_str));
      
      for (ln = 1; ln < d_hierarchy->getNumberOfLevels(); ln++) 
      {
         flevel = d_hierarchy->getPatchLevel(ln);
         clevel = d_hierarchy->getPatchLevel(ln-1);
         d_flux_coarsen_schedule[ln-1] = flux_coarsen_alg.createSchedule(clevel, flevel); 
      }
      
      d_schedules_initialized=true;

   }
}


void
PCDensityMultilevelOperator::setFlux(const int coarse_ln,
				     const int fine_ln,
				     const int *u_id,
				     const int *u_idx)
{
#ifdef DEBUG_CHECK_ASSERTIONS
      assert(d_flux_id>=0);
      assert(u_id!=NULL);
#endif   

   for(int ln=fine_ln; ln>=coarse_ln; ln--)
   {
	PCDensityLevelOperator *pcDensityLevelOp = dynamic_cast<PCDensityLevelOperator *>(this->getLevelOperator(ln));

#ifdef DEBUG_CHECK_ASSERTIONS
	assert(pcDensityLevelOp!=NULL);
#endif
      
     pcDensityLevelOp->setFlux(d_flux_id, u_id, u_idx);
     
     if(ln<fine_ln && d_coarsen_diffusive_fluxes)
       {
#ifdef DEBUG_CHECK_ASSERTIONS
         assert(d_schedules_initialized==true);
         assert(!d_flux_coarsen_schedule[ln].isNull());
#endif
         d_flux_coarsen_schedule[ln]->coarsenData();
      } 
   }
}

void
PCDensityMultilevelOperator::interpolate(const tbox::Pointer<hier::PatchLevel > &flevel,
					 const tbox::Pointer<hier::PatchLevel > &clevel,
					 const int *dst_id,
					 const int *src_id,
					 const int *scratch_id,
					 std::string &refine_op)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(dst_id!=NULL);
   assert(src_id!=NULL);
   assert(scratch_id!=NULL);
#endif
   hier::VariableDatabase* variable_db = hier::VariableDatabase::getDatabase();
   tbox::Pointer< hier::Variable > dstvar;
   variable_db->mapIndexToVariable(scratch_id[0], dstvar);

   tbox::Pointer< hier::GridGeometry > geometry = d_hierarchy->getGridGeometry();

   xfer::RefineAlgorithm refine_alg(d_hierarchy->getDim());

   refine_alg.registerRefine(dst_id[0], src_id[0], scratch_id[0],
                             geometry->lookupRefineOperator(dstvar, 
                                                            refine_op));
   
   const int ln = flevel->getLevelNumber();

   tbox::Pointer< xfer::RefinePatchStrategy > refine_strategy = d_set_boundary_ghosts;

   PCDensityRefinePatchStrategy *ptr=NULL;

   if(!refine_strategy.isNull())
   {

#ifdef DEBUG_CHECK_ASSERTIONS
      assert(!refine_strategy.isNull());      
#endif
      // set the index to fill for the refine patch strategy
      ptr=dynamic_cast<PCDensityRefinePatchStrategy*>(refine_strategy.getPointer());

#ifdef DEBUG_CHECK_ASSERTIONS
      assert(ptr!=NULL);
#endif

      ptr->setIndexToFill(scratch_id[0], 0);

   }

   if(d_interpolate_schedule[ln].isNull()||(!refine_alg.checkConsistency(d_interpolate_schedule[ln])))
   {
      if(ln==0)
      {
         d_interpolate_schedule[ln]=refine_alg.createSchedule(flevel,
                                                             ptr);
      }
      else
	{
	  tbox::Pointer<hier::PatchLevel> nullLevelPtr;
         d_interpolate_schedule[ln]=refine_alg.createSchedule(flevel,
                                                              nullLevelPtr,
                                                              ln-1,
                                                              d_hierarchy,
                                                              ptr);
      }
   }
   else
   {
      refine_alg.resetSchedule(d_interpolate_schedule[ln]);
   }

#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!d_interpolate_schedule[ln].isNull());
#endif
   d_interpolate_schedule[ln]->fillData(0.0);
}
  
}
}
