#include "SAMRAI/geom/CartesianGridGeometry.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/xfer/RefineAlgorithm.h"
#include "SAMRAI/pdat/FaceData.h"
#include "SAMRAI/pdat/CellVariable.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/geom/CartesianPatchGeometry.h"

#include "source/AMRUtilities.h"

#include "PCDiagonalLevelOperator.h"

namespace SAMRAI {

namespace SAMRSolvers {

PCDiagonalLevelOperator::PCDiagonalLevelOperator()
{
}

PCDiagonalLevelOperator::PCDiagonalLevelOperator(LevelOperatorParameters *parameters):LevelOperator(parameters)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(parameters!=NULL);
#endif

   d_interpolate_ghost_values     = false;
   d_sibling_fill_cached          = false;
   d_variable_order_interpolation = false;
   d_coefficients_changed         = true;
   d_reset_ghost_cells            = true;
   
   d_sibling_fill_schedule.setNull();

   const int dim = d_level->getDim().getValue();

   d_bdry_types = new int[2*dim];

   // verify that any routines coming after this one do
   // not overwrite values initialized in getFromInput
   // unintentionally
   PCDiagonalLevelOperator::getFromInput(parameters->d_db);

   initializeInternalVariableData();

}

PCDiagonalLevelOperator::~PCDiagonalLevelOperator()
{

   hier::VariableDatabase* variable_db = hier::VariableDatabase::getDatabase();

   if(d_level->checkAllocated(d_flux_id))
   {
      d_level->deallocatePatchData(d_flux_id);
   }

   if(!d_isPartOfMultilevelOp)
     {
       // the calls to remove should only be after the
       // associated memory has been freed
       // the variable will be removed by the multilevel
       // operator if it exists
       variable_db->removePatchDataIndex(d_flux_id);
     }
   
   d_sibling_fill_schedule.setNull();

   delete [] d_bdry_types;
   d_bdry_types = NULL;
   
}

void
PCDiagonalLevelOperator::reset(DiscreteOperatorParameters *params)
{
  PCDiagonalLevelOperator::getFromInput(params->d_db);
}


void
PCDiagonalLevelOperator::initializeInternalVariableData()
{
   hier::VariableDatabase* variable_db = hier::VariableDatabase::getDatabase();

   const hier::IntVector zero_ghosts(hier::IntVector::getZero(d_level->getDim()));

   const tbox::Pointer< hier::VariableContext > scratch_cxt = variable_db->getContext("SCRATCH");

   std::ostringstream ibuffer;
   ibuffer<<d_object_id;
   std::string object_str=ibuffer.str();
   
   // if this operator is a part of a multilevel operator the multilevel operator will
   // supply storage for the flux. The implicit assumption is that both variables will
   // not need to use the flux storage at the same time
   if(!d_isPartOfMultilevelOp)
     {
       std::string cellFlux("PCDiagonalOperator_InternalFlux");
       cellFlux+=object_str;
       
#ifdef DEBUG_CHECK_ASSERTIONS
       if(d_debug_print_info_level>4)
	 {
	   const int ln = (d_level->inHierarchy())?d_level->getLevelNumber():-1;
	   tbox::pout << "PCDiagonalLevelOperator::level number " << ln << std::endl;
	   tbox::pout << "PCDiagonalLevelOperator::Cell Flux variable name " << cellFlux << std::endl;
	   
	 }
#endif
       
       d_flux = variable_db->getVariable(cellFlux);
       
       if (!d_flux)
	 {
	   d_flux = new pdat::FaceVariable<double>(d_level->getDim(), cellFlux,1);
	 }
       
       d_flux_id = variable_db->registerVariableAndContext(d_flux,
							   scratch_cxt,
							   zero_ghosts);
     }
   
   // allocate storage space
   if(!d_level->checkAllocated(d_flux_id))
   {
      d_level->allocatePatchData(d_flux_id);
   }
}

void
PCDiagonalLevelOperator::getFromInput(const tbox::Pointer<tbox::Database> &db)
{
   LevelOperator::getFromInput(db);

#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!db.isNull());
#endif

   // if this operator is part of a multilevel operator the multilevel operator
   // will supply space for the flux
   if(d_isPartOfMultilevelOp)
     {
       if (db->keyExists("flux_id"))
	 {
	   d_flux_id = db->getInteger("flux_id");
	 }
       else
	 {
	   TBOX_ERROR( "PCDiagonalLevelOperator"
		       << " -- Required key `flux_id'"
		       << " missing in input.");
	 }
     }
   
   if (db->keyExists("tangent_interp_scheme"))
   {
      d_tangent_interp_scheme = SAMRAI::RefinementBoundaryInterpolation::lookupInterpolationScheme(db->getString("tangent_interp_scheme"));
   }
   else
   {
      TBOX_ERROR( "PCDiagonalLevelOperator"
                 << " -- Required key `tangent_interp_scheme'"
                 << " missing in input.");
   }

   if (db->keyExists("normal_interp_scheme"))
   {
      d_normal_interp_scheme = SAMRAI::RefinementBoundaryInterpolation::lookupInterpolationScheme(db->getString("normal_interp_scheme"));
   }
   else
   {
      TBOX_ERROR( "PCDiagonalLevelOperator"
                 << " -- Required key `normal_interp_scheme'"
                 << " missing in input.");
   }

   if (db->keyExists("variable_order_interpolation"))
   {
      d_variable_order_interpolation = db->getBool("variable_order_interpolation");
   }
   else
   {
      TBOX_ERROR("PCDiagonalLevelOperator"
                 << " -- Required key `variable_order_interpolation'"
                 << " missing in input.");
   }
   // this key needs to be explicitly specified though on first appearance
   // it may appear that it is set to true only when adjust_cf_coefficients is false.
   // However, it is needed when ghost cell data
   // is already aligned from the coarser level and we do not
   // want to disturb those values. The situation where this arises
   // is on restricted patches underlying fine patch levels during
   // the AFAC or AFACx solution process
   if (db->keyExists("interpolate_ghost_values"))
   {
      d_interpolate_ghost_values = db->getBool("interpolate_ghost_values");

   }
   else
   {
      TBOX_ERROR("PCDiagonalLevelOperator"
                 << " -- Required key `interpolate_ghost_values'"
                 << " missing in input.");
   }

   if (db->keyExists("extrapolation_order"))
     {
       d_extrapolation_order = db->getInteger("extrapolation_order");
     }
   else
     {
       TBOX_ERROR("PCDiagonalLevelOperator"
		  << " -- Required key `extrapolation_order'"
		  << " missing in input.");
     }
   
   if (db->keyExists("boundary_conditions"))
     {
       // get the database object for boundary conditions
       db->getIntegerArray("boundary_conditions", d_bdry_types, 2*(d_level->getDim().getValue()));
       
       tbox::Pointer<geom::CartesianGridGeometry > grid_geometry = d_level->getGridGeometry();
       hier::IntVector shift = grid_geometry->getPeriodicShift(hier::IntVector::getOne(d_level->getDim()));
       
       for (int d=0; d<d_level->getDim().getValue(); d++)
	 {
	   if (shift[d]!=0)
	     {
	       d_bdry_types[2*d+0] = PERIODIC;
	       d_bdry_types[2*d+1] = PERIODIC;
	     }
	 }
     }
   else
     {
       TBOX_ERROR("PCDiagonalLevelOperator"
		  << " -- Required key `boundary_conditions'"
		  << " missing in input.");
     }
   
}

void
PCDiagonalLevelOperator::setExtrapolationOrder(const int extrapolation_order)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(extrapolation_order>=0);
#endif

   bool reset=(extrapolation_order==d_extrapolation_order)?false:true;

   d_extrapolation_order=extrapolation_order;

   if(reset)
   {
      // if any patch touches a physical boundary and the extrapolation
      // order has changed then the stencils at that boundary will be
      // wrong and will need to be recomputed. For now we recompute
      // the stencil at all points
      for (hier::PatchLevel::Iterator p(d_level); p; p++)
      {
         const tbox::Pointer<hier::Patch > &patch = *p;
#ifdef DEBUG_CHECK_ASSERTIONS
         assert(!patch.isNull());
#endif
	 const hier::Box &mappedBox = patch->getBox();

         if(d_level->patchTouchesRegularBoundary(mappedBox.getId()))
         {
	   //            d_stencil_initialized=false; // reset this to false so that stencils are recomputed
         }
      }
   }
}

void
PCDiagonalLevelOperator::setInterpolationScheme(SAMRAI::RefinementBoundaryInterpolation::InterpolationScheme scheme,
					       const int dir)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(dir==0||dir==1);
#endif

   if(dir==0)
   {
      d_tangent_interp_scheme=scheme;
   }
   else
   {
     d_normal_interp_scheme=scheme;
   }

   // BP: we need to do the equivalent for the commented statement below
   // reset this to false so that stencils are recomputed
   //   d_stencil_initialized=false;
}

std::vector<int>
PCDiagonalLevelOperator::getStencilOffsets(const int i,
					  const int j,
					  const int k)
{
   std::vector<int> offsets;

   return offsets;
}

tbox::Pointer< hier::PatchData >
PCDiagonalLevelOperator::getStencilBlock(const tbox::Pointer<hier::Patch> patch,
					const int i,
					const int j,
					const int k)
{
   tbox::Pointer< hier::PatchData > stencilBlock;
   stencilBlock.setNull();

   return stencilBlock;
}

void
PCDiagonalLevelOperator::applyBoundaryCondition(const int *var_id,
					       const int *var_idx,
					       const int *var_components,
					       const int number_of_variables)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(var_id!=NULL);
   assert(var_components==NULL);
#endif

   // if this operator was created using createOperator
   // the coefficients need to be initialized from a level in
   // the hierarchy before we proceed
   if((!d_copy_schedule.isNull())&&d_coefficients_changed)
   {
      d_copy_schedule->fillData(0.0);
      d_coefficients_changed=false;
   }

   tbox::Pointer< hier::RefineOperator > nullRefineOp;

   if(d_reset_ghost_cells)
   {
      xfer::RefineAlgorithm ghost_cell_fill(d_level->getDim());
      ghost_cell_fill.registerRefine(var_id[0], var_id[0], var_id[0], nullRefineOp);

      PCDiagonalRefinePatchStrategy *ptr=NULL;

      if(!d_set_boundary_ghosts.isNull())
      {
         ptr=dynamic_cast<PCDiagonalRefinePatchStrategy*>(d_set_boundary_ghosts.getPointer());

#ifdef DEBUG_CHECK_ASSERTIONS
         assert(ptr!=NULL);
#endif

         ptr->setIndexToFill(var_id[0]);
      }

      if(!d_sibling_fill_cached)
      {
         d_sibling_fill_schedule = ghost_cell_fill.createSchedule(d_level, ptr);
         d_sibling_fill_cached=true;
      }
      else
      {
#ifdef DEBUG_CHECK_ASSERTIONS
         assert(ghost_cell_fill.checkConsistency(d_sibling_fill_schedule));
#endif
         ghost_cell_fill.resetSchedule(d_sibling_fill_schedule);
      }

      d_sibling_fill_schedule->fillData(0.0);
   }

   const int ln = d_level->inHierarchy()?d_level->getLevelNumber():0;
   // it is necessary to allow the case with d_cf_interpolant==NULL
   // without aborting to handle restricted levels in AFAC for example
   if(d_level->inHierarchy()
         &&(ln>0)
         &&(d_interpolate_ghost_values==true)
         &&(d_cf_interpolant!=NULL))
   {
      SAMRAI::RefinementBoundaryInterpolation::InterpolationScheme normal_scheme=d_normal_interp_scheme;

      bool cachedIsVariableInterpolationOrder=d_cf_interpolant->getVariableOrderInterpolation();
      d_cf_interpolant->setVariableOrderInterpolation(d_variable_order_interpolation);

      const int idx=(var_idx==NULL)?0:var_idx[0];
      d_cf_interpolant->interpolateGhostValues(ln,
            d_tangent_interp_scheme,
            normal_scheme,
            var_id[0],
            idx);

      d_cf_interpolant->setVariableOrderInterpolation(cachedIsVariableInterpolationOrder);
   }

}

void
PCDiagonalLevelOperator::setFlux(const int flux_id,
				const int *u_id,
				const int *u_idx)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(u_id!=NULL);
#endif


   // make sure variables are aligned properly for calculating fluxes
   bool cached_interpolate_ghost_values = d_interpolate_ghost_values;
   // for calculating fluxes we need to ensure data is properly
   // aligned in coarse fine ghost cells
   bool cached_reset_ghost_cells        = d_reset_ghost_cells;
   d_interpolate_ghost_values           = true;
   d_reset_ghost_cells                  = true;

   applyBoundaryCondition(u_id, u_idx );

   d_interpolate_ghost_values           = cached_interpolate_ghost_values;
   d_reset_ghost_cells                  = cached_reset_ghost_cells;

   const hier::IntVector no_ghosts(hier::IntVector::getZero(d_level->getDim()));

   const int dim = d_level->getDim().getValue();

   for (hier::PatchLevel::Iterator p(d_level); p; p++)
   {
      tbox::Pointer<hier::Patch > patch = *p;

      const hier::Box box = patch->getBox();

      tbox::Pointer< pdat::FaceData<double> > b_data;
      tbox::Pointer< pdat::CellData<double> > u_data;
      tbox::Pointer< pdat::FaceData<double> > flux_data;

      u_data = patch->getPatchData(u_id[0]);

      hier::IntVector ghost_cell_width = u_data->getGhostCellWidth();
      const int gcw=ghost_cell_width(0);
#ifdef DEBUG_CHECK_ASSERTIONS
      assert(gcw>=0);
#endif
      flux_data = patch->getPatchData(flux_id);

      const tbox::Pointer<geom::CartesianPatchGeometry >
         geometry = patch->getPatchGeometry();

      const hier::Index ifirst = box.lower();
      const hier::Index ilast = box.upper();

      const int idx=(u_idx==NULL)?0:u_idx[0];
   }
}

void
PCDiagonalLevelOperator::apply(const int *f_id,
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


   setFlux(d_flux_id,
           u_id,
           u_idx );

   apply(d_flux_id,
         f_id, u_id, r_id,
         f_idx, u_idx, r_idx,
         a, b);
}

void
PCDiagonalLevelOperator::apply(const int flux_id,
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

   const int dim = d_level->getDim().getValue();
   
   const int fidx=(f_idx==NULL)?0:f_idx[0];
   const int uidx=(u_idx==NULL)?0:u_idx[0];
   const int ridx=(r_idx==NULL)?0:r_idx[0];
   
   for (hier::PatchLevel::Iterator p(d_level); p; p++)
     {
       tbox::Pointer<hier::Patch > patch = *p;
       const hier::Box interior(patch->getBox());
       const hier::Index ifirst  = interior.lower();
       const hier::Index ilast   = interior.upper();
       
#ifdef DEBUG_CHECK_ASSERTIONS
       assert(patch->checkAllocated(flux_id));
       assert(patch->checkAllocated(r_id[0]));
#endif
       
       tbox::Pointer< pdat::CellData<double> > f_data;
       
       if(f_id[0]>=0)
	 {
	   f_data = patch->getPatchData(f_id[0]);
	   
#ifdef DEBUG_CHECK_ASSERTIONS
	   assert(patch->checkAllocated(f_id[0]));
#endif
	 }
       
       tbox::Pointer< pdat::FaceData<double> > flux_data = patch->getPatchData(flux_id);
       tbox::Pointer< pdat::CellData<double> > r_data = patch->getPatchData(r_id[0]);
       
       
       hier::IntVector ghost_cell_width(d_level->getDim());
       int fgcw;
       
       if(!f_data.isNull())
	 {
	   ghost_cell_width = f_data->getGhostCellWidth();
	   fgcw=ghost_cell_width(0);
	   
#ifdef DEBUG_CHECK_ASSERTIONS
	   assert(fgcw>=0);
#endif
	 }
       
       ghost_cell_width = r_data->getGhostCellWidth();
       const int rgcw=ghost_cell_width(0);
       
#ifdef DEBUG_CHECK_ASSERTIONS
       assert(rgcw>=0);
#endif
       
       const tbox::Pointer<geom::CartesianPatchGeometry >
         patch_geometry = patch->getPatchGeometry();
     }
}
  
void
PCDiagonalLevelOperator::initializeDatabase(tbox::Pointer<tbox::Database> db)
  {
    // initialize any required base class entries
   LevelOperator::initializeDatabase(db);

   if (!db->keyExists("tangent_interp_scheme"))
   {
      db->putString("tangent_interp_scheme", "CONSTANT");
   }

   if (!db->keyExists("normal_interp_scheme"))
   {
      db->putString("normal_interp_scheme", "CONSTANT");
   }

   if (!db->keyExists("adjust_cf_coefficients"))
   {
      db->putBool("adjust_cf_coefficients", false);
   }

   if (!db->keyExists("variable_order_interpolation"))
   {
      db->putBool("variable_order_interpolation", false);
   }

   if (!db->keyExists("interpolate_ghost_values"))
   {
      db->putBool("interpolate_ghost_values", false);

   }

   if (!db->keyExists("extrapolation_order"))
   {
      db->putInteger("extrapolation_order", d_extrapolation_order);
   }

   if (!db->keyExists("boundary_conditions"))
   {
      const int dim = d_level->getDim().getValue();
      // put the database object for boundary conditions
      db->putIntegerArray("boundary_conditions", d_bdry_types, 2*dim);
   }
}

LevelOperator *
PCDiagonalLevelOperator::constructOperator(LevelOperatorParameters *parameters)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!parameters->d_level.isNull());
#endif

   parameters->d_cf_interpolant = NULL;
   parameters->d_set_boundary_ghosts = d_set_boundary_ghosts;

   // populate the database with entries that are required for initialization
   // of a PCDiagonalLevelOperator, inserting only entries that do not already
   // exist
   initializeDatabase(parameters->d_db);

   tbox::Pointer<hier::PatchLevel > level = parameters->d_level;

   xfer::RefineAlgorithm level_copy_alg(d_level->getDim());

   tbox::Pointer< hier::RefineOperator > nullRefineOp;

#if 0
   // initialize d_a_index, d_b_index on the level by copies if necessary
   if(d_b_index>=0)
   {
      level->allocatePatchData(d_b_index, 0.0);
      level_copy_alg.registerRefine(d_b_index, d_b_index, d_b_index, nullRefineOp);
   }

   if(d_a_index>=0)
   {
      level->allocatePatchData(d_a_index, 0.0);
      level_copy_alg.registerRefine(d_a_index, d_a_index, d_a_index, nullRefineOp);
   }

   // we create a copy schedule to use whenever the coefficients need to be updated
   // by the created operator
   if((d_a_index>=0)||(d_b_index>=0))
   {
      parameters->d_copy_schedule = level_copy_alg.createSchedule(level, d_level, NULL);
   }
#endif
   
   // create the new operator and cache a pointer to it
   d_constructed_op = new PCDiagonalLevelOperator(parameters);

   return d_constructed_op;
}

int
PCDiagonalLevelOperator::getVariableIndex(std::string &name,
					 tbox::Pointer<hier::VariableContext> &context,
					 tbox::Pointer<hier::Variable > &var,
					 hier::IntVector nghosts,
					 int depth,
					 bool bOverride,
					 std::string centering)
{
   int var_id = -1;

   if(!bOverride)
     {
       hier::VariableDatabase* variable_db = hier::VariableDatabase::getDatabase();
       
       var = variable_db->getVariable(name);
       
       if(!var)
	 {
	   var = new pdat::CellVariable< double>(d_level->getDim(), name, depth);
	 }
       
       var_id = variable_db->registerVariableAndContext(var,
						       context,
							nghosts);
     }
   else
     {
       abort();
     }
   
   return var_id;
}
  
}
}
