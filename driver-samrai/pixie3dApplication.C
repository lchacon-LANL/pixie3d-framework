//
// $Id: pixie3dApplication.C 3036 2007-05-29 17:30:55Z bphilip $
// $Revision: 3036 $
// $Date: 2007-05-29 11:30:55 -0600 (Tue, 29 May 2007) $
//
//
// File:  pixie3dApplication.C
// Copyright:  (c) 2005 The Regents of the University of California
// Description:  Abstract class that defines interfaces that must be 
//               provided by applications.
//


// Maximum ghost cell width
// If it does not equal 1 there will be extra copies during the refinement
#ifndef GHOST
   #define GHOST 2
#endif

// This forces the ghost cell width to be 1 in the z-direction
// If we are using a refinment ratio of 1 is the z-direction set 
// this to 1 otherwise set it to 0.  Note: the refine operator must
// support a stencil width of 0 in the z-direction if this is set
#ifndef IS2D
   #define IS2D 1
#endif


#include "pixie3dApplication.h"

#include "CartesianGridGeometry.h"
#include "CartesianPatchGeometry.h"
#include "Geometry.h"
#include "HierarchyCellDataOpsReal.h"
#include "BoundaryConditionStrategy.h"
#include "RefineOperator.h"

#include "AMRUtilities.h"
//#include "test_utilities.h"
#include "RefinementBoundaryInterpolation.h"
#include "CartesianCellDoubleWeightedAverage.h"
#include "CartesianCellDoubleCubicCoarsen.h"
#include "CartesianCellDoubleInjectionCoarsen.h"
#include "CartesianFaceDoubleCubicCoarsen.h"
#include "CartesianCellDoubleCubicRefine.h"
#include "CoarsenAlgorithm.h"
#include "CartesianGridGeometry.h"
#include "pixie3dRefinePatchStrategy.h"

#include "fortran.h"
#include "math.h"
#include <iostream>
#include <sstream>
#include <string>
#include "CCellVariable.h"
#include "varrayContainer.h"
#include "CartesianCellDoubleCubicRefine.h"
#include "CartesianCellDoubleLinearRefine.h"

extern "C" {
#include "assert.h"
}

extern "C"{
#include <assert.h>

}



/************************************************************************
*                                                                       *
* External declarations for FORTRAN 77 routines                         *
*                                                                       *
************************************************************************/
/*


#if (NDIM == 1)
#include "fortran/1d/prototypes.h"
#elif (NDIM == 2)
#include "fortran/2d/prototypes.h"
#endif
}*/
extern "C" {
#ifdef absoft
extern void FORTRAN_NAME(GETBE)(void*, int*, int**);
extern void FORTRAN_NAME(FORTRANDESTROY) ();
extern void FORTRAN_NAME(READINPUTFILE) (input_CTX*);
extern void FORTRAN_NAME(FORMEQUILIBRIUM) (void*);
extern void FORTRAN_NAME(INITIALIZE_U_N) (void*);
extern void FORTRAN_NAME(FORMINITIALCONDITION) (void*, int*, double*);
extern void FORTRAN_NAME(EVALUATENONLINEARRESIDUAL) (void*, int*, double*, void*);
#else
  extern void FORTRAN_NAME(getbc)(void*, int*, int**, int*);
extern void FORTRAN_NAME(fortrandestroy) ();
extern void FORTRAN_NAME(readinputfile) (input_CTX*);
extern void FORTRAN_NAME(formequilibrium) (void*);
extern void FORTRAN_NAME(initialize_u_n) (void*);
extern void FORTRAN_NAME(forminitialcondition) (void*, int*, double*);
extern void FORTRAN_NAME(evaluatenonlinearresidual) (void*, int*, double*, void*);
#endif
}

/***********************************************************************
*                                                                      *
* Constructor sets some bad values.                                    *
*                                                                      *
***********************************************************************/
pixie3dApplication::pixie3dApplication()
{
   /*d_number_solution_components = 0;

   d_hierarchy.setNull();
   d_solution.setNull();
   d_application_ctx.setNull();

   d_initial_conditions_id  = -1;
   d_test_problem_source_id = -1;
   d_laplacian_result_id    = -1;

   d_data_server = (TagAndInitDataServer*) NULL;
   d_laplace_operator = (LaplaceOperator*) NULL;
   d_test_problem = (TestProblemStrategy*) NULL;*/
}


/***********************************************************************
*                                                                      *
* Construct from parameter list.  Calls initialize.                    *
*                                                                      *
***********************************************************************/
pixie3dApplication::pixie3dApplication( const pixie3dApplicationParameters* parameters )
{
   //d_coarsen_op_str = "CONSERVATIVE_COARSEN";
   d_coarsen_op_str = "CELL_DOUBLE_INJECTION_COARSEN";

   //d_refine_op_str  = "CONSTANT_REFINE";
   //d_refine_op_str  = "LINEAR_REFINE";   
   d_refine_op_str  = "CELL_DOUBLE_CUBIC_REFINE";

   initialize( parameters );
}


/***********************************************************************
*                                                                      *
* Destructor.                                                          *
*                                                                      *
***********************************************************************/
pixie3dApplication::~pixie3dApplication()
{
   // Delete all fortran variables
   tbox::pout << "Call fortrandestroy\n";
   //#ifdef absoft
   //  FORTRAN_NAME(FORTRANDESTROY) ();
   //#else
   //  FORTRAN_NAME(fortrandestroy) ();
   //#endif
   // Delete the level containers
   LevelContainer *level_container;
   for ( int ln=0; ln<d_hierarchy->getNumberOfLevels(); ln++ )
     {
       level_container = (LevelContainer *) level_container_array[ln];
       delete level_container;
     }
   
   delete [] level_container_array;
}


/***********************************************************************
*                                                                      *
* Initialize object.                                                   *
*                                                                      *
***********************************************************************/
void
pixie3dApplication::initialize( const pixie3dApplicationParameters* parameters )
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert( parameters != (pixie3dApplicationParameters*) NULL );
#endif
   d_object_name = "pixie3d";
   d_hierarchy = parameters->d_hierarchy;

   tbox::pout << "Initializing\n";

   // Get input file name (not needed, it loads pixie3d.in from operating directory)
   /*string input_file;
   if ( parameters->d_db->keyExists("input_file") ) {
      input_file = parameters->d_db->getString("input_file");
   } else {
      TBOX_ERROR( " Key data 'file' missing in input " );
   }*/
   
   // Load input file
   data.npx = 0;
   data.npy = 0;
   data.npz = 0;
   //readInputFile_(&data);
   #ifdef absoft
      FORTRAN_NAME(READINPUTFILE) (&data);
   #else
      FORTRAN_NAME(readinputfile) (&data);
   #endif
   time = 0.0;
   
   assert(data.nvar>0);
   assert(data.nauxs>=0);
   assert(data.nauxv>=0);

   u0_id       = new int[data.nvar];
   u_id        = new int[data.nvar];
   u_tmp_id    = new int[data.nvar];
   
   auxs_id     = new int[data.nauxs];
   auxv_id     = new int[data.nauxv];
   
   auxs_tmp_id = new int[data.nauxs];
   auxv_tmp_id = new int[data.nauxv];

   f_id = new int[data.nvar];
   
   // Check for consistency between domain sizes
   tbox::Pointer<geom::CartesianGridGeometry<NDIM> > grid_geometry = d_hierarchy->getGridGeometry();

   const SAMRAI::hier::BoxArray<NDIM> &physicalDomain = grid_geometry->getPhysicalDomain() ;

   assert(grid_geometry->getDomainIsSingleBox());
   
   int nbox[NDIM];
   
   for (int i=0; i<NDIM; i++)
     {
       nbox[i] = physicalDomain[0].numberCells(i);
     }
   
   if ( data.nxd != nbox[0] )
     {
       TBOX_ERROR("The domain size does not match pixie3d in x-direction");
     }
   
   if ( data.nyd != nbox[1] )
     {
       TBOX_ERROR("The domain size does not match pixie3d in y-direction");
     }
   
   if ( data.nzd != nbox[2] )
     {
       TBOX_ERROR("The domain size does not match pixie3d in z-direction");
     }
   
   // Check for consistency between periodic dimensions
   
   // Check processors?

   // Allocate data for u, u_0
   hier::IntVector<NDIM> ghost0 = hier::IntVector<NDIM>::IntVector(0);
   hier::IntVector<NDIM> ghost1 = hier::IntVector<NDIM>::IntVector(1);
   hier::IntVector<NDIM> ghost2 = hier::IntVector<NDIM>::IntVector(GHOST);
   if ( IS2D == 1 )
      ghost2(2) = 1;
   hier::VariableDatabase<NDIM>* var_db = hier::VariableDatabase<NDIM>::getDatabase();
   tbox::Pointer<hier::VariableContext> context_x = var_db->getContext("pixie3d-x");
   tbox::Pointer<hier::VariableContext> context_f = var_db->getContext("pixie3d-f");
   tbox::Pointer<hier::VariableContext> context_xt = var_db->getContext("pixie3d-x_tmp");
   tbox::Pointer<hier::VariableContext> context_in = var_db->getContext("pixie3d-initial");
   tbox::Pointer<hier::VariableContext> context_eq = var_db->getContext("pixie3d-equilibrium");
   tbox::Pointer< pdat::CellVariable<NDIM,double> > var;
   int var_id;
   std::string var_name;
   d_x = new solv::SAMRAIVectorReal<NDIM,double>("xVec",d_hierarchy,0,d_hierarchy->getFinestLevelNumber());
   d_f = new solv::SAMRAIVectorReal<NDIM,double>("fVec",d_hierarchy,0,d_hierarchy->getFinestLevelNumber());
   d_x_tmp = new solv::SAMRAIVectorReal<NDIM,double>("xTmpVec",d_hierarchy,0,d_hierarchy->getFinestLevelNumber());
   d_initial = new solv::SAMRAIVectorReal<NDIM,double>("xInitialVec",d_hierarchy,0,d_hierarchy->getFinestLevelNumber());
   d_equilibrium = new solv::SAMRAIVectorReal<NDIM,double>("xEquilibriumVec",d_hierarchy,0,d_hierarchy->getFinestLevelNumber());

   std::stringstream stream;
   for (int i=0; i<data.nvar; i++)
     {
       stream << "x(" << i << ")"; 
       var_name = stream.str();
       stream.str("");
       var = new pdat::CellVariable<NDIM,double>( var_name, 1 );
       var_id = var_db->registerVariableAndContext(var, context_x, ghost1);
       d_x->addComponent( var, var_id );
       var_id = var_db->registerVariableAndContext(var, context_f, ghost1);
       d_f->addComponent( var, var_id );
       var_id = var_db->registerVariableAndContext(var, context_xt, ghost2);
       d_x_tmp->addComponent( var, var_id );
       var_id = var_db->registerVariableAndContext(var, context_in, ghost1);
       d_initial->addComponent( var, var_id );
       var_id = var_db->registerVariableAndContext(var, context_eq, ghost1);
      d_equilibrium->addComponent( var, var_id );
     }

   // allocate data for all variables on all levels
   d_x->allocateVectorData();
   d_f->allocateVectorData();
   d_x_tmp->allocateVectorData();
   d_initial->allocateVectorData();
   d_equilibrium->allocateVectorData();

   // set to negative values
   d_x->setToScalar( -1.0 );
   d_f->setToScalar( -1.0 );
   d_initial->setToScalar( -1.0 );
   d_equilibrium->setToScalar( -1.0 );

   // Allocate data for auxillary variables
   tbox::Pointer<hier::VariableContext> context_aux = var_db->getContext("pixie3d-aux");
   d_aux_scalar = new solv::SAMRAIVectorReal<NDIM,double>("auxs",d_hierarchy,0,d_hierarchy->getFinestLevelNumber());
   d_aux_vector = new solv::SAMRAIVectorReal<NDIM,double>("auxv",d_hierarchy,0,d_hierarchy->getFinestLevelNumber());
   d_aux_scalar_tmp = new solv::SAMRAIVectorReal<NDIM,double>("auxs",d_hierarchy,0,d_hierarchy->getFinestLevelNumber());
   d_aux_vector_tmp = new solv::SAMRAIVectorReal<NDIM,double>("auxv",d_hierarchy,0,d_hierarchy->getFinestLevelNumber());
   
   for (int i=0; i<data.nauxs; i++)
     {
       stream << "auxs(" << i << ")"; 
       var_name = stream.str();
       stream.str("");
       var = new pdat::CellVariable<NDIM,double>( var_name, 1 );
       var_id = var_db->registerVariableAndContext(var, context_x, ghost1);
       d_aux_scalar->addComponent( var, var_id );
       var_id = var_db->registerVariableAndContext(var, context_xt, ghost2);
       d_aux_scalar_tmp->addComponent( var, var_id );
     }
   
   for (int i=0; i<data.nauxv; i++)
     {
       stream << "auxv(" << i << ")"; 
       var_name = stream.str();
       stream.str("");
       var = new pdat::CellVariable<NDIM,double>( var_name, NDIM );
       var_id = var_db->registerVariableAndContext(var, context_x, ghost1);
       d_aux_vector->addComponent( var, var_id );
       var_id = var_db->registerVariableAndContext(var, context_xt, ghost2);
       d_aux_vector_tmp->addComponent( var, var_id );
     }
   
   d_aux_scalar->allocateVectorData();
   d_aux_vector->allocateVectorData();
   
   d_aux_scalar->setToScalar( -1.0 );
   d_aux_vector->setToScalar( -1.0 );

   d_aux_scalar_tmp->allocateVectorData();
   d_aux_vector_tmp->allocateVectorData();

   d_aux_scalar_tmp->setToScalar( -1.0 );
   d_aux_vector_tmp->setToScalar( -1.0 );
   
   // Allocate data for f_src
   d_f_src = new pdat::CellVariable<NDIM,double>( "fsrc", data.nvar );
   f_src_id = var_db ->registerVariableAndContext(d_f_src, context_f, ghost0);
   for (int ln=0; ln<d_hierarchy->getNumberOfLevels(); ln++)
     {
       tbox::Pointer<hier::PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);       
       level->allocatePatchData(f_src_id);
     }
   
   // Create patch variables
   int nx, ny, nz;
   for (int i=0; i<data.nvar; i++)
     {
       u0_id[i] = d_initial->getComponentDescriptorIndex(i);
     }
   
   for (int i=0; i<data.nvar; i++)
     {
       u_id[i]  = d_x->getComponentDescriptorIndex(i);
     }
   
   for (int i=0; i<data.nvar; i++)
     {
       u_tmp_id[i]  = d_x_tmp->getComponentDescriptorIndex(i);
     }
   
   for (int i=0; i<data.nauxs; i++)
     {
       auxs_id[i] = d_aux_scalar->getComponentDescriptorIndex(i);
     }
   
   for (int i=0; i<data.nauxv; i++)
     {
       auxv_id[i] = d_aux_vector->getComponentDescriptorIndex(i);
     }
   
   for (int i=0; i<data.nauxs; i++)
     {
       auxs_tmp_id[i] = d_aux_scalar_tmp->getComponentDescriptorIndex(i);
     }
   
   for (int i=0; i<data.nauxv; i++)
     {
       auxv_tmp_id[i] = d_aux_vector_tmp->getComponentDescriptorIndex(i);
     }
   
   for (int i=0; i<data.nvar; i++)
     {
       f_id[i] = d_f->getComponentDescriptorIndex(i);
     }
   
   LevelContainer *level_container;
   level_container_array = new void *[d_hierarchy->getNumberOfLevels()];
   hier::IntVector<NDIM> gcwc;
   for ( int ln=0; ln<d_hierarchy->getNumberOfLevels(); ln++ )
     {
       tbox::Pointer<hier::PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
       const hier::IntVector<NDIM> ratio = level->getRatio();
       nx = nbox[0]*ratio(0);
       ny = nbox[1]*ratio(1);
       nz = nbox[2]*ratio(2);
       // Create the level container
       level_container_array[ln] = new LevelContainer(level->getNumberOfPatches(),nx,ny,nz);
       level_container = (LevelContainer *) level_container_array[ln];
       // Create each patch object
       for (hier::PatchLevel<NDIM>::Iterator p(level); p; p++) {
         tbox::Pointer< hier::Patch<NDIM> > patch = level->getPatch(p());
         level_container->CreatePatch(p(),
				      patch,
				      data.nvar,
				      u0_id,
				      u_id,
				      data.nauxs, auxs_id,
				      data.nauxv,auxv_id);
       }
     }

   // Setup transfer schedule
   generateTransferSchedules();
   
}


/***********************************************************************
*                                                                      *
* Evaluate initial conditions.                                         *
*                                                                      *
***********************************************************************/
void pixie3dApplication::setInitialConditions( const double initial_time )
{
   d_initial_time = initial_time;

   // Form equilibrium
   // Loop through hierarchy
   for ( int ln=0; ln<d_hierarchy->getNumberOfLevels(); ln++ )
     {
       tbox::Pointer<hier::PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
       // Get the Level container
       LevelContainer *level_container = (LevelContainer *) level_container_array[ln];
       // Loop through the different patches
       for (hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
	 {
	   // Form Equilibium
#ifdef absoft
	   FORTRAN_NAME(FORMEQUILIBRIUM)(level_container->getPtr(p()));
#else
	   FORTRAN_NAME(formequilibrium)(level_container->getPtr(p()));
#endif
	 }
     }

   // Initialize u_n
   // Loop through hierarchy
   for ( int ln=0; ln<d_hierarchy->getNumberOfLevels(); ln++ )
     {
       tbox::Pointer<hier::PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
       // Get the Level container
       LevelContainer *level_container = (LevelContainer *) level_container_array[ln];
       // Loop through the different patches
       for (hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
	 {
	   // Form Equilibium
#ifdef absoft
	   FORTRAN_NAME(INITIALIZE_U_N)(level_container->getPtr(p()));
#else
	   FORTRAN_NAME(initialize_u_n)(level_container->getPtr(p()));
#endif
	 }
     }

   // Get d_BoundaryConditionSequence
   // Loop through hierarchy
   int it=0;
   for ( int ln=0; ln<d_hierarchy->getNumberOfLevels(); ln++ )
     {
       tbox::Pointer<hier::PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
       // Get the Level container
       LevelContainer *level_container = (LevelContainer *) level_container_array[ln];
       // Loop through the different patches
       for (hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
	 {
	   // Get d_BoundaryConditionSequence
#ifdef absoft
	   FORTRAN_NAME(GETBC)(level_container->getPtr(p()),&d_NumberOfBoundaryConditions,&d_BoundaryConditionSequence, &it);
#else
	   FORTRAN_NAME(getbc)(level_container->getPtr(p()),&d_NumberOfBoundaryConditions,&d_BoundaryConditionSequence, &it);
#endif
	 }
     }
   
   for (int i=0; i<d_NumberOfBoundaryConditions; i++)
     {
       tbox::pout << "BC(" << i+1 << ",:) = " << d_BoundaryConditionSequence[i] << ", " << d_BoundaryConditionSequence[i+d_NumberOfBoundaryConditions] << "\n";
     }
   
   // Apply boundary conditions
   synchronizeVariables();
   
   // Form initial conditions
   // Loop through hierarchy
   for ( int ln=0; ln<d_hierarchy->getNumberOfLevels(); ln++ )
     {
       tbox::Pointer<hier::PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
       // Get the Level container
       LevelContainer *level_container = (LevelContainer *) level_container_array[ln];
       // Loop through the different patches
       for (hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
	 {
	   tbox::Pointer< hier::Patch<NDIM> > patch = level->getPatch(p());
	   tbox::Pointer< pdat::CellData<NDIM,double> > tmp = patch->getPatchData(f_src_id);
	   double *fsrc = tmp->getPointer();
	   int n_elem = patch->getBox().size()*tmp->getDepth();
	   // Form Initial Conditions
#ifdef absoft
	   FORTRAN_NAME(FORMINITIALCONDITION)(level_container->getPtr(p()),&n_elem,fsrc);
#else
	   FORTRAN_NAME(forminitialcondition)(level_container->getPtr(p()),&n_elem,fsrc);
#endif
	 }
     }
}


/***********************************************************************
*                                                                      *
* Evaluate initial conditions on a single level                        *
*                                                                      *
***********************************************************************/
void pixie3dApplication::setInitialConditions( const double initial_time, tbox::Pointer< hier::PatchLevel<NDIM> > level )
{
   d_initial_time = initial_time;
}


/***********************************************************************
*                                                                      *
* Register vector for setting initial conditions.                      *
*                                                                      *
***********************************************************************/
void
pixie3dApplication::setInitialConditions( tbox::Pointer< solv::SAMRAIVectorReal<NDIM,double> > ic )
{
}


/***********************************************************************
*                                                                      *
* Provide descriptor indices managed here to be allocated.             *
*                                                                      *
***********************************************************************/
hier::ComponentSelector 
pixie3dApplication::getDataToAllocate()
{
   return( d_application_data );
}


/***********************************************************************
*                                                                      *
* Provide descriptor indices managed here to be timestamped.           *
*                                                                      *
***********************************************************************/
hier::ComponentSelector 
pixie3dApplication::getDataToTimeStamp()
{
   return( d_application_data );
}


/***********************************************************************
*                                                                      *
* Empty implementation.                                                *
*                                                                      *
***********************************************************************/
void pixie3dApplication::setValuesOnNewLevel( tbox::Pointer< hier::PatchLevel<NDIM> > level )
{
}


/***********************************************************************
*                                                                      *
* Provide list of variables used by this test problem.                 *
*                                                                      *
***********************************************************************/
tbox::Array< tbox::Pointer< hier::Variable<NDIM> > >
pixie3dApplication::getVariables()
{
   return(d_variable_list);
}

void
pixie3dApplication::apply(const int *f_id,
			  const int *u_id, 
			  const int *r_id,
			  const int *f_idx=NULL,
			  const int *u_idx=NULL,
			  const int *r_idx=NULL,
			  const double a,
			  const double b )
{
}

/***********************************************************************
 *                                                                      *
 * Evaluate right-hand side of IVP being solved.                        *
 *                                                                      *
 ***********************************************************************/
void
pixie3dApplication::apply( tbox::Pointer< solv::SAMRAIVectorReal<NDIM,double> >  f,
			   tbox::Pointer< solv::SAMRAIVectorReal<NDIM,double> >  x,
			   tbox::Pointer< solv::SAMRAIVectorReal<NDIM,double> >  r,
			   double a, double b)
{
  // Copy x
  d_x->copyVector(x);
  
  // Coarsen and Refine x
  
  // Apply boundary conditions
  synchronizeVariables();
  
  // Call EvaluateFunction
  // Loop through hierarchy
  for ( int ln=0; ln<d_hierarchy->getNumberOfLevels(); ln++ ) 
    {
      tbox::Pointer<hier::PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
      // Get the Level container
      LevelContainer *level_container = (LevelContainer *) level_container_array[ln];
      // Loop through the different patches
      for (hier::PatchLevel<NDIM>::Iterator p(level); p; p++) {
	tbox::Pointer< hier::Patch<NDIM> > patch = level->getPatch(p());
	// Get fsrc
	tbox::Pointer< pdat::CellData<NDIM,double> > tmp = patch->getPatchData(f_src_id);
	double *fsrc = tmp->getPointer();
	int n_elem = patch->getBox().size()*tmp->getDepth();
	// Create varray on fortran side
	varrayContainer *varray = new varrayContainer(patch,data.nvar,f_id);
	// Create f
#ifdef absoft
	FORTRAN_NAME(EVALUATENONLINEARRESIDUAL)(level_container->getPtr(p()),&n_elem,fsrc,varray->getPtr());
#else
	FORTRAN_NAME(evaluatenonlinearresidual)(level_container->getPtr(p()),&n_elem,fsrc,varray->getPtr());
#endif
	// Delete varray
	delete varray;
      }
    }
  
  // Copy r
  r->copyVector(d_f);
  r->scale(-1.0, r);
}

/***********************************************************************
 *                                                                      *
 * Print identifying string.                                            *
 *                                                                      *
 ***********************************************************************/
void
pixie3dApplication::printObjectName( std::ostream& os ) 
{
  os << d_object_name;
}

/***********************************************************************
* Print the data for variables on a patch hierarchy. This function is  *
* just to ensure that I am accessing the data correctly.               *
***********************************************************************/
void pixie3dApplication::printVector( const tbox::Pointer< solv::SAMRAIVectorReal<NDIM,double> > vector)
{
   tbox::pout << vector->getName() << "\n";
   const tbox::Pointer< hier::PatchHierarchy<NDIM> > hierarchy = vector->getPatchHierarchy();
   // Find the depth of the vector
   tbox::Pointer<hier::PatchLevel<NDIM> > level = hierarchy->getPatchLevel(0);
   tbox::Pointer<hier::Patch<NDIM> > patch = level->getPatch(0);
   tbox::Pointer< pdat::CellData<NDIM,double> > tmp = vector->getComponentPatchData(0,*patch);
   // Loop through the different levels
   for (int ln = 0; ln < hierarchy->getNumberOfLevels(); ln++) {
      tbox::pout << "Level # = " << ln << "\n";
      tbox::Pointer<hier::PatchLevel<NDIM> > level = hierarchy->getPatchLevel(ln);
      // Loop through the different patched (this is currently only 1 processor!!)
      for (hier::PatchLevel<NDIM>::Iterator p(level); p; p++) {
         tbox::Pointer<hier::Patch<NDIM> > patch = level->getPatch(p());
         const hier::Index<NDIM> ifirst = patch->getBox().lower();
         const hier::Index<NDIM> ilast  = patch->getBox().upper();
         const tbox::Pointer<geom::CartesianPatchGeometry<NDIM> > patch_geometry = patch->getPatchGeometry();
         const double* dx = patch_geometry->getDx();
         const double* xlo = patch_geometry->getXLower();
         // Get a pointer to the data
         tbox::Pointer< pdat::CellData<NDIM,double> > tmp = vector->getComponentPatchData(0,*patch);
         const int depth = tmp->getDepth();
         const hier::IntVector<NDIM> ghost = tmp->getGhostCellWidth();
         double *data = tmp->getPointer();
#if NDIM == 1
	 // 1d version
	 tbox::pout << "  Patch # = " << p << ", ifirst = " << ifirst << ", ilast = " << ilast 
		    << ", ghost = " << ghost << ", xlo = " << xlo[0] << ", dx = " << dx[0] << "\n";
	 int num1 = ilast(0)-ifirst(0)+1+2*ghost(0);
	 for (int k=0; k<depth; k++) {
	   tbox::pout << "    (1) ";
	   for (int i=0; i<num1; i++)
	     tbox::pout << scientific << setprecision(4) << setw(12) << data[i+k*num1] << " ";
	   tbox::pout << "\n";
	 }
#elif NDIM == 2
	 // Put 2d version here
	 tbox::pout << "  Patch # = " << p << ", ifirst = " << ifirst << ", ilast = " << ilast 
		    << ", ghost = " << ghost << ", xlo = (" << xlo[0] << "," << xlo[1] << ")" 
		    << ", dx = (" << dx[0] << "," << dx[1] << ")\n";
	 int num1 = ilast(0)-ifirst(0)+1+2*ghost(0);
	 int num2 = ilast(1)-ifirst(1)+1+2*ghost(1);
	 for (int k=0; k<depth; k++)
	   {
	     for (int j=0; j<num2; j++)
	       {
		 if ( j==0 )
		   tbox::pout << "    (" << k << ") ";
		 else
		   tbox::pout << "        ";
		 for (int i=0; i<num1; i++)
		      tbox::pout << scientific << setprecision(4) << setw(12) << data[i+j*num1+k*num1*num2] << " ";
		 tbox::pout << "\n";
	       }
	   }
#else
	 // Throw error
	 TBOX_ERROR( "printVector not yet programmed for 3d" );
#endif
      }
   }
   tbox::pout << "\n";
}

/***********************************************************************
* Coarsen the data                                                     *
***********************************************************************/
void pixie3dApplication::coarsenVariables(void)
{
   for ( int ln=d_hierarchy->getFinestLevelNumber()-1; ln>=0; ln-- ) 
   {
      d_cell_coarsen_schedules[ln]->coarsenData();
   }

}


/***********************************************************************
* Refine the data                                                      *
***********************************************************************/
void  pixie3dApplication::refineVariables(void)
{
   tbox::Pointer< geom::CartesianGridGeometry <NDIM> > grid_geometry = d_hierarchy->getGridGeometry();
   tbox::Pointer< hier::Variable< NDIM > > var;
   
   // Copy the data from the single ghost cell width variables to the multiple ghost cell width variables
   if ( GHOST != 1 )
     {
       // copy the vectors componentwise
       d_x_tmp->copyVector(d_x, false);
       d_aux_scalar_tmp->copyVector(d_aux_scalar, false);
       d_aux_vector_tmp->copyVector(d_aux_vector, false);
     }
   
   // Fill ghost cells
   // moving from coarser to finer levels fill boundary conditions
   for ( int ln=0; ln<d_hierarchy->getNumberOfLevels(); ln++ )
     {
       for( int i=0; i<d_NumberOfBoundaryConditions; i++)
	 {
	   ((pixie3dRefinePatchStrategy*) d_refine_strategy)->setRefineStrategyDataId(i);            
	   d_refine_schedules[ln]->fillData(0.0);
	   // call applybc for interior patches
	   tbox::Pointer<hier::PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
	   for (hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
	     {
	       tbox::Pointer< hier::Patch<NDIM> > patch = level->getPatch(p());
	       bool touches_boundary = ((pixie3dRefinePatchStrategy*) d_refine_strategy)->checkPhysicalBoundary(*patch);
	       if ( !touches_boundary )
		 {
		   ((pixie3dRefinePatchStrategy*) d_refine_strategy)->setPhysicalBoundaryConditions(*patch,0.0,hier::IntVector<NDIM>::IntVector(1));
		 }
	     }
	 }
       d_level_schedules[ln]->fillData(0.0);
     }

   // Copy the data from the multiple ghost cell width variables to the single ghost cell width variables
   if ( GHOST != 1 )
     {
       // copy the vectors componentwise
       d_x->copyVector(d_x_tmp, false);
       d_aux_scalar->copyVector(d_aux_scalar_tmp, false);
       d_aux_vector->copyVector(d_aux_vector_tmp, false);
     }  
}


/***********************************************************************
* Apply the coarse and refine operators                                *
***********************************************************************/
void pixie3dApplication::synchronizeVariables(void)
{
   coarsenVariables();
   refineVariables();
}


/***********************************************************************
* Create the refinement and coarsen schedules                          *
***********************************************************************/
void
pixie3dApplication::generateTransferSchedules(void)
{

   // Add refinement operator
   tbox::Pointer< geom::CartesianGridGeometry <NDIM> > grid_geometry = d_hierarchy->getGridGeometry();
   tbox::Pointer< xfer::RefineOperator<NDIM> > refine_op;
   if(d_refine_op_str=="CELL_DOUBLE_CUBIC_REFINE")
     {
       // Create the CartesianCellDoubleCubicRefine operator
       CartesianCellDoubleCubicRefine* temp = new CartesianCellDoubleCubicRefine();
       // manually set the refinement ratio and stencil width for the CartesianCellDoubleCubicRefine
      hier::IntVector<NDIM> width = hier::IntVector<NDIM>::IntVector(2);
      hier::IntVector<NDIM> ratio = hier::IntVector<NDIM>::IntVector(3);
      if ( IS2D == 1 )
	{
	  width(2) = 0;
	  ratio(2) = 1;
	}
      
      temp->setStencilWidth(width);
      temp->setRefinementRatio(ratio);
      // Add the refinement operator
      refine_op = temp;
      grid_geometry->addSpatialRefineOperator ( refine_op ) ;
     } 
   
   // Add coarsen operator
   tbox::Pointer< xfer::CoarsenOperator<NDIM> > coarsen_op;
   if(d_coarsen_op_str=="CELL_DOUBLE_INJECTION_COARSEN")
     {
       coarsen_op = new CartesianCellDoubleInjectionCoarsen();
       grid_geometry->addSpatialCoarsenOperator ( coarsen_op ) ;
     }
   else if(d_coarsen_op_str=="CELL_DOUBLE_CUBIC_COARSEN")
     {
       coarsen_op = new CartesianCellDoubleCubicCoarsen();
       grid_geometry->addSpatialCoarsenOperator ( coarsen_op ) ;      
     } 

   // Setup pixie3dRefinePatchStrategy
   d_refine_strategy = (xfer::RefinePatchStrategy<NDIM>*)new pixie3dRefinePatchStrategy();
   ((pixie3dRefinePatchStrategy*) d_refine_strategy)->setGridGeometry(grid_geometry);
   ((pixie3dRefinePatchStrategy*) d_refine_strategy)->setPixie3dHierarchyData(level_container_array);
   if ( GHOST == 1 )
     {
     ((pixie3dRefinePatchStrategy*) d_refine_strategy)->setPixie3dDataIDs(0, data.nvar, data.nauxs, data.nauxv,
									  u_id, u_tmp_id, auxs_id, auxs_tmp_id, auxv_id, auxv_tmp_id );
     }
   else
     {
       ((pixie3dRefinePatchStrategy*) d_refine_strategy)->setPixie3dDataIDs(1, data.nvar, data.nauxs, data.nauxv,
									    u_id, u_tmp_id, auxs_id, auxs_tmp_id, auxv_id, auxv_tmp_id );
     }

   // do registerRefine for a scalar refine algorithm
   tbox::Pointer< hier::Variable<NDIM> > var0;
   int data_id;
   
   // u_n
   if ( data.nvar > 0 )
     {
       for(int i=0; i<data.nvar; i++)
	 {
	   if ( GHOST == 1 )
	     {
	       var0 = d_x->getComponentVariable(i);
	       data_id = d_x->getComponentDescriptorIndex(i);
	     }
	   else
	     {
	       var0 = d_x_tmp->getComponentVariable(i);
	       data_id = d_x_tmp->getComponentDescriptorIndex(i);
	     }
	   d_refine_algorithm.registerRefine( data_id, data_id, data_id,
					 grid_geometry->lookupRefineOperator(var0,d_refine_op_str) );
	 }
     }

   // auxilliary scalars
   if ( data.nauxs > 0 )
     {
       for(int i=0; i<data.nauxs; i++)
	 {
	   if ( GHOST == 1 )
	     {
	       var0 = d_aux_scalar->getComponentVariable(i);
	       data_id = d_aux_scalar->getComponentDescriptorIndex(i);
	     }
	   else
	     {
	       var0 = d_aux_scalar_tmp->getComponentVariable(i);
	       data_id = d_aux_scalar_tmp->getComponentDescriptorIndex(i);
	     }
	   
	   d_refine_algorithm.registerRefine( data_id, data_id, data_id,
					      grid_geometry->lookupRefineOperator(var0,d_refine_op_str) );
	 }
       
     }
   
   // auxilliary vectors
   if ( data.nauxv > 0 )
     {
       var0 = d_aux_vector->getComponentVariable(0);
       for (int i=0; i<data.nauxv; i++)
	 {
	   if ( GHOST == 1 )
	     {
	       var0 = d_aux_vector->getComponentVariable(i);
	       data_id = d_aux_vector->getComponentDescriptorIndex(i);
	     }
	   else
	     {
	       var0 = d_aux_vector_tmp->getComponentVariable(i);
	       data_id = d_aux_vector_tmp->getComponentDescriptorIndex(i);
	     }
	   d_refine_algorithm.registerRefine( data_id, data_id, data_id,
					      grid_geometry->lookupRefineOperator(var0,d_refine_op_str) );
	 }
     }

   const int hierarchy_size = d_hierarchy->getNumberOfLevels();
   d_refine_schedules.resizeArray(hierarchy_size);
   d_level_schedules.resizeArray(hierarchy_size);
   d_cell_coarsen_schedules.resizeArray(hierarchy_size);
   
   for (int ln=0; ln<hierarchy_size; ln++)
     {
       tbox::Pointer<hier::PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
       d_refine_schedules[ln] = d_refine_algorithm.createSchedule( level, ln-1, 
								   d_hierarchy, (xfer::RefinePatchStrategy<NDIM>*)d_refine_strategy );
       d_level_schedules[ln] = d_refine_algorithm.createSchedule( 
								 level, (xfer::RefinePatchStrategy<NDIM>*) NULL );
       if(ln<hierarchy_size-1)
	 {
	   tbox::Pointer<hier::PatchLevel<NDIM> > flevel = d_hierarchy->getPatchLevel(ln+1);         
	   d_cell_coarsen_schedules[ln]=d_cell_coarsen_alg.createSchedule(level,flevel);
	 }
     }   
}

int
pixie3dApplication::getNumberOfDependentVariables()
{
  //  assert(data!=NULL);
  return data.nvar;
}
  
