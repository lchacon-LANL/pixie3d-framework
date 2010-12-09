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
//   #define GHOST 2
   #define GHOST 1
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
#include <string>
#include "CCellVariable.h"
#include "varrayContainer.h"
#include "CartesianCellDoubleCubicRefine.h"
#include "CartesianCellDoubleLinearRefine.h"
#include "SiblingGhostAlgorithm.h"

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
extern void FORTRAN_NAME(GETBC)(void*, int*, int**);
extern void FORTRAN_NAME(READINPUTFILE) (input_CTX*);
extern void FORTRAN_NAME(CREATEGRIDSTRUCTURES) (void*);
extern void FORTRAN_NAME(FORMEQUILIBRIUM) (void*);
extern void FORTRAN_NAME(INITIALIZE_U_N) (void*);
extern void FORTRAN_NAME(FORMINITIALCONDITION) (void*, int*, double*);
extern void FORTRAN_NAME(EVALUATENONLINEARRESIDUAL) (void*, int*, double*, void*);
#else
extern void FORTRAN_NAME(readinputfile) (input_CTX*);
extern void FORTRAN_NAME(creategridstructures) (void*);
extern void FORTRAN_NAME(formequilibrium) (void*);
extern void FORTRAN_NAME(initialize_u_n) (void*);
extern void FORTRAN_NAME(forminitialcondition) (void*, int*, double*);
extern void FORTRAN_NAME(evaluatenonlinearresidual) (void*, int*, double*, void*);
extern void FORTRAN_NAME(setupvarinitseq)(void *, int*);
extern void FORTRAN_NAME(getnumberofbcgroups)(void*, int*);
extern void FORTRAN_NAME(getbc)(void*, int*, int*, int**);
extern void FORTRAN_NAME(initializeauxvar)(void *, int*);
extern void FORTRAN_NAME(get_var_names)(void *, char*, char*, char*);
#endif
}

namespace SAMRAI{
/***********************************************************************
*                                                                      *
* Constructor sets some bad values.                                    *
*                                                                      *
***********************************************************************/
pixie3dApplication::pixie3dApplication()
{
    TBOX_ERROR("The empty constructor is not supported");
}


/***********************************************************************
*                                                                      *
* Construct from parameter list.  Calls initialize.                    *
*                                                                      *
***********************************************************************/
pixie3dApplication::pixie3dApplication(  pixie3dApplicationParameters* parameters ): SAMRSolvers::DiscreteOperator(parameters)
{
   //d_coarsen_op_str = "CONSERVATIVE_COARSEN";
   d_coarsen_op_str = "CELL_DOUBLE_INJECTION_COARSEN";

   //d_refine_op_str  = "CONSTANT_REFINE";
   //   d_refine_op_str  = "LINEAR_REFINE";   
   d_refine_op_str  = "CELL_DOUBLE_CUBIC_REFINE";

   d_RefineSchedulesGenerated=false;

   d_VizWriter = parameters->d_VizWriter;
   
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
    tbox::pout << "Deleting application\n";
    // Delete the level containers
    LevelContainer *level_container;
    for ( int ln=0; ln<d_hierarchy->getNumberOfLevels(); ln++ ) {
        level_container = (LevelContainer *) level_container_array[ln];
        delete level_container;
    }
    delete [] level_container_array;
    // Delete misc variables
    delete input_data;
    delete [] u0_id;
    delete [] u_id;
    delete [] u_tmp_id;
    delete [] auxs_id;
    delete [] auxv_id;
    delete [] auxs_tmp_id;
    delete [] auxv_tmp_id;
    delete [] f_id;
    delete d_refine_strategy;
    d_initial->freeVectorComponents();
    d_x_tmp->freeVectorComponents();
    d_x->freeVectorComponents();
    d_x_r->freeVectorComponents();
    d_x_ic->freeVectorComponents();
    d_aux_scalar->freeVectorComponents();
    d_aux_vector->freeVectorComponents();;
    d_aux_scalar_tmp->freeVectorComponents();
    d_aux_vector_tmp->freeVectorComponents();
    delete [] depVarLabels;
    delete [] auxScalarLabels;
    delete [] auxVectorLabels;
}


/***********************************************************************
*                                                                      *
* Initialize object.                                                   *
*                                                                      *
***********************************************************************/
void
pixie3dApplication::initialize( pixie3dApplicationParameters* parameters )
{
    #ifdef DEBUG_CHECK_ASSERTIONS
        assert( parameters != (pixie3dApplicationParameters*) NULL );
    #endif
    d_object_name = "pixie3d";
    d_hierarchy = parameters->d_hierarchy;
   
    tbox::pout << "Initializing\n";
    input_data = new input_CTX;
   
    // Load input file
    input_data->npx = 0;
    input_data->npy = 0;
    input_data->npz = 0;
    //readInputFile_(&input_data);
    #ifdef absoft
        FORTRAN_NAME(READINPUTFILE) (input_data);
    #else
        FORTRAN_NAME(readinputfile) (input_data);
    #endif
   
    assert(input_data->nvar>0);
    assert(input_data->nauxs>=0);
    assert(input_data->nauxv>=0);

    u0_id       = new int[input_data->nvar];
    u_id        = new int[input_data->nvar];
    u_tmp_id    = new int[input_data->nvar];
   
    auxs_id     = new int[input_data->nauxs];
    auxv_id     = new int[input_data->nauxv];
   
    auxs_tmp_id = new int[input_data->nauxs];
    auxv_tmp_id = new int[input_data->nauxv];

    f_id = new int[input_data->nvar];
   
    // Check for consistency between domain sizes
    tbox::Pointer<geom::CartesianGridGeometry<NDIM> > grid_geometry = d_hierarchy->getGridGeometry();
    const SAMRAI::hier::BoxArray<NDIM> &physicalDomain = grid_geometry->getPhysicalDomain() ;

    assert(grid_geometry->getDomainIsSingleBox());
   
    int nbox[NDIM];
    for (int i=0; i<NDIM; i++)
        nbox[i] = physicalDomain[0].numberCells(i);
   
    if ( input_data->nxd != nbox[0] )
        TBOX_ERROR("The domain size does not match pixie3d in x-direction");
    if ( input_data->nyd != nbox[1] )
        TBOX_ERROR("The domain size does not match pixie3d in y-direction");
    if ( input_data->nzd != nbox[2] )
        TBOX_ERROR("The domain size does not match pixie3d in z-direction");
   
    // Check for consistency between periodic dimensions
   
    // Check processors?

    // first allocate a weight or control volume variable

    /*
     * hier::Variable<NDIM> to weight solution vector entries on a composite grid.
     */
    
    hier::VariableDatabase<NDIM>* var_db = hier::VariableDatabase<NDIM>::getDatabase();
    tbox::Pointer<hier::VariableContext> context_x = var_db->getContext("pixie3d-x");
    tbox::Pointer<hier::VariableContext> context_xr = var_db->getContext("pixie3d-x_r");
    tbox::Pointer<hier::VariableContext> context_ic = var_db->getContext("pixie3d-x_ic");
    tbox::Pointer<hier::VariableContext> context_xt = var_db->getContext("pixie3d-x_tmp");
    tbox::Pointer<hier::VariableContext> context_in = var_db->getContext("pixie3d-initial");
    tbox::Pointer<hier::VariableContext> context_f = var_db->getContext("pixie3d-source");
    tbox::Pointer< pdat::CellVariable<NDIM,double> > var;
    hier::IntVector<NDIM> ghost0 = hier::IntVector<NDIM>::IntVector(0);
    hier::IntVector<NDIM> ghost1 = hier::IntVector<NDIM>::IntVector(1);
    hier::IntVector<NDIM> ghost2 = hier::IntVector<NDIM>::IntVector(GHOST);
   
    var = new pdat::CellVariable<NDIM,double>("weight", 1);
    const int weight_id = var_db->registerVariableAndContext(var, context_xt, ghost0);

    for (int ln=0; ln<d_hierarchy->getNumberOfLevels(); ln++) {
        tbox::Pointer<hier::PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);       
        level->allocatePatchData(weight_id);
    }

    AMRUtilities::setVectorWeights(d_hierarchy, weight_id);

    // Allocate data for u, u_0, u_ic
    if ( IS2D == 1 )
        ghost2(2) = 1;
    int var_id;
    std::string var_name;
    d_x = new solv::SAMRAIVectorReal<NDIM,double>("xVec",d_hierarchy,0,d_hierarchy->getFinestLevelNumber());
    d_x_r = new solv::SAMRAIVectorReal<NDIM,double>("xVec_r",d_hierarchy,0,d_hierarchy->getFinestLevelNumber());
    d_x_ic = new solv::SAMRAIVectorReal<NDIM,double>("xVecIC",d_hierarchy,0,d_hierarchy->getFinestLevelNumber());
    d_x_tmp = new solv::SAMRAIVectorReal<NDIM,double>("xTmpVec",d_hierarchy,0,d_hierarchy->getFinestLevelNumber());
    d_initial = new solv::SAMRAIVectorReal<NDIM,double>("xInitialVec",d_hierarchy,0,d_hierarchy->getFinestLevelNumber());
    std::stringstream stream;
    for (int i=0; i<input_data->nvar; i++) {
       stream << "x(" << i << ")"; 
       var_name = stream.str();
       stream.str("");
       var = new pdat::CellVariable<NDIM,double>( var_name, 1 );
       var_id = var_db->registerVariableAndContext(var, context_x, ghost1);
       d_x->addComponent( var, var_id, weight_id );
       
       var_id = var_db->registerVariableAndContext(var, context_xt, ghost2);
       d_x_tmp->addComponent( var, var_id, weight_id );
       var_id = var_db->registerVariableAndContext(var, context_in, ghost1);
       d_initial->addComponent( var, var_id, weight_id );
       var_id = var_db->registerVariableAndContext(var, context_ic, ghost0);
       d_x_ic->addComponent( var, var_id, weight_id );
       var_id = var_db->registerVariableAndContext(var, context_xr, ghost1);
       d_x_r->addComponent( var, var_id, weight_id );
    }
    // allocate data for all variables on all levels
    d_x->allocateVectorData();
    d_x_r->allocateVectorData();
    d_x_ic->allocateVectorData();
    d_x_tmp->allocateVectorData();
    d_initial->allocateVectorData();
    
    // set to negative values
    //   d_x->setToScalar( -1.0, false );
    //   d_initial->setToScalar( -1.0, false );

    // Allocate data for auxillary variables
    tbox::Pointer<hier::VariableContext> context_aux = var_db->getContext("pixie3d-aux");
    d_aux_scalar = new solv::SAMRAIVectorReal<NDIM,double>("auxs",d_hierarchy,0,d_hierarchy->getFinestLevelNumber());
    d_aux_vector = new solv::SAMRAIVectorReal<NDIM,double>("auxv",d_hierarchy,0,d_hierarchy->getFinestLevelNumber());
    d_aux_scalar_tmp = new solv::SAMRAIVectorReal<NDIM,double>("auxs",d_hierarchy,0,d_hierarchy->getFinestLevelNumber());
    d_aux_vector_tmp = new solv::SAMRAIVectorReal<NDIM,double>("auxv",d_hierarchy,0,d_hierarchy->getFinestLevelNumber());
    for (int i=0; i<input_data->nauxs; i++) {
        stream << "auxs(" << i << ")"; 
        var_name = stream.str();
        stream.str("");
        var = new pdat::CellVariable<NDIM,double>( var_name, 1 );
        var_id = var_db->registerVariableAndContext(var, context_x, ghost1);
        d_aux_scalar->addComponent( var, var_id, weight_id );
        var_id = var_db->registerVariableAndContext(var, context_xt, ghost2);
        d_aux_scalar_tmp->addComponent( var, var_id, weight_id );
    }
    for (int i=0; i<input_data->nauxv; i++) {
        stream << "auxv(" << i << ")"; 
        var_name = stream.str();
        stream.str("");
        var = new pdat::CellVariable<NDIM,double>( var_name, NDIM );
        var_id = var_db->registerVariableAndContext(var, context_x, ghost1);
        d_aux_vector->addComponent( var, var_id, weight_id );
        var_id = var_db->registerVariableAndContext(var, context_xt, ghost2);
        d_aux_vector_tmp->addComponent( var, var_id, weight_id );
    }
    d_aux_scalar->allocateVectorData();
    d_aux_vector->allocateVectorData();
    d_aux_scalar->setToScalar( 0.0, false );
    d_aux_vector->setToScalar( 0.0, false );
    d_aux_scalar_tmp->allocateVectorData();
    d_aux_vector_tmp->allocateVectorData();
    d_aux_scalar_tmp->setToScalar( 0.0, false );
    d_aux_vector_tmp->setToScalar( 0.0, false );
   
    // Allocate data for f_src
    d_f_src = new pdat::CellVariable<NDIM,double>( "fsrc", input_data->nvar );
    f_src_id = var_db ->registerVariableAndContext(d_f_src, context_f, ghost0);

    for (int ln=0; ln<d_hierarchy->getNumberOfLevels(); ln++) {
        tbox::Pointer<hier::PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);       
        level->allocatePatchData(f_src_id);
    }
   
    // Create patch variables
    for (int i=0; i<input_data->nvar; i++)
        u0_id[i] = d_initial->getComponentDescriptorIndex(i);
    for (int i=0; i<input_data->nvar; i++)
        u_id[i]  = d_x->getComponentDescriptorIndex(i);
    for (int i=0; i<input_data->nvar; i++)
        u_tmp_id[i]  = d_x_tmp->getComponentDescriptorIndex(i);
    for (int i=0; i<input_data->nauxs; i++)
        auxs_id[i] = d_aux_scalar->getComponentDescriptorIndex(i);
    for (int i=0; i<input_data->nauxv; i++)
        auxv_id[i] = d_aux_vector->getComponentDescriptorIndex(i);
    for (int i=0; i<input_data->nauxs; i++)
        auxs_tmp_id[i] = d_aux_scalar_tmp->getComponentDescriptorIndex(i);
    for (int i=0; i<input_data->nauxv; i++)
        auxv_tmp_id[i] = d_aux_vector_tmp->getComponentDescriptorIndex(i);

    //tbox::Pointer< hier::pdat::pixie3dData > pixie_data = new pdat::pixie3dData<NDIM>( "fsrc", input_data->nvar );

    LevelContainer *level_container;
    hier::IntVector<NDIM> gcwc;
    int N_levels = d_hierarchy->getNumberOfLevels();
    level_container_array = new LevelContainer*[N_levels];
    for ( int ln=0; ln<N_levels; ln++ ) {
        tbox::Pointer<hier::PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        // Create the level container
        level_container_array[ln] = new LevelContainer(level->getNumberOfPatches(),d_hierarchy,
            input_data->nvar,u0_id,u_id,input_data->nauxs,auxs_id,input_data->nauxv,auxv_id);
        level_container = (LevelContainer *) level_container_array[ln];
        // Create each patch object
        for (hier::PatchLevel<NDIM>::Iterator p(level); p; p++) {
            tbox::Pointer< hier::Patch<NDIM> > patch = level->getPatch(p());
            level_container->CreatePatch(p(),patch);
        }
    }

    // Setup pixie3dRefinePatchStrategy
    d_refine_strategy = new pixie3dRefinePatchStrategy();
    (d_refine_strategy)->setHierarchy(d_hierarchy);
    (d_refine_strategy)->setGridGeometry(grid_geometry);
    (d_refine_strategy)->setPixie3dHierarchyData((void **)level_container_array);
    if ( GHOST == 1 ) {
        (d_refine_strategy)->setPixie3dDataIDs(0, input_data->nvar, input_data->nauxs, input_data->nauxv,
									  u0_id, u_id, u_tmp_id, auxs_id, auxs_tmp_id, auxv_id, auxv_tmp_id );
    } else {
        (d_refine_strategy)->setPixie3dDataIDs(1, input_data->nvar, input_data->nauxs, input_data->nauxv,
									    u0_id, u_id, u_tmp_id, auxs_id, auxs_tmp_id, auxv_id, auxv_tmp_id );
    }
   
    // Copy the data from u0 to u and set the boundary conditions (needed to set the auxillary variable names)
    d_x->copyVector(d_initial,false);
    setBoundarySchedules(true);
    synchronizeVariables();

    // Get the variable names
    tbox::Pointer<hier::PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(0);
    level_container = (LevelContainer *) level_container_array[0];
    char *tmp_depVarLabels = new char[21*input_data->nvar];
    char *tmp_auxScalarLabels = new char[21*input_data->nauxs];
    char *tmp_auxVectorLabels = new char[21*input_data->nauxv];
    for (hier::PatchLevel<NDIM>::Iterator p(level); p; p++) {
        FORTRAN_NAME(get_var_names)(level_container->getPtr(p()),tmp_depVarLabels,tmp_auxScalarLabels,tmp_auxVectorLabels);
    }
    depVarLabels = new std::string[input_data->nvar];
    auxScalarLabels = new std::string[input_data->nauxs];
    auxVectorLabels = new std::string[input_data->nauxv];
    for (int i=0; i<input_data->nvar; i++)
        depVarLabels[i] = &tmp_depVarLabels[21*i];
    for (int i=0; i<input_data->nauxs; i++)
        auxScalarLabels[i] = &tmp_auxScalarLabels[21*i];
    for (int i=0; i<input_data->nauxv; i++)
        auxVectorLabels[i] = &tmp_auxVectorLabels[21*i];
    delete [] tmp_depVarLabels;
    delete [] tmp_auxScalarLabels;
    delete [] tmp_auxVectorLabels;

    // Register the plot quantities
    for (int i=0; i<input_data->nvar; i++) {
        if ( depVarLabels[i].empty() )
            continue;
        d_VizWriter->registerPlotQuantity(depVarLabels[i], "SCALAR", u_id[i]);
    }
    for (int i=0; i<input_data->nauxs; i++) {
        if ( auxScalarLabels[i].empty() )
            continue;
        if ( auxScalarLabels[i].compare("zeros")==0 || auxScalarLabels[i].compare("ones")==0 )
            continue;
        d_VizWriter->registerPlotQuantity(auxScalarLabels[i], "SCALAR", auxs_id[i]);
    }
    for (int i=0; i<input_data->nauxv; i++) {
        if ( auxVectorLabels[i].empty() )
            continue;
        if ( auxVectorLabels[i].compare("zeros")==0 || auxVectorLabels[i].compare("ones")==0 )
            continue;
        //d_VizWriter->registerPlotQuantity(auxVectorLabels[i], "VECTOR", auxv_id);
        d_VizWriter->registerPlotQuantity(auxVectorLabels[i]+"_x", "SCALAR", auxv_id[i], 0);
        d_VizWriter->registerPlotQuantity(auxVectorLabels[i]+"_y", "SCALAR", auxv_id[i], 1);
        d_VizWriter->registerPlotQuantity(auxVectorLabels[i]+"_z", "SCALAR", auxv_id[i], 2);
    }

}


/***********************************************************************
*                                                                      *
* Evaluate initial conditions.                                         *
*                                                                      *
***********************************************************************/
void pixie3dApplication::setInitialConditions( const double initial_time )
{
    // We have already set u0 and created the patch container object when we created the level container

    // Copy the data from u0 to u 
    d_x->copyVector(d_initial,false);

    //set the boundary schedules before forming an equilibrium
    setBoundarySchedules(true);
    // Apply boundary conditions
    synchronizeVariables();
   
    // Form initial conditions
    // Loop through hierarchy
    for ( int ln=0; ln<d_hierarchy->getNumberOfLevels(); ln++ ) {
        tbox::Pointer<hier::PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        // Get the Level container
        LevelContainer *level_container = (LevelContainer *) level_container_array[ln];
        // Loop through the different patches
        for (hier::PatchLevel<NDIM>::Iterator p(level); p; p++) {
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

    // Copy the data from u to u_ic
    d_x_ic->copyVector(d_x,false);

    //tbox::pout << "max(d_x) (initial) = " << d_x->max(true) << ", " << d_x->max(false) <<std::endl; 
    //tbox::pout << "max(d_aux_scalar) (initial) = " << d_aux_scalar->max(true) << ", " << d_aux_scalar->max(false) <<std::endl; 
    //tbox::pout << "max(d_aux_vector) (initial) = " << d_aux_vector->max(true) << ", " << d_aux_vector->max(false) <<std::endl; 


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
* Empty implementation.                                                *
*                                                                      *
***********************************************************************/
void pixie3dApplication::setValuesOnNewLevel( tbox::Pointer< hier::PatchLevel<NDIM> > level )
{
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
pixie3dApplication::apply( tbox::Pointer< solv::SAMRAIVectorReal<NDIM,double> >  &f,
			   tbox::Pointer< solv::SAMRAIVectorReal<NDIM,double> >  &x,
			   tbox::Pointer< solv::SAMRAIVectorReal<NDIM,double> >  &r,
			   double a, double b)
{

    for ( int ln=0; ln<d_hierarchy->getNumberOfLevels(); ln++ ) {
        tbox::Pointer<hier::PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        for (hier::PatchLevel<NDIM>::Iterator p(level); p; p++) {
            tbox::Pointer< hier::Patch<NDIM> > patch = level->getPatch(p());
            tbox::Pointer< pdat::CellData<NDIM,double> > tmp_x = patch->getPatchData(x->getComponentDescriptorIndex(0));
            tbox::Pointer< pdat::CellData<NDIM,double> > tmp_r = patch->getPatchData(r->getComponentDescriptorIndex(0));
            double tmp2 = 0.0;
        }
    }

    // Copy x
    //d_x->setToScalar( -1.0e20, false );
    //d_aux_scalar->setToScalar( -1.0e20, false );         // ones and zeros vectors cause problems
    //d_aux_vector->setToScalar( -1.0e20, false );
    //d_aux_scalar_tmp->setToScalar( -1.0e20, false );
    //d_aux_vector_tmp->setToScalar( -1.0e20, false );
    d_x->copyVector(x);

    // Coarsen and Refine x
  
    // Apply boundary conditions
    synchronizeVariables();
  

    //tbox::pout << "max(d_x) = " << d_x->max(true) << ", " << d_x->max(false) <<std::endl; 
    //tbox::pout << "max(d_aux_scalar) = " << d_aux_scalar->max(true) << ", " << d_aux_scalar->max(false) <<std::endl; 
    //tbox::pout << "max(d_aux_vector) = " << d_aux_vector->max(true) << ", " << d_aux_vector->max(false) <<std::endl; 
    //tbox::pout << "min(d_x) = " << d_x->min(true) << ", " << d_x->min(false) <<std::endl; 
    //tbox::pout << "min(d_aux_scalar) = " << d_aux_scalar->min(true) << ", " << d_aux_scalar->min(false) <<std::endl; 
    //tbox::pout << "min(d_aux_vector) = " << d_aux_vector->min(true) << ", " << d_aux_vector->min(false) <<std::endl; 
    

    if(d_debug_print_info_level>5) {
        tbox::pout << "*****************************************" << std::endl; 
        tbox::pout << "pixie3dApplication::apply() d_x vector after synchronizeVariables(): " << std::endl; 
        d_x->print(tbox::pout, false);
        tbox::pout << "*****************************************" << std::endl; 
    }

    // Call EvaluateFunction
    for (int i=0; i<input_data->nvar; i++)
        f_id[i] = d_x_r->getComponentDescriptorIndex(i);
    tbox::Pointer<geom::CartesianGridGeometry<NDIM> > grid_geometry = d_hierarchy->getGridGeometry();
    const SAMRAI::hier::BoxArray<NDIM> &physicalDomain = grid_geometry->getPhysicalDomain();
    int nbox[NDIM];
    for (int i=0; i<NDIM; i++)
        nbox[i] = physicalDomain[0].numberCells(i);
    // Loop through hierarchy
    for ( int ln=0; ln<d_hierarchy->getNumberOfLevels(); ln++ ) {
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
        	varrayContainer *varray = new varrayContainer(patch,input_data->nvar,f_id);
            // Get the patch container
            const hier::IntVector<NDIM> ratio = level->getRatio();
            int ng[NDIM];
            for (int i=0; i<NDIM; i++)
                ng[i] = nbox[i]*ratio(i);
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
    
//    tbox::pout << "min(d_x) = " << d_x->min(false) << std::endl; 
//    tbox::pout << "min(d_aux_scalar) = " << d_aux_scalar->min(false) << std::endl; 
//    tbox::pout << "min(d_aux_vector) = " << d_aux_vector->min(false) << std::endl; 
//    tbox::pout << "min(d_x) (interior) = " << d_x->min(true) << std::endl; 

    // Copy r
    r->copyVector(d_x_r);
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
    tbox::Pointer< hier::Variable<NDIM> > var0;
    int data_id;
    tbox::Pointer< geom::CartesianGridGeometry <NDIM> > grid_geometry = d_hierarchy->getGridGeometry();
    // Loop through hierarchy
    tbox::Pointer< hier::Variable< NDIM > > var;
   
    // Copy the data from the single ghost cell width variables to the multiple ghost cell width variables
    if ( GHOST != 1 ) {
        // copy the vectors componentwise
        d_x_tmp->copyVector(d_x, false);
        d_aux_scalar_tmp->copyVector(d_aux_scalar, false);
        d_aux_vector_tmp->copyVector(d_aux_vector, false);
    }

    // Fill ghost cells
    // moving from coarser to finer levels fill boundary conditions
    for ( int ln=0; ln<d_hierarchy->getNumberOfLevels(); ln++ ) {
        tbox::Pointer<hier::PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        hier::PatchLevel<NDIM>::Iterator p(level);
        // Get the Level container
        LevelContainer *level_container = (LevelContainer *) level_container_array[ln];
        int numberOfBoundarySequenceGroups = -1;
        void *pixiePatchData = level_container->getPtr(p());    // This is an arbitrary patch to give us the number of boundary sequency groups
        assert(pixiePatchData!=NULL);
        FORTRAN_NAME(getnumberofbcgroups)(pixiePatchData, &numberOfBoundarySequenceGroups);
        assert(numberOfBoundarySequenceGroups>0);
        // process the boundary sequence groups in order
        for( int iSeq=1; iSeq<=numberOfBoundarySequenceGroups; iSeq++) {
            // initialize the aux variable on all patches on the level before interpolating
            // coarse values up and sync-ing periodic/sibling boundaries
            for (hier::PatchLevel<NDIM>::Iterator ip(level); ip; ip++) {
                tbox::Pointer< hier::Patch<NDIM> > patch_ip = level->getPatch(ip());
                pixiePatchData = level_container->getPtr(ip());
                FORTRAN_NAME(initializeauxvar)(pixiePatchData, &iSeq);
	        }
            (d_refine_strategy)->setRefineStrategyDataId(iSeq);            
#ifdef absoft
	        FORTRAN_NAME(GETBC)(pixiePatchData, &iSeq, &d_NumberOfBoundaryConditions,&d_BoundaryConditionSequence);
#else
	        FORTRAN_NAME(getbc)(pixiePatchData, &iSeq, &d_NumberOfBoundaryConditions,&d_BoundaryConditionSequence);
#endif
	        xfer::RefineAlgorithm<NDIM> refineVariableAlgorithm;
            for( int i=0; i<d_NumberOfBoundaryConditions; i++) {
                int cId = d_BoundaryConditionSequence[i];
                int varType = d_BoundaryConditionSequence[i+d_NumberOfBoundaryConditions];
                int idx = abs(cId)-1;
       	        // call applybc for interior patches
       	        // we are forced to do this on the coarsest level for individual patches because if a patch touches
                // no physical boundary conditions then the apply bc might not be called for it
                // which might result in some boundaries not being set correctly
                // on finer levels this detail is taken care of (need to confirm)
                // by the postprocessRefine call
	            if(varType==0) {
                    if ( GHOST == 1 ) {
                        var0 = (cId>0)? d_x->getComponentVariable(idx) : d_aux_scalar->getComponentVariable(idx);
                        data_id = (cId>0) ? d_x->getComponentDescriptorIndex(idx):d_aux_scalar->getComponentDescriptorIndex(idx) ;
                    } else {
                        var0 = (cId>0)? d_x_tmp->getComponentVariable(idx) : d_aux_scalar_tmp->getComponentVariable(idx);
                        data_id = (cId>0) ? d_x_tmp->getComponentDescriptorIndex(idx):d_aux_scalar_tmp->getComponentDescriptorIndex(idx) ;
                    }
                    refineVariableAlgorithm.registerRefine( data_id, data_id, data_id,
							   grid_geometry->lookupRefineOperator(var0,d_refine_op_str) );
                }
                // next do the registerRefine for a vector of components
	            if((varType==1)&&(cId>0)) {
                    // a vector will be composed of NDIM scalar components for the dependent variables
                    for(int j=0; j<NDIM; j++) {
                        if ( GHOST == 1 ) {
                            var0 = d_x->getComponentVariable(idx+j);
                            data_id = d_x->getComponentDescriptorIndex(idx+j);
                        } else {
                            var0 = d_x_tmp->getComponentVariable(idx+j);
                            data_id = d_x_tmp->getComponentDescriptorIndex(idx+j);
                        }
                        refineVariableAlgorithm.registerRefine( data_id, data_id, data_id,
						           grid_geometry->lookupRefineOperator(var0,d_refine_op_str) );
                    }
                }
                // next do the registerRefine for an auxillary vector
                if((varType==1)&&(cId<0)) {
                    if ( GHOST == 1 ) {
                        var0 = d_aux_vector->getComponentVariable(idx);
                        data_id = d_aux_vector->getComponentDescriptorIndex(idx);
                    } else {
                        var0 = d_aux_vector_tmp->getComponentVariable(idx);
                        data_id = d_aux_vector_tmp->getComponentDescriptorIndex(idx);
                    }
                    refineVariableAlgorithm.registerRefine( data_id, data_id, data_id,
						        grid_geometry->lookupRefineOperator(var0,d_refine_op_str) );
                }
            }
            tbox::Pointer< xfer::RefineSchedule<NDIM> > schedule = refineVariableAlgorithm.createSchedule(
                                level, ln-1, d_hierarchy, d_refine_strategy);
            schedule->fillData(0.0);
            schedule.setNull();

            // Fill corners and edges
            for( int i=0; i<d_NumberOfBoundaryConditions; i++) {
                int cId = d_BoundaryConditionSequence[i];
                int varType = d_BoundaryConditionSequence[i+d_NumberOfBoundaryConditions];
                int idx = abs(cId)-1;
                // Refine interior patches
	            if(varType==0) {
                    if ( GHOST == 1 ) {
                        var0 = (cId>0)? d_x->getComponentVariable(idx) : d_aux_scalar->getComponentVariable(idx);
                        data_id = (cId>0) ? d_x->getComponentDescriptorIndex(idx):d_aux_scalar->getComponentDescriptorIndex(idx) ;
                    } else {
                        var0 = (cId>0)? d_x_tmp->getComponentVariable(idx) : d_aux_scalar_tmp->getComponentVariable(idx);
                        data_id = (cId>0) ? d_x_tmp->getComponentDescriptorIndex(idx):d_aux_scalar_tmp->getComponentDescriptorIndex(idx) ;
                    }
                    xfer::SiblingGhostAlgorithm<NDIM> siblingGhostAlg;
                    siblingGhostAlg.registerSiblingGhost(data_id,data_id,data_id);
                    tbox::Pointer< xfer::SiblingGhostSchedule<NDIM> > siblingSchedule = siblingGhostAlg.createSchedule(level);
                    siblingSchedule->fillData(0.0);
                }
                // Refine for a vector of components
	            if((varType==1)&&(cId>0)) {
                    // a vector will be composed of NDIM scalar components for the dependent variables
                    for(int j=0; j<NDIM; j++) {
                        if ( GHOST == 1 ) {
                            var0 = d_x->getComponentVariable(idx+j);
                            data_id = d_x->getComponentDescriptorIndex(idx+j);
                        } else {
                            var0 = d_x_tmp->getComponentVariable(idx+j);
                            data_id = d_x_tmp->getComponentDescriptorIndex(idx+j);
                        }
                        xfer::SiblingGhostAlgorithm<NDIM> siblingGhostAlg;
                        siblingGhostAlg.registerSiblingGhost(data_id,data_id,data_id);
                        tbox::Pointer< xfer::SiblingGhostSchedule<NDIM> > siblingSchedule = siblingGhostAlg.createSchedule(level);
                        siblingSchedule->fillData(0.0);
                    }
                }
                // Refine for an auxillary vector
                if((varType==1)&&(cId<0)) {
                    if ( GHOST == 1 ) {
                        var0 = d_aux_vector->getComponentVariable(idx);
                        data_id = d_aux_vector->getComponentDescriptorIndex(idx);
                    } else {
                        var0 = d_aux_vector_tmp->getComponentVariable(idx);
                        data_id = d_aux_vector_tmp->getComponentDescriptorIndex(idx);
                    }
                    xfer::SiblingGhostAlgorithm<NDIM> siblingGhostAlg;
                    siblingGhostAlg.registerSiblingGhost(data_id,data_id,data_id);
                    tbox::Pointer< xfer::SiblingGhostSchedule<NDIM> > siblingSchedule = siblingGhostAlg.createSchedule(level);
                    siblingSchedule->fillData(0.0);
                }
            }


            /*// Fill corners and edges for all variables (only for testing purposes)
            tbox::Pointer<hier::PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
            for (int i=0; i<input_data->nvar; i++) {
                xfer::SiblingGhostAlgorithm<NDIM> siblingGhostAlg;
                siblingGhostAlg.registerSiblingGhost(u_id[i],u_id[i],u_id[i]);
                tbox::Pointer< xfer::SiblingGhostSchedule<NDIM> > siblingSchedule = siblingGhostAlg.createSchedule(level);
                siblingSchedule->fillData(0.0);
            }
            for (int i=0; i<input_data->nauxs; i++) {
                xfer::SiblingGhostAlgorithm<NDIM> siblingGhostAlg;
                siblingGhostAlg.registerSiblingGhost(auxs_id[i],auxs_id[i],auxs_id[i]);
                tbox::Pointer< xfer::SiblingGhostSchedule<NDIM> > siblingSchedule = siblingGhostAlg.createSchedule(level);
                siblingSchedule->fillData(0.0);
            }
            for (int i=0; i<input_data->nauxv; i++) {
                xfer::SiblingGhostAlgorithm<NDIM> siblingGhostAlg;
                siblingGhostAlg.registerSiblingGhost(auxv_id[i],auxv_id[i],auxv_id[i]);
                tbox::Pointer< xfer::SiblingGhostSchedule<NDIM> > siblingSchedule = siblingGhostAlg.createSchedule(level);
                siblingSchedule->fillData(0.0);
            }*/


        }
    }
   

     // Copy the data from the multiple ghost cell width variables to the single ghost cell width variables
     if ( GHOST != 1 ) {
         // copy the vectors componentwise
         d_x->copyVector(d_x_tmp, false);
         d_aux_scalar->copyVector(d_aux_scalar_tmp, false);
         d_aux_vector->copyVector(d_aux_vector_tmp, false);
     }  
}


void
pixie3dApplication::setBoundarySchedules(bool bIsInitialTime)
{
    d_bIsInitialTime = bIsInitialTime;
  
   int it=(bIsInitialTime)?0:1;
   it = 0;

    (d_refine_strategy)->setPixie3DTime(it);

    // Get d_BoundaryConditionSequence
    // we need to call this before the call to refine as the bc schedules have to be set correctly
  
    // Loop through hierarchy
    for ( int ln=0; ln<d_hierarchy->getNumberOfLevels(); ln++ ) {
        tbox::Pointer<hier::PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        // Get the Level container
        LevelContainer *level_container = (LevelContainer *) level_container_array[ln];
        // Loop through the different patches
        for (hier::PatchLevel<NDIM>::Iterator p(level); p; p++) {
            #if 0	   
                // Get d_BoundaryConditionSequence
                #ifdef absoft
                    FORTRAN_NAME(GETBC)(level_container->getPtr(p()),&d_NumberOfBoundaryConditions,&d_BoundaryConditionSequence, &it);
                #else
                    FORTRAN_NAME(getbc)(level_container->getPtr(p()),&d_NumberOfBoundaryConditions,&d_BoundaryConditionSequence, &it);
                #endif
            #else
                #ifdef absoft
                    FORTRAN_NAME()(level_container->getPtr(p()),&d_NumberOfBoundaryConditions,&d_BoundaryConditionSequence, &it);
                #else
                    FORTRAN_NAME(setupvarinitseq)(level_container->getPtr(p()),&it);
                #endif   
            #endif
        }
    }

#if 0   
   for (int i=0; i<d_NumberOfBoundaryConditions; i++)
     {
       tbox::pout << "BC(" << i+1 << ",:) = " << d_BoundaryConditionSequence[i] << ", " << d_BoundaryConditionSequence[i+d_NumberOfBoundaryConditions] << "\n";
     }
#endif
   
}

/***********************************************************************
* Apply the coarse and refine operators                                *
***********************************************************************/
void pixie3dApplication::synchronizeVariables(void)
{
  if(!d_RefineSchedulesGenerated)
    {
      generateTransferSchedules();
      d_RefineSchedulesGenerated=true;
    }
  
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

   // do registerRefine for a scalar refine algorithm
   tbox::Pointer< hier::Variable<NDIM> > var0;

   const int hierarchy_size = d_hierarchy->getNumberOfLevels();
   
   d_cell_coarsen_schedules.resizeArray(hierarchy_size-1);

   int ln=0;
   tbox::Pointer<hier::PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);

   for (ln=0; ln<hierarchy_size-1; ln++)
     {
       tbox::Pointer<hier::PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
       
       if(ln<hierarchy_size-1)
	 {
	   tbox::Pointer<hier::PatchLevel<NDIM> > flevel = d_hierarchy->getPatchLevel(ln+1);         
	   d_cell_coarsen_schedules[ln]=d_cell_coarsen_alg.createSchedule(level,flevel);
	 }
     }   
}

int pixie3dApplication::getNumberOfDependentVariables()
{
    //  assert(data!=NULL);
    return input_data->nvar;
}



// Write a single variable to a file
void pixie3dApplication::writeCellData( FILE *fp, int var_id ) {
    #if NDIM != 3
        #error not porgramed for dimensions other than 3
    #endif
    int rank = tbox::SAMRAI_MPI::getRank();
    if ( rank==0 && fp==NULL )
        TBOX_ERROR( "File pointer is NULL" );
    // Get the ghost cell width and depth of the variable
    hier::IntVector<NDIM> gcw = hier::IntVector<NDIM>(0);
    int depth = 0;
    for ( int ln=0; ln<d_hierarchy->getNumberOfLevels(); ln++ ) {
        tbox::Pointer<hier::PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        for (hier::PatchLevel<NDIM>::Iterator p(level); p; p++) {
            tbox::Pointer< hier::Patch<NDIM> > patch = level->getPatch(p());
            tbox::Pointer< pdat::CellData<NDIM,double> > pdat = patch->getPatchData(var_id);
            gcw = pdat->getGhostCellWidth();
            depth = pdat->getDepth();
            break;
        }
        if ( depth!=0 )
            break;
    }
    if ( depth==0 && rank==0 ) {
        // No local patches and rank 0
        TBOX_ERROR( "Error saving data" );
    }
    // Loop through the patches and levels, saving the data
    for ( int ln=0; ln<d_hierarchy->getNumberOfLevels(); ln++ ) {
        tbox::Pointer<hier::PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        const hier::IntVector<NDIM> ratio = level->getRatio();
        if ( rank==0 )
            fprintf(fp,"level = %i, ratio = (%i,%i,%i), n_patch = %i\n",ln,ratio(0),ratio(1),ratio(2),level->getNumberOfPatches());
        for (int p=0; p<level->getNumberOfPatches(); p++) {
            hier::Box<NDIM> box = level->getBoxForPatch(p);
            hier::Index<NDIM> ifirst = box.lower();
            hier::Index<NDIM> ilast  = box.upper();
            int N = (ilast(0)-ifirst(0)+1+2*gcw(0))*(ilast(1)-ifirst(1)+1+2*gcw(1))*(ilast(2)-ifirst(2)+1+2*gcw(2))*depth;
            int size = N*sizeof(double)/sizeof(int);    // The number of int's needed to send the data
            if ( rank==0 ) {
                // Rank 0 is writing the data
                fprintf(fp,"patch_num = %i, ifirst = (%i,%i,%i), ilast = (%i,%i,%i), gcw = (%i,%i,%i), depth = %i\n",
                    p,ifirst(0),ifirst(1),ifirst(2),ilast(0),ilast(1),ilast(2),gcw(0),gcw(1),gcw(2),depth);
                if ( level->getMappingForPatch(p)==0 ) {
                    // The patch is on processor 0, write the data directly
                	tbox::Pointer< hier::Patch<NDIM> > patch = level->getPatch(p);
                    tbox::Pointer< pdat::CellData<NDIM,double> > pdat = patch->getPatchData(var_id);
                    double *data = pdat->getPointer();
                    fwrite(data,sizeof(double),N,fp);
                } else {
                    // The patch is on another processor, we need to communicate it, then save
                    double *data = new double[N];
                    int proc = level->getMappingForPatch(p);
                    tbox::SAMRAI_MPI::recv((int*)data,size,proc,false,p);
                    fwrite(data,sizeof(double),N,fp);
                    delete [] data;
                }
                fprintf(fp,"\n");
            } else {
                // All other ranks send the data to processor 0
                if ( level->getMappingForPatch(p)==rank ) {
                    // The patch is on current processor, send the data to processor 0
                	tbox::Pointer< hier::Patch<NDIM> > patch = level->getPatch(p);
                    tbox::Pointer< pdat::CellData<NDIM,double> > pdat = patch->getPatchData(var_id);
                    double *data = pdat->getPointer();
                    tbox::SAMRAI_MPI::send((int*)data,size,0,false,p);
                }
            }
        }
        tbox::SAMRAI_MPI::barrier();
    }
}


// Write the primary and auxillary variables (including ghost cells to binary file that can be read by MATLAB)
void pixie3dApplication::writeDebugData( FILE *fp, const int it, const double time ) {
    // Print some information about the time and the domain size
    tbox::Pointer<geom::CartesianGridGeometry<NDIM> > grid_geometry = d_hierarchy->getGridGeometry();
    const double *lower = grid_geometry->getXLower();
    const double *upper = grid_geometry->getXUpper();
    const SAMRAI::hier::BoxArray<NDIM> &physicalDomain = grid_geometry->getPhysicalDomain() ;
    int nbox[NDIM];
    for (int i=0; i<NDIM; i++)
        nbox[i] = physicalDomain[0].numberCells(i);
    if ( fp != NULL ) {
        fprintf(fp,"iteration = %i\n",it);
        fprintf(fp,"time = %e\n",time);
        fprintf(fp,"lower = (%e,%e,%e)\n",lower[0],lower[1],lower[2]);
        fprintf(fp,"upper = (%e,%e,%e)\n",upper[0],upper[1],upper[2]);
        fprintf(fp,"nbox = (%i,%i,%i)\n",nbox[0],nbox[1],nbox[2]);
        fprintf(fp,"N_levels = %i\n",d_hierarchy->getNumberOfLevels());
        fprintf(fp,"N_vars = %i\n",input_data->nvar+input_data->nauxs+input_data->nauxv);
    }
    // Save the dependent variables
    for (int i=0; i<input_data->nvar; i++) {
        if ( fp != NULL )
            fprintf(fp,"var_name = %s\n",depVarLabels[i].c_str());
        writeCellData(fp,u_id[i]);
    }
    // Save the scalar auxillary variables
    for (int i=0; i<input_data->nauxs; i++) {
        if ( fp != NULL )
            fprintf(fp,"var_name = %s\n",auxScalarLabels[i].c_str());
        writeCellData(fp,auxs_id[i]);
    }
    // Save the scalar vector variables
    for (int i=0; i<input_data->nauxv; i++) {
        if ( fp != NULL )
            fprintf(fp,"var_name = %s\n",auxVectorLabels[i].c_str());
        writeCellData(fp,auxv_id[i]);
    }
    if ( fp != NULL )
        fprintf(fp,"\n");
}


}


