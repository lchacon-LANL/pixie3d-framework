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
extern void FORTRAN_NAME(setupvarinitseq)(void*, int*);
extern void FORTRAN_NAME(getnumberofbcgroups)(void*, int*);
extern void FORTRAN_NAME(getbcsequence)(void*, int*, int*, int**);
extern void FORTRAN_NAME(initializeauxvar)(void*, int*);
extern void FORTRAN_NAME(get_var_names)(void*, char*, char*, char*);
extern void FORTRAN_NAME(findexplicitdt)(void*, double*);
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
    d_hierarchy = NULL;
    for (int i=0; i<MAX_LEVELS; i++)
        level_container_array[i] = NULL;
    d_NumberOfBoundarySequenceGroups = 0;
    d_BoundarySequenceGroups = NULL;
    for (int i=0; i<MAX_LEVELS; i++) {
        refineSchedule[i] = NULL;
        siblingSchedule[i] = NULL;
    }
}


/***********************************************************************
*                                                                      *
* Construct from parameter list.  Calls initialize.                    *
*                                                                      *
***********************************************************************/
pixie3dApplication::pixie3dApplication(  pixie3dApplicationParameters* parameters ): SAMRSolvers::DiscreteOperator(parameters)
{
    d_hierarchy = NULL;
    for (int i=0; i<MAX_LEVELS; i++)
        level_container_array[i] = NULL;
    d_NumberOfBoundarySequenceGroups = 0;
    d_BoundarySequenceGroups = NULL;
    for (int i=0; i<MAX_LEVELS; i++) {
        refineSchedule[i] = NULL;
        siblingSchedule[i] = NULL;
    }
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
        if ( level_container_array[ln] != NULL ) {
            level_container = (LevelContainer *) level_container_array[ln];
            delete level_container;
        }
    }
    // Delete the refine/coarsen schedules
    for (int i=0; i<MAX_LEVELS; i++) {
        if ( refineSchedule[i] != NULL ) {
            for (int j=0; j<d_NumberOfBoundarySequenceGroups; j++)
                refineSchedule[i][j].setNull();
            delete [] refineSchedule[i];
        }
        if ( siblingSchedule[i] != NULL ) {
            for (int j=0; j<d_NumberOfBoundarySequenceGroups; j++)
                siblingSchedule[i][j].setNull();
            delete [] siblingSchedule[i];
        }
    }
    if ( d_BoundarySequenceGroups != NULL ) {
        delete [] d_BoundarySequenceGroups;
        d_BoundarySequenceGroups = NULL;
    }
    d_NumberOfBoundarySequenceGroups = 0;
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
   d_coarsen_op_str = "CONSERVATIVE_COARSEN";
   //d_coarsen_op_str = "CELL_DOUBLE_INJECTION_COARSEN";
   //d_refine_op_str  = "CONSTANT_REFINE";
   d_refine_op_str  = "LINEAR_REFINE";   
   //d_refine_op_str  = "CELL_DOUBLE_CUBIC_REFINE";
   d_RefineSchedulesGenerated=false;

    // Load basic information from the parameters
    #ifdef DEBUG_CHECK_ASSERTIONS
        assert( parameters != (pixie3dApplicationParameters*) NULL );
    #endif
    d_object_name = "pixie3d";
    d_hierarchy = parameters->d_hierarchy;
    d_VizWriter = parameters->d_VizWriter;
   
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
   
/*    if ( input_data->nxd != nbox[0] )
        TBOX_ERROR("The domain size does not match pixie3d in x-direction");
    if ( input_data->nyd != nbox[1] )
        TBOX_ERROR("The domain size does not match pixie3d in y-direction");
    if ( input_data->nzd != nbox[2] )
        TBOX_ERROR("The domain size does not match pixie3d in z-direction");
*/   

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
        u_tmp_id[i] = d_x_tmp->getComponentDescriptorIndex(i);
    for (int i=0; i<input_data->nauxs; i++)
        auxs_id[i] = d_aux_scalar->getComponentDescriptorIndex(i);
    for (int i=0; i<input_data->nauxv; i++)
        auxv_id[i] = d_aux_vector->getComponentDescriptorIndex(i);
    for (int i=0; i<input_data->nauxs; i++)
        auxs_tmp_id[i] = d_aux_scalar_tmp->getComponentDescriptorIndex(i);
    for (int i=0; i<input_data->nauxv; i++)
        auxv_tmp_id[i] = d_aux_vector_tmp->getComponentDescriptorIndex(i);

    //tbox::Pointer< hier::pdat::pixie3dData > pixie_data = new pdat::pixie3dData<NDIM>( "fsrc", input_data->nvar );

   
    /*LevelContainer *level_container;
    hier::IntVector<NDIM> gcwc;
    int N_levels = d_hierarchy->getNumberOfLevels();
    if ( N_levels > MAX_LEVELS )
        TBOX_ERROR("Maximum number of levels exceeded");
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
    }*/

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
   
    // Reset the Hierarchy Configuration (this will create the level container and initialize the communication schedules)
    for ( int ln=0; ln<d_hierarchy->getNumberOfLevels(); ln++ ) {
        const tbox::Pointer<hier::BasePatchHierarchy<NDIM> > hierarchy = d_hierarchy;
        resetHierarchyConfiguration(hierarchy,0,d_hierarchy->getFinestLevelNumber());
    }

    // Copy the data from u0 to u and set the boundary conditions (needed to set the auxillary variable names)
    d_x->copyVector(d_initial,false);
    synchronizeVariables();

    // Get the variable names
    tbox::Pointer<hier::PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(0);
    LevelContainer *level_container = (LevelContainer *) level_container_array[0];
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
        if ( auxScalarLabels[i].compare("NULL")==0 )
            continue;
        d_VizWriter->registerPlotQuantity(auxScalarLabels[i], "SCALAR", auxs_id[i]);
    }
    for (int i=0; i<input_data->nauxv; i++) {
        if ( auxVectorLabels[i].empty() )
            continue;
        if ( auxVectorLabels[i].compare("zeros")==0 || auxVectorLabels[i].compare("ones")==0 )
            continue;
        if ( auxVectorLabels[i].compare("NULL")==0 )
            continue;
        //d_VizWriter->registerPlotQuantity(auxVectorLabels[i], "VECTOR", auxv_id);
        d_VizWriter->registerPlotQuantity(auxVectorLabels[i]+"_x", "SCALAR", auxv_id[i], 0);
        d_VizWriter->registerPlotQuantity(auxVectorLabels[i]+"_y", "SCALAR", auxv_id[i], 1);
        d_VizWriter->registerPlotQuantity(auxVectorLabels[i]+"_z", "SCALAR", auxv_id[i], 2);
    }
    for (int i=0; i<input_data->nvar; i++) {
        if ( depVarLabels[i].empty() )
            continue;
        d_VizWriter->registerPlotQuantity(depVarLabels[i]+"_r", "SCALAR", d_x_r->getComponentDescriptorIndex(i));
    }

}


/***********************************************************************
*                                                                      *
* Evaluate initial conditions.                                         *
*                                                                      *
***********************************************************************/
void pixie3dApplication::setInitialConditions( const double initial_time )
{
    // We have already set u0 we created the level container on the resetHierarchyConfiguration

    // Copy the data from u0 to u 
    d_x->copyVector(d_initial,false);

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

    // Apply boundary conditions
    synchronizeVariables();   

    // Copy the data from u to u_ic
    d_x_ic->copyVector(d_x,false);

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

    // Copy x
    d_x->copyVector(x);

    // Coarsen and Refine x
  
    // Apply boundary conditions
    synchronizeVariables();
  
    if(d_debug_print_info_level>5) {
        tbox::pout << "*****************************************" << std::endl; 
        tbox::pout << "pixie3dApplication::apply() d_x vector after synchronizeVariables(): " << std::endl; 
        d_x->print(tbox::pout, false);
        tbox::pout << "*****************************************" << std::endl; 
    }

    // Call EvaluateFunction
    double dt_exp = 1e10;
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
        	// Create f
            #ifdef absoft
            	FORTRAN_NAME(EVALUATENONLINEARRESIDUAL)(level_container->getPtr(p()),&n_elem,fsrc,varray->getPtr());
            #else
            	FORTRAN_NAME(evaluatenonlinearresidual)(level_container->getPtr(p()),&n_elem,fsrc,varray->getPtr());
            #endif
            // Comupute the timestep required for an explicit method, computed by pixie3d
            double dt_patch;
            #ifdef absoft
            	FORTRAN_NAME(FINDEXPLICITDT)(level_container->getPtr(p()),&n_elem,fsrc,varray->getPtr());
            #else
              	FORTRAN_NAME(findexplicitdt)(level_container->getPtr(p()),&dt_patch);
            #endif
            if ( dt_patch < dt_exp )
                dt_exp = dt_patch;
        	// Delete varray
        	delete varray;
        }
    }
    // Get the global minimum timestep
    dt_exp = tbox::SAMRAI_MPI::minReduction(dt_exp);

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
#if NDIM == 1
         // 1d version
         const int depth = tmp->getDepth();
         const hier::IntVector<NDIM> ghost = tmp->getGhostCellWidth();
         double *data = tmp->getPointer();
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
         const int depth = tmp->getDepth();
         const hier::IntVector<NDIM> ghost = tmp->getGhostCellWidth();
         double *data = tmp->getPointer();
         tbox::pout << "  Patch # = " << p << ", ifirst = " << ifirst << ", ilast = " << ilast 
            << ", ghost = " << ghost << ", xlo = (" << xlo[0] << "," << xlo[1] << ")" 
            << ", dx = (" << dx[0] << "," << dx[1] << ")\n";
         int num1 = ilast(0)-ifirst(0)+1+2*ghost(0);
         int num2 = ilast(1)-ifirst(1)+1+2*ghost(1);
         for (int k=0; k<depth; k++) {
            for (int j=0; j<num2; j++) {
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
    for ( int ln=d_hierarchy->getFinestLevelNumber()-1; ln>=0; ln-- ) {
        coarsenSchedule[ln]->coarsenData();
    }
}


/***********************************************************************
* Refine the data                                                      *
***********************************************************************/
void  pixie3dApplication::refineVariables(void)
{
    tbox::Pointer< hier::Variable<NDIM> > var0;
    tbox::Pointer< geom::CartesianGridGeometry <NDIM> > grid_geometry = d_hierarchy->getGridGeometry();
    tbox::Pointer< hier::Variable< NDIM > > var;
    tbox::Pointer<hier::PatchLevel<NDIM> > level;
    LevelContainer *level_container;
    void *pixiePatchData;

    // Copy the data from the single ghost cell width variables to the multiple ghost cell width variables
    if ( GHOST != 1 ) {
        // copy the vectors componentwise
        d_x_tmp->copyVector(d_x, false);
        d_aux_scalar_tmp->copyVector(d_aux_scalar, false);
        d_aux_vector_tmp->copyVector(d_aux_vector, false);
    }

  /*  // Check the number of boundary condition groups and the sequence for each group
    level = d_hierarchy->getPatchLevel(0);
    hier::PatchLevel<NDIM>::Iterator p(level);
    level_container = (LevelContainer *) level_container_array[0];
    pixiePatchData = level_container->getPtr(p());    // This is an arbitrary patch to give us the number of boundary sequency groups
    assert(pixiePatchData!=NULL);
    int tmp_NumberOfBoundarySequenceGroups;
    FORTRAN_NAME(getnumberofbcgroups)(pixiePatchData,&tmp_NumberOfBoundarySequenceGroups);
    if ( tmp_NumberOfBoundarySequenceGroups != d_NumberOfBoundarySequenceGroups )
        TBOX_ERROR( "The number of boundary sequency groups changed" );
    for( int iSeq=0; iSeq<d_NumberOfBoundarySequenceGroups; iSeq++) {
        int tmp_NumberOfBoundaryConditions;
        int *tmp_BoundaryConditionSequence;
        int iSeq2 = iSeq+1;     // The Fortran code starts indexing at 1
        FORTRAN_NAME(getbc)(pixiePatchData, &iSeq2, &tmp_NumberOfBoundaryConditions,&tmp_BoundaryConditionSequence);
        if ( tmp_NumberOfBoundaryConditions != d_NumberOfBoundaryConditions[iSeq] )
            TBOX_ERROR( "The number of boundary conditions in a boundary sequency group changed" );
        for (int i=0; i<2*d_NumberOfBoundaryConditions[iSeq]; i++) {
            if ( d_BoundaryConditionSequence[iSeq][i] != tmp_BoundaryConditionSequence[i] )
                TBOX_ERROR( "The boundary conditions in a boundary sequency group changed" );
        }
    }*/

    // Fill ghost cells
    // moving from coarser to finer levels fill boundary conditions
    for ( int ln=0; ln<d_hierarchy->getNumberOfLevels(); ln++ ) {
        level = d_hierarchy->getPatchLevel(ln);
        level_container = (LevelContainer *) level_container_array[ln];
        // process the boundary sequence groups in order
        for( int iSeq=0; iSeq<d_NumberOfBoundarySequenceGroups; iSeq++) {
            // initialize the aux variable on all patches on the level before interpolating
            // coarse values up and sync-ing periodic/sibling boundaries
            for (hier::PatchLevel<NDIM>::Iterator ip(level); ip; ip++) {
                tbox::Pointer< hier::Patch<NDIM> > patch_ip = level->getPatch(ip());
                pixiePatchData = level_container->getPtr(ip());
                assert(pixiePatchData!=NULL);
                int iSeq2 = iSeq+1;     // The Fortran code starts indexing at 1
                FORTRAN_NAME(initializeauxvar)(pixiePatchData, &iSeq2);
	        }
            // Fill the ghost cells and apply the boundary conditions
            (d_refine_strategy)->setRefineStrategySequence(d_BoundarySequenceGroups[iSeq]);
            refineSchedule[ln][iSeq]->fillData(0.0);
            // Fill corners and edges
            siblingSchedule[ln][iSeq]->fillData(0.0);
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
 
}

int pixie3dApplication::getNumberOfDependentVariables()
{
    //  assert(data!=NULL);
    return input_data->nvar;
}



// Write a single variable to a file for all patches, all levels, including ghostcells
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


// Write a single variable to a file using only the coarse level as a single patch
void pixie3dApplication::writeGlobalCellData( FILE *fp, int var_id ) {
    #if NDIM != 3
        #error not porgramed for dimensions other than 3
    #endif
    int rank = tbox::SAMRAI_MPI::getRank();
    if ( rank==0 && fp==NULL )
        TBOX_ERROR( "File pointer is NULL" );
    // Get the ghost cell width and depth of the variable
    hier::IntVector<NDIM> gcw = hier::IntVector<NDIM>(0);
    int depth = 0;
    tbox::Pointer<hier::PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(0);
    for (hier::PatchLevel<NDIM>::Iterator p(level); p; p++) {
        tbox::Pointer< hier::Patch<NDIM> > patch = level->getPatch(p());
        tbox::Pointer< pdat::CellData<NDIM,double> > pdat = patch->getPatchData(var_id);
        gcw = pdat->getGhostCellWidth();
        depth = pdat->getDepth();
        break;
    }
    if ( depth==0 && rank==0 ) {
        // No local patches and rank 0
        TBOX_ERROR( "Error saving data" );
    }
    // Get the physical domain
    tbox::Pointer< geom::CartesianGridGeometry <NDIM> > grid_geometry = d_hierarchy->getGridGeometry();
    const hier::BoxArray<NDIM> domain = grid_geometry->getPhysicalDomain();
    if ( domain.size() != 1 )
        TBOX_ERROR( "Domain must be a single box" );
    const hier::Box<NDIM> domain_box = domain.getBox(0);
    const hier::Index<NDIM> lower = domain_box.lower();
    const hier::Index<NDIM> upper = domain_box.upper();
    // Loop through the patches, saving all data on processor 0
    int Ngx = upper(0)-lower(0)+1;
    int Ngy = upper(1)-lower(1)+1;
    int Ngz = upper(2)-lower(2)+1;
    double *data = NULL;
    if ( rank == 0 )
        data = new double[Ngx*Ngy*Ngz*depth];
    for (int p=0; p<level->getNumberOfPatches(); p++) {
        hier::Box<NDIM> box = level->getBoxForPatch(p);
        hier::Index<NDIM> ifirst = box.lower();
        hier::Index<NDIM> ilast  = box.upper();
        int Nx = ilast(0)-ifirst(0)+1+2*gcw(0);
        int Ny = ilast(1)-ifirst(1)+1+2*gcw(1);
        int Nz = ilast(2)-ifirst(2)+1+2*gcw(2);
        int size = Nx*Ny*Nz*depth*sizeof(double)/sizeof(int);
        int proc = level->getMappingForPatch(p);
        if ( rank==0 ) {
            // Processor 0 will either contain the data, or recieve it from another processor and copy it to data
            double *tmp_data;
            if ( level->getMappingForPatch(p)==0 ) {
                tbox::Pointer< hier::Patch<NDIM> > patch = level->getPatch(p);
                tbox::Pointer< pdat::CellData<NDIM,double> > pdat = patch->getPatchData(var_id);
                tmp_data = pdat->getPointer();
            } else {
                tmp_data = new double[Nx*Ny*Nz*depth];
                tbox::SAMRAI_MPI::recv((int*)tmp_data,size,proc,false,p);
            }
            for (int i=0; i<Nx-2*gcw(0); i++) {
                for (int j=0; j<Ny-2*gcw(1); j++) {
                    for (int k=0; k<Nz-2*gcw(2); k++) {
                        for (int d=0; d<depth; d++) {
                            data[(i+ifirst(0))+(j+ifirst(1))*Ngx+(k+ifirst(2))*Ngx*Ngy+d*Ngx*Ngy*Ngz] = 
                                tmp_data[(i+gcw(0))+(j+gcw(1))*Nx+(k+gcw(2))*Nx*Ny+d*Nx*Ny*Nz];
                        }
                    }
                }
            }
            if ( level->getMappingForPatch(p)!=0 )
                delete [] tmp_data;
        } else if ( level->getMappingForPatch(p)==rank ) {
            // The current processor has the data, send it to processor 0
            tbox::Pointer< hier::Patch<NDIM> > patch = level->getPatch(p);
            tbox::Pointer< pdat::CellData<NDIM,double> > pdat = patch->getPatchData(var_id);
            double *tmp_data = pdat->getPointer();
            tbox::SAMRAI_MPI::send((int*)tmp_data,size,0,false,p);
        } else { 
            // The current processor does not contain the data and is not processor 0, do nothing
        }
    }
    // Save the results
    if ( rank == 0 ) {
        fprintf(fp," depth = %i\n",depth);
        fwrite(data,sizeof(double),Ngx*Ngy*Ngz*depth,fp);
        fprintf(fp,"\n");
        delete [] data;
    }
    tbox::SAMRAI_MPI::barrier();
}


/**************************************************************************
* Write the primary and auxillary variables.                              *
* type = 1:  Write each patch with ghost cells for all levels             *
* type = 2:  Write the coarse level without ghost cells as a single patch *
**************************************************************************/
void pixie3dApplication::writeDebugData( FILE *fp, const int it, const double time, int type ) {
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
        if ( type == 1 )
            fprintf(fp,"N_levels = %i\n",d_hierarchy->getNumberOfLevels());
        else if ( type == 2 )
            fprintf(fp,"N_levels = %i\n",-1);
        fprintf(fp,"N_vars = %i\n",input_data->nvar+input_data->nauxs+input_data->nauxv);
    }
    // Save the dependent variables
    for (int i=0; i<input_data->nvar; i++) {
        if ( fp != NULL )
            fprintf(fp,"var_name = %s\n",depVarLabels[i].c_str());
        if ( type == 1 )
            writeCellData(fp,u_id[i]);
        else if ( type == 2 ) 
            writeGlobalCellData(fp,u_id[i]);
    }
    // Save the scalar auxillary variables
    for (int i=0; i<input_data->nauxs; i++) {
        if ( fp != NULL )
            fprintf(fp,"var_name = %s\n",auxScalarLabels[i].c_str());
        if ( type == 1 )
            writeCellData(fp,auxs_id[i]);
        else if ( type == 2 ) 
            writeGlobalCellData(fp,auxs_id[i]);
    }
    // Save the scalar vector variables
    for (int i=0; i<input_data->nauxv; i++) {
        if ( fp != NULL )
            fprintf(fp,"var_name = %s\n",auxVectorLabels[i].c_str());
        if ( type == 1 )
            writeCellData(fp,auxv_id[i]);
        else if ( type == 2 ) 
            writeGlobalCellData(fp,auxv_id[i]);
    }
    if ( fp != NULL )
        fprintf(fp,"\n");
}



/***********************************************************************
* Initialize level data.                                               *
* Function overloaded from mesh::StandardTagAndInitStrategy<NDIM>.     *
***********************************************************************/
void pixie3dApplication::initializeLevelData( const tbox::Pointer<hier::BasePatchHierarchy<NDIM> > hierarchy,
    const int level_number, const double time, const bool can_be_refined, const bool initial_time,
    const tbox::Pointer<hier::BasePatchLevel<NDIM> > old_level, const bool allocate_data )
{
    // Check if the application has been initialized
    if ( d_hierarchy.isNull() )
        return;

    // Get the new level
    #ifdef DEBUG_CHECK_ASSERTIONS
        assert(!hierarchy.isNull());
        assert( (level_number >= 0) && (level_number <= hierarchy->getFinestLevelNumber()) );
        if ( !(old_level.isNull()) )
            assert( level_number == old_level->getLevelNumber() );
        assert(!(hierarchy->getPatchLevel(level_number)).isNull());
    #endif
    tbox::Pointer<hier::PatchLevel<NDIM> > level = hierarchy->getPatchLevel(level_number);

    // Check the new level for overlapping boxes
    const hier::BoxArray<NDIM> boxArray = level->getBoxes();
    tbox::Pointer<hier::BoxTree<NDIM> > boxTree = level->getBoxTree();
    for (int i=0; i<boxArray.size(); i++) {
        hier::Box<NDIM> box = boxArray.getBox(i);
        hier::BoxList<NDIM> overlap;
        boxTree->findOverlapBoxes(overlap,box);
        if ( overlap.getNumberOfBoxes() > 1 ) {
            TBOX_ERROR("Overlapping boxes were detected on new level, and are not supported");
        }
    }

    // Allocate data when called for, otherwise set timestamp on allocated data.
    hier::ComponentSelector d_problem_data = false;
    if ( !old_level.isNull() ) {
        // Use the old level to determine which components need to be allocated
        // Assume that the old level is a level on a rectangular domain (not a multiblock domain)
        const tbox::Pointer<hier::PatchLevel<NDIM> > old_level2 = old_level;
        for (int i=0; i<d_problem_data.getSize(); i++) {
            if ( old_level2->checkAllocated(i) )
                d_problem_data.setFlag(i);
        }
    } else { 
        // Use the coarsest level to determine which components need to be allocated
        tbox::Pointer<hier::PatchLevel<NDIM> > coarse_level = d_hierarchy->getPatchLevel(0);
        #ifdef DEBUG_CHECK_ASSERTIONS
            assert(!coarse_level.isNull());
        #endif
        for (int i=0; i<d_problem_data.getSize(); i++) {
            if ( coarse_level->checkAllocated(i) )
                d_problem_data.setFlag(i);
        }
    }


    if ( allocate_data )  {
        level->allocatePatchData(d_problem_data, time);
    } else  {
        level->setTime(time, d_problem_data);
    }

}



/***********************************************************************
* Reset cached information that depends on the hierarchy configuration.*
* Function overloaded from mesh::StandardTagAndInitStrategy<NDIM>.     *
***********************************************************************/
void pixie3dApplication::resetHierarchyConfiguration(
    const tbox::Pointer<hier::BasePatchHierarchy<NDIM> > hierarchy,
    const int coarsest_level, const int finest_level )
{
    // Check if the application has been initialized
    if ( d_hierarchy.isNull() )
        return;
    int N_levels = d_hierarchy->getNumberOfLevels();
    if ( N_levels > MAX_LEVELS ) 
        TBOX_ERROR("Maximum number of levels exceeded");

    // Reset the vectors
    d_initial->resetLevels(0,N_levels-1);
    d_x_tmp->resetLevels(0,N_levels-1);
    d_x->resetLevels(0,N_levels-1);
    d_x_r->resetLevels(0,N_levels-1);
    d_x_ic->resetLevels(0,N_levels-1);
    d_aux_scalar->resetLevels(0,N_levels-1);
    d_aux_vector->resetLevels(0,N_levels-1);
    d_aux_scalar_tmp->resetLevels(0,N_levels-1);
    d_aux_vector_tmp->resetLevels(0,N_levels-1);

    // Reset the level container
    LevelContainer *level_container;
    hier::IntVector<NDIM> gcwc;
    for ( int ln=0; ln<d_hierarchy->getNumberOfLevels(); ln++ ) {
        if ( level_container_array[ln] != NULL ) {
            level_container = (LevelContainer *) level_container_array[ln];
            delete level_container;
            level_container_array[ln] = NULL;
        }
    }
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

    // Reset the communication schedules
    for (int i=0; i<MAX_LEVELS; i++) {
        if ( refineSchedule[i] != NULL ) {
            for (int j=0; j<d_NumberOfBoundarySequenceGroups; j++)
                refineSchedule[i][j].setNull();
            delete [] refineSchedule[i];
        }
        if ( siblingSchedule[i] != NULL ) {
            for (int j=0; j<d_NumberOfBoundarySequenceGroups; j++)
                siblingSchedule[i][j].setNull();
            delete [] siblingSchedule[i];
        }
    }
    if ( d_BoundarySequenceGroups != NULL ) {
        delete [] d_BoundarySequenceGroups;
        d_BoundarySequenceGroups = NULL;
    }
    d_NumberOfBoundarySequenceGroups = 0;
    // Get the number of boundary condition groups and the sequence for each group
    tbox::Pointer<hier::PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(0);
    hier::PatchLevel<NDIM>::Iterator p(level);
    level_container = (LevelContainer *) level_container_array[0];
    void *pixiePatchData = level_container->getPtr(p());    // This is an arbitrary patch to give us the number of boundary sequency groups
    assert(pixiePatchData!=NULL);
    FORTRAN_NAME(getnumberofbcgroups)(pixiePatchData,&d_NumberOfBoundarySequenceGroups);
    d_BoundarySequenceGroups = new pixie3dRefinePatchStrategy::bcgrp_struct[d_NumberOfBoundarySequenceGroups];
    for( int iSeq=0; iSeq<d_NumberOfBoundarySequenceGroups; iSeq++) {
        int iSeq2 = iSeq+1;     // The Fortran code starts indexing at 1
        int N_sequence;
        int *tmp = NULL;
        FORTRAN_NAME(getbcsequence)(pixiePatchData,&iSeq2,&N_sequence,&tmp);
        pixie3dRefinePatchStrategy::bcgrp_struct tmp2(N_sequence);
        d_BoundarySequenceGroups[iSeq] = tmp2;
        for (int i=0; i<N_sequence; i++) {
            d_BoundarySequenceGroups[iSeq].bc_seq[i] = tmp[i];
            d_BoundarySequenceGroups[iSeq].vector[i] = tmp[i+N_sequence];
            d_BoundarySequenceGroups[iSeq].fillBC[i] = tmp[i+2*N_sequence];
        }
    }

    // Create the refineSchedule and siblingSchedule
    int data_id;
    tbox::Pointer< hier::Variable<NDIM> > var0;
    tbox::Pointer< geom::CartesianGridGeometry <NDIM> > grid_geometry = d_hierarchy->getGridGeometry();
    for ( int ln=0; ln<d_hierarchy->getNumberOfLevels(); ln++ ) {
        level = d_hierarchy->getPatchLevel(ln);
        refineSchedule[ln] = new tbox::Pointer< xfer::RefineSchedule<NDIM> >[d_NumberOfBoundarySequenceGroups];
        siblingSchedule[ln] = new tbox::Pointer< xfer::SiblingGhostSchedule<NDIM> >[d_NumberOfBoundarySequenceGroups];
        for( int iSeq=0; iSeq<d_NumberOfBoundarySequenceGroups; iSeq++) {
	        xfer::RefineAlgorithm<NDIM> refineVariableAlgorithm;
            xfer::SiblingGhostAlgorithm<NDIM> siblingGhostAlgorithm;
            // Register the variables in the current squence
            for (int i=0; i<d_BoundarySequenceGroups[iSeq].nbc_seq; i++) {
                int idx = abs(d_BoundarySequenceGroups[iSeq].bc_seq[i])-1;
                if ( d_BoundarySequenceGroups[iSeq].bc_seq[i]>0 && d_BoundarySequenceGroups[iSeq].vector[i]==0 ) {
                    // Dependent scalar variable
                    if ( GHOST == 1 ) {
                        var0 = d_x->getComponentVariable(idx);
                        data_id = d_x->getComponentDescriptorIndex(idx);
                    } else {
                        var0 = d_x_tmp->getComponentVariable(idx);
                        data_id = d_x_tmp->getComponentDescriptorIndex(idx);
                    }
                } else if ( d_BoundarySequenceGroups[iSeq].bc_seq[i]<0 && d_BoundarySequenceGroups[iSeq].vector[i]==0 ) {
                    // Auxillary scalar variable
                    for(int j=0; j<NDIM; j++) {
                        if ( GHOST == 1 ) {
                            var0 = d_aux_scalar->getComponentVariable(idx);
                            data_id = d_aux_scalar->getComponentDescriptorIndex(idx);
                        } else {
                            var0 = d_aux_scalar_tmp->getComponentVariable(idx);
                            data_id = d_aux_scalar_tmp->getComponentDescriptorIndex(idx);
                        }
                    }
                } else if ( d_BoundarySequenceGroups[iSeq].bc_seq[i]<0 && d_BoundarySequenceGroups[iSeq].vector[i]==1 ) {
                    // Auxillary vector variable
                    for(int j=0; j<NDIM; j++) {
                        if ( GHOST == 1 ) {
                            var0 = d_aux_vector->getComponentVariable(idx);
                            data_id = d_aux_vector->getComponentDescriptorIndex(idx);
                        } else {
                            var0 = d_aux_vector_tmp->getComponentVariable(idx);
                            data_id = d_aux_vector_tmp->getComponentDescriptorIndex(idx);
                        }
                    }
                } else {
                    // Unknown variable type
                    TBOX_ERROR("Bad boundary group member");
                }
                refineVariableAlgorithm.registerRefine( data_id, data_id, data_id,
                        grid_geometry->lookupRefineOperator(var0,d_refine_op_str) );
                siblingGhostAlgorithm.registerSiblingGhost(data_id,data_id,data_id);
            }
            // Create the schedules
            refineSchedule[ln][iSeq] = refineVariableAlgorithm.createSchedule(level, ln-1, d_hierarchy, d_refine_strategy);
            siblingSchedule[ln][iSeq] = siblingGhostAlgorithm.createSchedule(level);
        }
    }
    // Create the coarsenSchedule
    for ( int ln=0; ln<d_hierarchy->getFinestLevelNumber(); ln++ ) {
        xfer::CoarsenAlgorithm<NDIM> coarsenAlgorithm;
        for (int i=0; i<input_data->nvar; i++) {
            var0 = d_x->getComponentVariable(i);
            coarsenAlgorithm.registerCoarsen( u_id[i], u_id[i],
                grid_geometry->lookupCoarsenOperator(var0,d_coarsen_op_str) );
        }
        for (int i=0; i<input_data->nauxs; i++) {
            var0 = d_aux_scalar->getComponentVariable(i);
            coarsenAlgorithm.registerCoarsen( auxs_id[i], auxs_id[i],
                grid_geometry->lookupCoarsenOperator(var0,d_coarsen_op_str) );
        }
        for (int i=0; i<input_data->nauxv; i++) {
            var0 = d_aux_vector->getComponentVariable(i);
            coarsenAlgorithm.registerCoarsen( auxv_id[i], auxv_id[i],
                grid_geometry->lookupCoarsenOperator(var0,d_coarsen_op_str) );
        }
  	    tbox::Pointer<hier::PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);    
  	    tbox::Pointer<hier::PatchLevel<NDIM> > flevel = d_hierarchy->getPatchLevel(ln+1);    
        coarsenSchedule[ln] = coarsenAlgorithm.createSchedule(level,flevel);
    }

}



}


