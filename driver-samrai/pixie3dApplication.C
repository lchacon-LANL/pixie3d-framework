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


#include <iostream>
#include <sstream>
#include <string>
#include <string>
#include "math.h"

#include "pixie3dApplication.h"
#include "pixie3dRefinePatchStrategy.h"
#include "varrayContainer.h"
#include "fortran.h"

// SAMRAI headers
#include "SAMRAI/geom/CartesianGridGeometry.h"
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/geom/CartesianCellDoubleWeightedAverage.h"
#include "SAMRAI/hier/CoarsenOperator.h"
#include "SAMRAI/hier/NeighborhoodSet.h"
#include "SAMRAI/hier/MappedBox.h"
#include "SAMRAI/hier/MappedBoxSet.h"
#include "SAMRAI/xfer/PatchLevelFillPattern.h"
#include "SAMRAI/xfer/PatchLevelBorderFillPattern.h"
#include "SAMRAI/xfer/PatchLevelFullFillPattern.h"

//#include "HierarchyCellDataOpsReal.h"
//#include "BoundaryConditionStrategy.h"
//#include "RefineOperator.h"

#include "source/AMRUtilities.h"
//#include "test_utilities.h"

//#include "RefinementBoundaryInterpolation.h"
//#include "CoarsenAlgorithm.h"
//#include "CCellVariable.h"
//#include "CartesianCellDoubleCubicRefine.h"
//#include "CartesianCellDoubleLinearRefine.h"




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
pixie3dApplication::pixie3dApplication():
    dim(0)
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
    dt_exp = 1.0;
}


/***********************************************************************
*                                                                      *
* Construct from parameter list.  Calls initialize.                    *
*                                                                      *
***********************************************************************/
pixie3dApplication::pixie3dApplication(  pixie3dApplicationParameters* parameters ): 
    SAMRSolvers::DiscreteOperator(parameters),
    dim(parameters-> d_hierarchy->getDim())
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
    dt_exp = 1.0;
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
    d_hierarchy = parameters->d_hierarchy;
    if ( dim.getValue()==0 )
        dim.setValue(d_hierarchy->getDim().getValue());
    if ( dim!=d_hierarchy->getDim() )
        TBOX_ERROR("Error in dimension");
    if ( dim.getValue() != 3 )
        TBOX_ERROR("Only programmed for dimension == 3");
    d_object_name = "pixie3d";
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
    tbox::Pointer<hier::GridGeometry> grid_geometry = d_hierarchy->getGridGeometry();
    if ( grid_geometry->getNumberBlocks() != 1 )
        TBOX_ERROR("Multiblock domains are not supported");
    const SAMRAI::hier::BoxList &physicalDomainList = grid_geometry->getPhysicalDomain(0);
    if ( physicalDomainList.size() != 1 )
        TBOX_ERROR("Multiple box domains are not supported");
    const SAMRAI::hier::Box physicalDomain = physicalDomainList.getBoundingBox();
    int nbox[3];
    for (int i=0; i<dim.getValue(); i++)
        nbox[i] = physicalDomain.numberCells(i);
   
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
     * hier::Variable to weight solution vector entries on a composite grid.
     */
    
    hier::VariableDatabase* var_db = hier::VariableDatabase::getDatabase();
    tbox::Pointer<hier::VariableContext> context_x = var_db->getContext("pixie3d-x");
    tbox::Pointer<hier::VariableContext> context_xr = var_db->getContext("pixie3d-x_r");
    tbox::Pointer<hier::VariableContext> context_ic = var_db->getContext("pixie3d-x_ic");
    tbox::Pointer<hier::VariableContext> context_xt = var_db->getContext("pixie3d-x_tmp");
    tbox::Pointer<hier::VariableContext> context_in = var_db->getContext("pixie3d-initial");
    tbox::Pointer<hier::VariableContext> context_f = var_db->getContext("pixie3d-source");
    tbox::Pointer< pdat::CellVariable<double> > var;
    hier::IntVector ghost0 = hier::IntVector(dim,0);
    hier::IntVector ghost1 = hier::IntVector(dim,1);
    hier::IntVector ghost2 = hier::IntVector(dim,GHOST);
   
    var = new pdat::CellVariable<double>(dim, "weight", 1);
    const int weight_id = var_db->registerVariableAndContext(var, context_xt, ghost0);

    for (int ln=0; ln<d_hierarchy->getNumberOfLevels(); ln++) {
        tbox::Pointer<hier::PatchLevel> level = d_hierarchy->getPatchLevel(ln);       
        level->allocatePatchData(weight_id);
    }

    AMRUtilities::setVectorWeights(d_hierarchy, weight_id);

    // Allocate data for u, u_0, u_ic
    if ( IS2D == 1 )
        ghost2(2) = 1;
    int var_id;
    std::string var_name;
    d_x = new solv::SAMRAIVectorReal<double>("xVec",d_hierarchy,0,d_hierarchy->getFinestLevelNumber());
    d_x_r = new solv::SAMRAIVectorReal<double>("xVec_r",d_hierarchy,0,d_hierarchy->getFinestLevelNumber());
    d_x_ic = new solv::SAMRAIVectorReal<double>("xVecIC",d_hierarchy,0,d_hierarchy->getFinestLevelNumber());
    d_x_tmp = new solv::SAMRAIVectorReal<double>("xTmpVec",d_hierarchy,0,d_hierarchy->getFinestLevelNumber());
    d_initial = new solv::SAMRAIVectorReal<double>("xInitialVec",d_hierarchy,0,d_hierarchy->getFinestLevelNumber());
    std::stringstream stream;
    for (int i=0; i<input_data->nvar; i++) {
       stream << "x(" << i << ")"; 
       var_name = stream.str();
       stream.str("");
       var = new pdat::CellVariable<double>( dim, var_name, 1 );
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
    d_aux_scalar = new solv::SAMRAIVectorReal<double>("auxs",d_hierarchy,0,d_hierarchy->getFinestLevelNumber());
    d_aux_vector = new solv::SAMRAIVectorReal<double>("auxv",d_hierarchy,0,d_hierarchy->getFinestLevelNumber());
    d_aux_scalar_tmp = new solv::SAMRAIVectorReal<double>("auxs",d_hierarchy,0,d_hierarchy->getFinestLevelNumber());
    d_aux_vector_tmp = new solv::SAMRAIVectorReal<double>("auxv",d_hierarchy,0,d_hierarchy->getFinestLevelNumber());
    for (int i=0; i<input_data->nauxs; i++) {
        stream << "auxs(" << i << ")"; 
        var_name = stream.str();
        stream.str("");
        var = new pdat::CellVariable<double>( dim, var_name, 1 );
        var_id = var_db->registerVariableAndContext(var, context_x, ghost1);
        d_aux_scalar->addComponent( var, var_id, weight_id );
        var_id = var_db->registerVariableAndContext(var, context_xt, ghost2);
        d_aux_scalar_tmp->addComponent( var, var_id, weight_id );
    }
    for (int i=0; i<input_data->nauxv; i++) {
        stream << "auxv(" << i << ")"; 
        var_name = stream.str();
        stream.str("");
        var = new pdat::CellVariable<double>( dim, var_name, dim.getValue() );
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
    d_f_src = new pdat::CellVariable<double>( dim, "fsrc", input_data->nvar );
    f_src_id = var_db ->registerVariableAndContext(d_f_src, context_f, ghost0);

    for (int ln=0; ln<d_hierarchy->getNumberOfLevels(); ln++) {
        tbox::Pointer<hier::PatchLevel> level = d_hierarchy->getPatchLevel(ln);       
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

    //tbox::Pointer< hier::pdat::pixie3dData > pixie_data = new pdat::pixie3dData( "fsrc", input_data->nvar );

   
    /*LevelContainer *level_container;
    hier::IntVector gcwc;
    int N_levels = d_hierarchy->getNumberOfLevels();
    if ( N_levels > MAX_LEVELS )
        TBOX_ERROR("Maximum number of levels exceeded");
    for ( int ln=0; ln<N_levels; ln++ ) {
        tbox::Pointer<hier::PatchLevel > level = d_hierarchy->getPatchLevel(ln);
        // Create the level container
        level_container_array[ln] = new LevelContainer(level->getNumberOfPatches(),d_hierarchy,
            input_data->nvar,u0_id,u_id,input_data->nauxs,auxs_id,input_data->nauxv,auxv_id);
        level_container = (LevelContainer *) level_container_array[ln];
        // Create each patch object
        for (hier::PatchLevel::Iterator p(level); p; p++) {
            tbox::Pointer<hier::Patch> patch = level->getPatch(p());
            level_container->CreatePatch(p(),patch);
        }
    }*/

    // Setup pixie3dRefinePatchStrategy
    d_refine_strategy = new pixie3dRefinePatchStrategy(dim);
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
        const tbox::Pointer<hier::PatchHierarchy> hierarchy = d_hierarchy;
        resetHierarchyConfiguration(hierarchy,0,d_hierarchy->getFinestLevelNumber());
    }

    // Copy the data from u0 to u and set the boundary conditions (needed to set the auxillary variable names)
    d_x->copyVector(d_initial,false);
    synchronizeVariables();

    // Get the variable names
    tbox::Pointer<hier::PatchLevel> level = d_hierarchy->getPatchLevel(0);
    LevelContainer *level_container = (LevelContainer *) level_container_array[0];
    char *tmp_depVarLabels = new char[21*input_data->nvar];
    char *tmp_auxScalarLabels = new char[21*input_data->nauxs];
    char *tmp_auxVectorLabels = new char[21*input_data->nauxv];
    for (hier::PatchLevel::Iterator p(level); p; p++) {
        FORTRAN_NAME(get_var_names)(level_container->getPtr(*p),tmp_depVarLabels,tmp_auxScalarLabels,tmp_auxVectorLabels);
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
void pixie3dApplication::setInitialConditions( const double )
{
    // We have already set u0 we created the level container on the resetHierarchyConfiguration

    // Copy the data from u0 to u 
    d_x->copyVector(d_initial,false);

    // Apply boundary conditions
    synchronizeVariables();
   
    // Form initial conditions
    // Loop through hierarchy
    for ( int ln=0; ln<d_hierarchy->getNumberOfLevels(); ln++ ) {
        tbox::Pointer<hier::PatchLevel> level = d_hierarchy->getPatchLevel(ln);
        // Get the Level container
        LevelContainer *level_container = (LevelContainer *) level_container_array[ln];
        // Loop through the different patches
        for (hier::PatchLevel::Iterator p(level); p; p++) {
            tbox::Pointer<hier::Patch> patch = *p;
            tbox::Pointer< pdat::CellData<double> > tmp = patch->getPatchData(f_src_id);
            double *fsrc = tmp->getPointer();
            int n_elem = patch->getBox().size()*tmp->getDepth();
            // Form Initial Conditions
            #ifdef absoft
                FORTRAN_NAME(FORMINITIALCONDITION)(level_container->getPtr(patch),&n_elem,fsrc);
            #else
                FORTRAN_NAME(forminitialcondition)(level_container->getPtr(patch),&n_elem,fsrc);
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
pixie3dApplication::setInitialConditions( tbox::Pointer< solv::SAMRAIVectorReal<double> > )
{
}

/***********************************************************************
*                                                                      *
* Empty implementation.                                                *
*                                                                      *
***********************************************************************/
void pixie3dApplication::setValuesOnNewLevel( tbox::Pointer<hier::PatchLevel> )
{
}

void
pixie3dApplication::apply(const int *,
              const int *, 
              const int *,
              const int *,
              const int *,
              const int *,
              const double ,
              const double  )
{
}

/***********************************************************************
 *                                                                      *
 * Evaluate right-hand side of IVP being solved.                        *
 *                                                                      *
 ***********************************************************************/
void
pixie3dApplication::apply( tbox::Pointer< solv::SAMRAIVectorReal<double> >  &,
               tbox::Pointer< solv::SAMRAIVectorReal<double> >  &x,
               tbox::Pointer< solv::SAMRAIVectorReal<double> >  &r,
               double, double)
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
    dt_exp = 1e10;
    for (int i=0; i<input_data->nvar; i++)
        f_id[i] = d_x_r->getComponentDescriptorIndex(i);
    // Loop through hierarchy
    for ( int ln=0; ln<d_hierarchy->getNumberOfLevels(); ln++ ) {
        tbox::Pointer<hier::PatchLevel> level = d_hierarchy->getPatchLevel(ln);
        // Get the Level container
        LevelContainer *level_container = (LevelContainer *) level_container_array[ln];
        // Loop through the different patches
        for (hier::PatchLevel::Iterator p(level); p; p++) {
            // Get fsrc
            tbox::Pointer< pdat::CellData<double> > tmp = (*p)->getPatchData(f_src_id);
            double *fsrc = tmp->getPointer();
            int n_elem = (*p)->getBox().size()*tmp->getDepth();
            // Create varray on fortran side
            varrayContainer *varray = new varrayContainer(*p,input_data->nvar,f_id);
            // Create f
            #ifdef absoft
                FORTRAN_NAME(EVALUATENONLINEARRESIDUAL)(level_container->getPtr(*p),&n_elem,fsrc,varray->getPtr());
            #else
                FORTRAN_NAME(evaluatenonlinearresidual)(level_container->getPtr(*p),&n_elem,fsrc,varray->getPtr());
            #endif
            // Comupute the timestep required for an explicit method, computed by pixie3d
            double dt_patch;
            #ifdef absoft
                FORTRAN_NAME(FINDEXPLICITDT)(level_container->getPtr(*p),&n_elem,fsrc,varray->getPtr());
            #else
                  FORTRAN_NAME(findexplicitdt)(level_container->getPtr(*p),&dt_patch);
            #endif
            if ( dt_patch<=0.0 || dt_patch!=dt_patch )
                TBOX_ERROR("Invalid timestep detected\n");
            if ( dt_patch < dt_exp )
                dt_exp = dt_patch;
            // Delete varray
            delete varray;
        }
    }
    // Get the global minimum timestep
    const tbox::SAMRAI_MPI comm = tbox::SAMRAI_MPI::getSAMRAIWorld();
    double dt_tmp = dt_exp;
    comm.Allreduce(&dt_tmp,&dt_exp,1,MPI_DOUBLE,MPI_MIN);

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
void pixie3dApplication::printVector( const tbox::Pointer< solv::SAMRAIVectorReal<double> >)
{
   TBOX_ERROR( "printVector not yet programmed" );
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
    tbox::Pointer< hier::Variable > var0;
    tbox::Pointer< geom::CartesianGridGeometry > grid_geometry = d_hierarchy->getGridGeometry();
    tbox::Pointer< hier::Variable > var;
    tbox::Pointer<hier::PatchLevel > level;
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
    hier::PatchLevel::Iterator p(level);
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
            for (hier::PatchLevel::Iterator ip(level); ip; ip++) {
                pixiePatchData = level_container->getPtr(*ip);
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
    tbox::Pointer< geom::CartesianGridGeometry > grid_geometry = d_hierarchy->getGridGeometry();
    tbox::Pointer< hier::RefineOperator > refine_op;
    if (d_refine_op_str=="CELL_DOUBLE_CUBIC_REFINE") {
        TBOX_ERROR("Not Implimented");
       /*// Create the CartesianCellDoubleCubicRefine operator
       CartesianCellDoubleCubicRefine* temp = new CartesianCellDoubleCubicRefine();
       // manually set the refinement ratio and stencil width for the CartesianCellDoubleCubicRefine
       hier::IntVector width = hier::IntVector(dim,2);
       hier::IntVector ratio = hier::IntVector(dim,3);
       if ( IS2D == 1 ){
          width(2) = 0;
          ratio(2) = 1;
       }
      
       temp->setStencilWidth(width);
       temp->setRefinementRatio(ratio);
       // Add the refinement operator
       refine_op = temp;
       grid_geometry->addSpatialRefineOperator ( refine_op ) ;*/
    } 
   
    // Add coarsen operator
    tbox::Pointer<hier::CoarsenOperator> coarsen_op;
    if (d_coarsen_op_str=="CELL_DOUBLE_INJECTION_COARSEN") {
        TBOX_ERROR("Not Implimented");
        //coarsen_op = new CartesianCellDoubleInjectionCoarsen();
        //grid_geometry->addSpatialCoarsenOperator ( coarsen_op ) ;
    } else if(d_coarsen_op_str=="CELL_DOUBLE_CUBIC_COARSEN") {
        TBOX_ERROR("Not Implimented");
        //coarsen_op = new CartesianCellDoubleCubicCoarsen();
        //grid_geometry->addSpatialCoarsenOperator ( coarsen_op ) ;      
    } 
 
}

int pixie3dApplication::getNumberOfDependentVariables()
{
    //  assert(data!=NULL);
    return input_data->nvar;
}



// Write a single variable to a file for all patches, all levels, including ghostcells
void pixie3dApplication::writeCellData( FILE *fp, int var_id ) {
    if ( dim.getValue() != 3 )
        TBOX_ERROR("Not porgramed for dimensions other than 3");
    const tbox::SAMRAI_MPI comm = tbox::SAMRAI_MPI::getSAMRAIWorld();
    int rank = comm.getRank();
    if ( rank==0 && fp==NULL )
        TBOX_ERROR( "File pointer is NULL" );
    // Loop through the levels
    for ( int ln=0; ln<d_hierarchy->getNumberOfLevels(); ln++ ) {
        tbox::Pointer<hier::PatchLevel > level = d_hierarchy->getPatchLevel(ln);
        const hier::IntVector ratio = level->getRatioToLevelZero();
        // Gather all data to processor 0
        std::vector<commPatchData> patch_data = collectAllPatchData(level,var_id,0);
        // Save the data
        if ( rank==0 ) {
            // Save the header info
            fprintf(fp,"level = %i, ratio = (%i,%i,%i), n_patch = %i\n",ln,ratio(0),ratio(1),ratio(2),level->getNumberOfPatches());
            // Loop through the patches
            for (size_t i=0; i<patch_data.size(); i++) {
                hier::Box box = patch_data[i].getBox();
                hier::IntVector gcw = patch_data[i].getGCW();
                int depth = patch_data[i].getDepth();
                hier::Index ifirst = box.lower();
                hier::Index ilast  = box.upper();
                size_t N = (ilast(0)-ifirst(0)+1+2*gcw(0))*(ilast(1)-ifirst(1)+1+2*gcw(1))*(ilast(2)-ifirst(2)+1+2*gcw(2))*depth;
                // Write the patch
                fprintf(fp,"patch_num = %i, ifirst = (%i,%i,%i), ilast = (%i,%i,%i), gcw = (%i,%i,%i), depth = %i\n",
                    (int) i,ifirst(0),ifirst(1),ifirst(2),ilast(0),ilast(1),ilast(2),gcw(0),gcw(1),gcw(2),depth);
                double *data = patch_data[i].getData();
                fwrite(data,sizeof(double),N,fp);
                fprintf(fp,"\n");
            }
        }
        patch_data.clear();
        comm.Barrier();
    }
}


// Write a single variable to a file using only the coarse level as a single patch
void pixie3dApplication::writeGlobalCellData( FILE *fp, int var_id ) {
    if ( dim.getValue() != 3 )
        TBOX_ERROR("Not programed for dimensions other than 3");
    const tbox::SAMRAI_MPI comm = tbox::SAMRAI_MPI::getSAMRAIWorld();
    int rank = comm.getRank();
    if ( rank==0 && fp==NULL )
        TBOX_ERROR( "File pointer is NULL" );
    // Get the physical domain
    tbox::Pointer<hier::GridGeometry> grid_geometry = d_hierarchy->getGridGeometry();
    if ( grid_geometry->getNumberBlocks() != 1 )
        TBOX_ERROR("Multiblock domains are not supported");
    const SAMRAI::hier::BoxList &physicalDomainList = grid_geometry->getPhysicalDomain(0);
    if ( physicalDomainList.size() != 1 )
        TBOX_ERROR("Multiple box domains are not supported");
    const SAMRAI::hier::Box physicalDomain = physicalDomainList.getBoundingBox();
    const hier::Index lower = physicalDomain.lower();
    const hier::Index upper = physicalDomain.upper();
    int Ngx = upper(0)-lower(0)+1;
    int Ngy = upper(1)-lower(1)+1;
    int Ngz = upper(2)-lower(2)+1;
    // Gather all data to processor 0
    tbox::Pointer<hier::PatchLevel > level = d_hierarchy->getPatchLevel(0);
    std::vector<commPatchData> patch_data = collectAllPatchData(level,var_id,0);
    if ( rank == 0 ) {
        // Create the global patch
        int depth = patch_data[0].getDepth();
        double *data = new double[Ngx*Ngy*Ngz*depth];
        for (size_t i=0; i<patch_data.size(); i++) {
            hier::Box box = patch_data[i].getBox();
            hier::IntVector gcw = patch_data[i].getGCW();
            hier::Index ifirst = box.lower();
            hier::Index ilast  = box.upper();
            int Nx = ilast(0)-ifirst(0)+1+2*gcw(0);
            int Ny = ilast(1)-ifirst(1)+1+2*gcw(1);
            int Nz = ilast(2)-ifirst(2)+1+2*gcw(2);
            // Copy the data to the local patch
            double *tmp_data = patch_data[i].getData();
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
        }
        // Save the results
        fprintf(fp," depth = %i\n",depth);
        fwrite(data,sizeof(double),Ngx*Ngy*Ngz*depth,fp);
        fprintf(fp,"\n");
        delete [] data;
    }
    patch_data.clear();
    // Syncronize the processors
    comm.Barrier();
}


/**************************************************************************
* Write the primary and auxillary variables.                              *
* type = 1:  Write each patch with ghost cells for all levels             *
* type = 2:  Write the coarse level without ghost cells as a single patch *
**************************************************************************/
void pixie3dApplication::writeDebugData( FILE *fp, const int it, const double time, int type ) {
    // Print some information about the time and the domain size
    tbox::Pointer<geom::CartesianGridGeometry> grid_geometry = d_hierarchy->getGridGeometry();
    if ( grid_geometry->getNumberBlocks() != 1 )
        TBOX_ERROR("Multiblock domains are not supported");
    const SAMRAI::hier::BoxList &physicalDomainList = grid_geometry->getPhysicalDomain(0);
    if ( physicalDomainList.size() != 1 )
        TBOX_ERROR("Multiple box domains are not supported");
    const SAMRAI::hier::Box physicalDomain = physicalDomainList.getBoundingBox();
    const hier::Index ilower = physicalDomain.lower();
    const hier::Index iupper = physicalDomain.upper();
    int nbox[3];
    for (int i=0; i<dim.getValue(); i++)
        nbox[i] = iupper(i)-ilower(i)+1;
    const double *lower = grid_geometry->getXLower();
    const double *upper = grid_geometry->getXUpper();
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
* Function overloaded from mesh::StandardTagAndInitStrategy.          *
***********************************************************************/
void pixie3dApplication::initializeLevelData( const tbox::Pointer<hier::PatchHierarchy > hierarchy,
    const int level_number, const double time, const bool can_be_refined, const bool initial_time,
    const tbox::Pointer<hier::PatchLevel > old_level, const bool allocate_data )
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
    tbox::Pointer<hier::PatchLevel> level = hierarchy->getPatchLevel(level_number);

    // Check the new level for overlapping boxes
    const hier::MappedBoxLevel box_level = level->getGlobalizedMappedBoxLevel();
    const hier::MappedBoxSet box_set = box_level.getMappedBoxes();
    hier::PersistentOverlapConnectors persistentOverlap = box_level.getPersistentOverlapConnectors();
    hier::IntVector zero(level->getDim(),0);
    SAMRAI::hier::Connector connector = persistentOverlap.createConnector(box_level,zero);
    hier::NeighborhoodSet neighborhood = connector.getNeighborhoodSets();
    for (std::map<hier::MappedBoxId,hier::MappedBoxSet>::iterator it=neighborhood.begin(); it!=neighborhood.end(); it++) {
        hier::BoxList neighbors;
        it->second.convertToBoxList(neighbors);
        if ( neighbors.size() != 1 )
            TBOX_ERROR("Overlapping boxes were detected on new level, and are not supported");
    }

    // Allocate data when called for, otherwise set timestamp on allocated data.
    hier::ComponentSelector d_problem_data(false);
    if ( !old_level.isNull() ) {
        // Use the old level to determine which components need to be allocated
        // Assume that the old level is a level on a rectangular domain (not a multiblock domain)
        const tbox::Pointer<hier::PatchLevel > old_level2 = old_level;
        for (int i=0; i<d_problem_data.getSize(); i++) {
            if ( old_level2->checkAllocated(i) )
                d_problem_data.setFlag(i);
        }
    } else { 
        // Use the coarsest level to determine which components need to be allocated
        tbox::Pointer<hier::PatchLevel> coarse_level = d_hierarchy->getPatchLevel(0);
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
* Function overloaded from mesh::StandardTagAndInitStrategy.           *
***********************************************************************/
void pixie3dApplication::resetHierarchyConfiguration(
    const tbox::Pointer<hier::PatchHierarchy> hierarchy,
    const int coarsest_level, const int finest_level )
{
    // Check if the application has been initialized
    if ( d_hierarchy.isNull() )
        return;
    TBOX_ASSERT(hierarchy->getDim()==d_hierarchy->getDim());
    int N_levels = d_hierarchy->getNumberOfLevels();
    if ( N_levels > MAX_LEVELS ) 
        TBOX_ERROR("Maximum number of levels exceeded");
    TBOX_ASSERT(coarsest_level>=0);
    TBOX_ASSERT(finest_level<N_levels);

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
    for ( int ln=0; ln<d_hierarchy->getNumberOfLevels(); ln++ ) {
        if ( level_container_array[ln] != NULL ) {
            LevelContainer *level_container = (LevelContainer *) level_container_array[ln];
            delete level_container;
            level_container_array[ln] = NULL;
        }
    }
    for ( int ln=0; ln<N_levels; ln++ ) {
        tbox::Pointer<hier::PatchLevel> level = d_hierarchy->getPatchLevel(ln);
        // Create the level container
        level_container_array[ln] = new LevelContainer(d_hierarchy,level,
            input_data->nvar,u0_id,u_id,input_data->nauxs,auxs_id,input_data->nauxv,auxv_id);
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
    // Get an arbitrary patch to give us the number of boundary sequency groups
    void *pixiePatchData = NULL;
    for ( int ln=0; ln<d_hierarchy->getNumberOfLevels(); ln++ ) {
        LevelContainer *level_container = (LevelContainer *) level_container_array[ln];
        tbox::Pointer<hier::PatchLevel> level = d_hierarchy->getPatchLevel(ln);
        for (hier::PatchLevel::Iterator p(level); p; p++) {
            pixiePatchData = level_container->getPtr(*p);
            break;
        }
        if ( pixiePatchData != NULL )
            break;
    }
    if ( pixiePatchData == NULL )
        TBOX_ERROR("One of the processors does not have any patches");
    // Get the number of boundary condition groups and the sequence for each group
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
    tbox::Pointer<hier::Variable> var0;
    tbox::Pointer<geom::CartesianGridGeometry> grid_geometry = d_hierarchy->getGridGeometry();
    //tbox::Pointer<xfer::PatchLevelFillPattern> fill_pattern(new xfer::PatchLevelBorderFillPattern());
    tbox::Pointer<xfer::PatchLevelFillPattern> fill_pattern(new xfer::PatchLevelFullFillPattern());     // This may reduce performance
    for ( int ln=0; ln<d_hierarchy->getNumberOfLevels(); ln++ ) {
        tbox::Pointer<hier::PatchLevel> level = d_hierarchy->getPatchLevel(ln);
        refineSchedule[ln] = new tbox::Pointer< xfer::RefineSchedule >[d_NumberOfBoundarySequenceGroups];
        siblingSchedule[ln] = new tbox::Pointer< xfer::SiblingGhostSchedule >[d_NumberOfBoundarySequenceGroups];
        for( int iSeq=0; iSeq<d_NumberOfBoundarySequenceGroups; iSeq++) {
            int data_id=-1;
            xfer::RefineAlgorithm refineVariableAlgorithm(d_hierarchy->getDim());
            // Register the variables in the current squence
            std::vector<int> ids(d_BoundarySequenceGroups[iSeq].nbc_seq,-1);
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
                    for(int j=0; j<dim.getValue(); j++) {
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
                    for(int j=0; j<dim.getValue(); j++) {
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
                ids[i] = data_id;
                refineVariableAlgorithm.registerRefine( data_id, data_id, data_id,
                        grid_geometry->lookupRefineOperator(var0,d_refine_op_str) );
            }
            // Create the schedules
            refineSchedule[ln][iSeq] = refineVariableAlgorithm.createSchedule(level, ln-1, d_hierarchy, d_refine_strategy);
            siblingSchedule[ln][iSeq] = tbox::Pointer< xfer::SiblingGhostSchedule >(new xfer::SiblingGhostSchedule(level,ids,ids,ids,fill_pattern));
        }
    }
    // Create the coarsenSchedule
    for ( int ln=0; ln<d_hierarchy->getFinestLevelNumber(); ln++ ) {
        xfer::CoarsenAlgorithm coarsenAlgorithm(d_hierarchy->getDim());
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
        tbox::Pointer<hier::PatchLevel> level = d_hierarchy->getPatchLevel(ln);    
        tbox::Pointer<hier::PatchLevel> flevel = d_hierarchy->getPatchLevel(ln+1);    
        coarsenSchedule[ln] = coarsenAlgorithm.createSchedule(level,flevel);
    }

}


/****************************************************************************
* Collect the data for a given id for all patches onto a single processor   *
****************************************************************************/
std::vector<commPatchData> pixie3dApplication::collectAllPatchData(tbox::Pointer<hier::PatchLevel> level, int id, int root)
{
    // First get the number of patches per processor
    const tbox::SAMRAI_MPI comm = tbox::SAMRAI_MPI::getSAMRAIWorld();
    int rank = comm.getRank();
    int size = comm.getSize();
    int *Npatches = new int[size];
    int NpatchesLocal = level->getLocalNumberOfPatches();
    comm.Allgather(&NpatchesLocal,1,MPI_INT,Npatches,1,MPI_INT);
    int NpatchesGlobal = 0;
    for (int i=0; i<size; i++)
        NpatchesGlobal += Npatches[i];
    TBOX_ASSERT(NpatchesGlobal==level->getGlobalNumberOfPatches());
    // Allocate space for the output
    std::vector<commPatchData> patch_data;
    if ( rank==root )
        patch_data.resize(NpatchesGlobal,commPatchData(level->getDim()));
    // Create the local patch data objects
    std::vector<commPatchData> patch_data_local(NpatchesLocal,commPatchData(level->getDim()));
    int k=0;
    for (hier::PatchLevel::Iterator p(level); p; p++) {
        patch_data_local[k] = commPatchData(*p,id);
        k++;
    }
    // Send all patches to root
    if ( rank==root ) {
        // We are recieving all data
        k=0;
        for (int i=0; i<size; i++) {
            if ( i==rank ) {
                // Perform a local copy
                for (int j=0; j<Npatches[i]; j++) {
                    patch_data[k] = patch_data_local[j];
                    k++;
                }
            } else {
                // Recieve the data from the remote processor
                MPI_Status status;
                comm.Probe(i,132,&status);
                int buffer_size;
                MPI_Get_count(&status,MPI_INT,&buffer_size);
                int *buffer = new int[buffer_size];
                comm.Recv(buffer,buffer_size,MPI_INT,i,132,&status);
                int index = 0;
                for (int j=0; j<Npatches[i]; j++) {
                    patch_data[k].getFromIntBuffer(&buffer[index+1]);
                    index += 1+buffer[index];
                    k++;
                }
                delete [] buffer;
            }
        }
    } else {
        // We are sending our data
        // Compute the necessary buffer size
        size_t buffer_size = NpatchesLocal;
        for (int j=0; j<Npatches[rank]; j++)
            buffer_size += patch_data_local[j].commBufferSize();
        // Create the buffer
        int *buffer = new int[buffer_size];
        // Fill the buffer
        size_t k = 0;
        for (int j=0; j<Npatches[rank]; j++) {
            size_t patch_size = patch_data_local[j].commBufferSize();
            buffer[k] = (int) patch_size;
            patch_data_local[j].putToIntBuffer(&buffer[k+1]);
            k += 1+patch_size;
        }
        // Send the buffer
        comm.Send(buffer,buffer_size,MPI_INT,0,132);
        delete [] buffer;
    }
    // Finished
    comm.Barrier();
    return patch_data;
}


}


