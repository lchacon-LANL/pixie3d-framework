//
// $Id: pixie3dApplication.h 2457 2006-04-18 14:21:27Z philipb $
// $Revision: 2457 $
// $Date: 2006-04-18 08:21:27 -0600 (Tue, 18 Apr 2006) $
//
// File:  pixie3dApplication.h
// Copyright:  (c) 2005 The Regents of the University of California
// Description:  Concrete version of ApplicationStrategy that provides
//               an interface between those interfaces and specific test
//               problems that satisfy the TestProblemStrategy interface.
//

#ifndef included_pixie3dApplication
#define included_pixie3dApplication

#include <string>
#include <vector>

// SAMRAI headers
#include "SAMRAI/tbox/Array.h"
#include "SAMRAI/tbox/Pointer.h"
#include "SAMRAI/pdat/CellVariable.h"
#include "SAMRAI/pdat/FaceVariable.h"
#include "SAMRAI/hier/ComponentSelector.h"
#include "SAMRAI/hier/PatchHierarchy.h"
#include "SAMRAI/xfer/RefinePatchStrategy.h"
#include "SAMRAI/xfer/RefineAlgorithm.h"
#include "SAMRAI/xfer/RefineSchedule.h"
#include "SAMRAI/xfer/CoarsenAlgorithm.h"
#include "SAMRAI/xfer/CoarsenSchedule.h"
#include "SAMRAI/hier/VariableContext.h"
#include "SAMRAI/pdat/FaceData.h"
#include "SAMRAI/appu/VisItDataWriter.h"
#include "SAMRAI/mesh/StandardTagAndInitStrategy.h"

// SAMRUTILS headers
#include "transfer/SiblingGhostSchedule.h"

// SAMRSOLVERS headers
#include "operators/DiscreteOperator.h"

// Local headers
#include "LevelContainer.h"
#include "PatchContainer.h"
#include "CommPatchData.h"
#include "pixie3dApplicationParameters.h"
#include "pixie3dRefinePatchStrategy.h"


//#include "pixie3dData.h"


extern "C"{
/* User-defined application context */
//#include "petscsnes.h"
//#include "petscda.h"
typedef struct {
  double  tolgm;
  double  rtol;
  double  atol;
  double  damp;
  double  dt;
  double  tmax;
  double  mf_eps;
  int        nvar;
  int        nauxs;
  int        nauxv;
  int        ilevel;
  int        nxd;
  int        nyd;
  int        nzd;
  int        npx;
  int        npy;
  int        npz;
  int        numtime;
  int        maxitnwt;
  int        maxksp;
  int        maxitgm;
  int        method;
  int        global;
  int        iguess;
  int        precpass;
  int        sm_flag;
  int        bcs[6];
  int asm_PC;
} input_CTX;
}



namespace SAMRAI{


/** \class pixie3dApplication
 *
 * This is a concrete class that provides an interface between
 * the interface defined by ApplicationStrategy and the test
 * problem for pixie3d.
 */

class pixie3dApplication : 
   public SAMRSolvers::DiscreteOperator,
   public mesh::StandardTagAndInitStrategy
{
public:

   // Empty constructor.
   pixie3dApplication();

   // Constructor that takes a parameter list.  Calls initialize.
   pixie3dApplication( pixie3dApplicationParameters* parameters );

   // Destructor.
   ~pixie3dApplication();

   // Initialize application using specified parameters.
   void initialize( pixie3dApplicationParameters* parameters );
   
   // Set initial conditions on all levels
   void setInitialConditions( const double initial_time );

   // Register location to write initial conditions.
   void setInitialConditions( tbox::Pointer< solv::SAMRAIVectorReal<double> > ic );

   // return the number of dependent variables
   int getNumberOfDependentVariables(void);
   
   // Set data values on new level.
   void setValuesOnNewLevel( tbox::Pointer< hier::PatchLevel> level );

   tbox::Pointer< solv::SAMRAIVectorReal<double> > get_ic() { return(d_x_ic); }

   // Evaluate IVP forcing term.
   void apply( tbox::Pointer< solv::SAMRAIVectorReal<double> >  &f,
	       tbox::Pointer< solv::SAMRAIVectorReal<double> >  &x,
	       tbox::Pointer< solv::SAMRAIVectorReal<double> >  &r,
               double a = -1.0, double b=1.0 );

   // Evaluate IVP forcing term.
   void apply( const int*, const int*, const int*, const int*, const int*, const int*, double a = -1.0, double b=1.0 );

   // Print identifying string.
   void printObjectName( std::ostream& os );

   void registerVizWriter( SAMRAI::appu::VisItDataWriter* visit_writer){d_VizWriter = visit_writer;}

   /*
    * Write the primary and auxillary variables to a binary file.
    * type = 1:  Write each patch with ghost cells for all levels in double precision
    * type = 2:  Write the coarse level without ghost cells as a single patch in double precision
    */
   void writeDebugData( FILE *fp, const int it, const double time, int type=1 );

   // Get the current explicit timestep
   double getExpdT() { return dt_exp; }


/***********************************************************************
* Functions inherited from mesh::StandardTagAndInitStrategy            *
***********************************************************************/

   /*
    * Initialize data on a new level after it is inserted into an AMR
    * patch hierarchy by the gridding algorithm. The level number
    * indicates that of the new level.
    *
    * In this aspect of regridding operations, data that represents
    * the state of the solution must be interpolated from one grid
    * configuration to the next.
    *
    * When assertion checking is enabled, passing in a null pointer
    * for either the hierarchy or the level will result in an
    * unrecoverable exception.
    *
    * Function overloaded from mesh::StandardTagAndInitStrategy.
    *
    */
   void initializeLevelData( const tbox::Pointer<hier::PatchHierarchy> hierarchy,
                             const int level_number,
                             const double time,
                             const bool can_be_refined,
                             const bool initial_time,
                             const tbox::Pointer<hier::PatchLevel> old_level=tbox::Pointer<hier::PatchLevel>(),
                             const bool allocate_data=true );

   /**
    * Reset cached information that depends on the hierarchy configuration.  
    *
    * Function overloaded from mesh::StandardTagAndInitStrategy.
    */
   void resetHierarchyConfiguration(
           const tbox::Pointer<hier::PatchHierarchy> hierarchy,
           const int coarsest_level,
           const int finest_level );

protected:

   void printVector( const tbox::Pointer< solv::SAMRAIVectorReal<double> > vector);

   void generateTransferSchedules( void );

   void coarsenVariables(void);

   void refineVariables(void);

   void synchronizeVariables(void);

   void writeCellData( FILE *fp, int var_id );
   void writeGlobalCellData( FILE *fp, int var_id );

   // Name of application
   std::string d_object_name;

   // Hierarchy
   tbox::Pointer<hier::PatchHierarchy> d_hierarchy;
   tbox::Dimension dim;
   
   tbox::Pointer< solv::SAMRAIVectorReal<double> > d_initial;
   tbox::Pointer< solv::SAMRAIVectorReal<double> > d_x_tmp;
   tbox::Pointer< solv::SAMRAIVectorReal<double> > d_x;
   tbox::Pointer< solv::SAMRAIVectorReal<double> > d_x_r;
   tbox::Pointer< solv::SAMRAIVectorReal<double> > d_x_ic;
   tbox::Pointer< solv::SAMRAIVectorReal<double> > d_aux_scalar;
   tbox::Pointer< solv::SAMRAIVectorReal<double> > d_aux_vector;
   tbox::Pointer< solv::SAMRAIVectorReal<double> > d_aux_scalar_tmp;
   tbox::Pointer< solv::SAMRAIVectorReal<double> > d_aux_vector_tmp;
   tbox::Pointer< pdat::CellVariable<double> > d_f_src;
   int f_src_id;

   bool d_RefineSchedulesGenerated;

   int *u0_id, *u_id, *auxs_id, *auxv_id;
   int *f_id, *u_tmp_id, *auxs_tmp_id, *auxv_tmp_id;

   bool d_bIsInitialTime;

   double dt_exp;   // The current maximum timestep for explicit integration

   pixie3dRefinePatchStrategy* d_refine_strategy; 

   std::string d_refine_op_str;

   //   xfer::SiblingGhostAlgorithm d_siblingGhostVectorAlgorithm;

   tbox::Array< tbox::Pointer<xfer::RefineSchedule> > d_refineScalarSchedules;
   tbox::Array< tbox::Pointer<xfer::RefineSchedule> > d_refineVectorComponentSchedules;
   tbox::Array< tbox::Pointer<xfer::RefineSchedule> > d_refineVectorSchedules;
   tbox::Array< tbox::Pointer<xfer::RefineSchedule> > d_levelSchedules;

   //   tbox::Array< tbox::Pointer<xfer::SiblingGhostSchedule> > d_siblingGhostVectorSchedules;


   SAMRAI::appu::VisItDataWriter* d_VizWriter;

   input_CTX *input_data;

   // Level container info
   static const int MAX_LEVELS = 20;
   LevelContainer *level_container_array[MAX_LEVELS];

   // Data for applying the boundary conditions and the coarsen/refine schedules
   std::string d_coarsen_op_str;
   int d_NumberOfBoundarySequenceGroups;
   pixie3dRefinePatchStrategy::bcgrp_struct *d_BoundarySequenceGroups;
   tbox::Pointer<xfer::RefineSchedule> *refineSchedule[MAX_LEVELS];
   tbox::Pointer<xfer::SiblingGhostSchedule> *siblingSchedule[MAX_LEVELS];
   tbox::Pointer<xfer::CoarsenSchedule> coarsenSchedule[MAX_LEVELS];

   // The names of the primary and auxillary variables
   std::string *depVarLabels;
   std::string *auxScalarLabels;
   std::string *auxVectorLabels;

   // Function to collect the data for a given id for all patches onto a single processor
   // Note: this is provided for writing debug files, requires global communication, and is not scalable
   std::vector<commPatchData> collectAllPatchData(tbox::Pointer<hier::PatchLevel> level, int id, int root);

};
}
#endif
