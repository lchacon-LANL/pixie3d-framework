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
#include "boost/shared_ptr.hpp"
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
#include "transfer/TriangleRefineSchedule.h"
#include "interpolation/RefinementBoundaryInterpolation.h"

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
 *
 * The input database specifies several key parameters:
 * @code
       refine_method - Specifies the method used for refinement interpolation
                       "CONSTANT" - Use constant interpolation from the coarse data
                       "COARSE_LINEAR" - Use linear interpolation using only coarse data
                       "FINE_LINEAR" - Use linear interpolation using the coarse and fine data
       coarsen_method - Specifies the method use for coarsening
                       "CONSERVATIVE" - Use conservative coarsening for all quantities
       print_info_level - Specifies the level of printing
 * @endcode
 */

class pixie3dApplication : 
   public SAMRSolvers::DiscreteOperator,
   public mesh::StandardTagAndInitStrategy
{
public:

   // Empty constructor.
   pixie3dApplication();

   // Constructor that takes a parameter list.  Calls initialize.
   pixie3dApplication( boost::shared_ptr<pixie3dApplicationParameters> parameters );

   // Destructor.
   virtual ~pixie3dApplication();

   // Initialize application using specified parameters.
   void initialize( boost::shared_ptr<pixie3dApplicationParameters> parameters );
   
   // Set initial conditions on all levels
   void setInitialConditions( const double initial_time );

   // Register location to write initial conditions.
   void setInitialConditions( boost::shared_ptr< solv::SAMRAIVectorReal<double> > ic );

   // return the number of dependent variables
   int getNumberOfDependentVariables(void);
   
   // Set data values on new level.
   void setValuesOnNewLevel( boost::shared_ptr< hier::PatchLevel> level );

   boost::shared_ptr< solv::SAMRAIVectorReal<double> > get_ic() { return(d_x_ic); }

   // Evaluate IVP forcing term.
   void apply( boost::shared_ptr< solv::SAMRAIVectorReal<double> >  &f,
           boost::shared_ptr< solv::SAMRAIVectorReal<double> >  &x,
           boost::shared_ptr< solv::SAMRAIVectorReal<double> >  &r,
               double a = -1.0, double b=1.0 );

   // Evaluate IVP forcing term.
   void apply( const int*, const int*, const int*, const int*, const int*, const int*, double a = -1.0, double b=1.0 );

   void registerVector( boost::shared_ptr< solv::SAMRAIVectorReal<double> > x );

   // Print identifying string.
   void printObjectName( std::ostream& os );

   void registerVizWriter( boost::shared_ptr<SAMRAI::appu::VisItDataWriter> visit_writer){d_VizWriter = visit_writer;}

   /*
    * Write the primary and auxillary variables to a binary file.
    * type = 1:  Write each patch with ghost cells for all levels using 14 digits of precision
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
   void initializeLevelData( const boost::shared_ptr<hier::PatchHierarchy>& hierarchy,
                             const int level_number,
                             const double time,
                             const bool can_be_refined,
                             const bool initial_time,
                             const boost::shared_ptr<hier::PatchLevel>& old_level=boost::shared_ptr<hier::PatchLevel>(),
                             const bool allocate_data=true );

   /**
    * Reset cached information that depends on the hierarchy configuration.  
    *
    * Function overloaded from mesh::StandardTagAndInitStrategy.
    */
   virtual void resetHierarchyConfiguration( const boost::shared_ptr<hier::PatchHierarchy>& hierarchy,
                         const int coarsest_level,
                         const int finest_level );

   /**
    * Tag cells where the gradient of the solution exceeds a user-specified
    * threshold.  
    *
    * Function overloaded from mesh::StandardTagAndInitStrategy.
    *
    */
   virtual void applyGradientDetector( const boost::shared_ptr<hier::PatchHierarchy>& hierarchy,
			    const int level_number,
			    const double time,
			    const int tag_index,
			    const bool initial_time,
			    const bool uses_richardson_extrapolation_too );

protected:

   void printVector( const boost::shared_ptr< solv::SAMRAIVectorReal<double> > vector);

   void generateTransferSchedules( void );

   void coarsenVariables(void);

   void refineVariables(void);

   void synchronizeVariables(void);

   void writeCellData( FILE *fp, int var_id );
   void writeGlobalCellData( FILE *fp, int var_id );

   // Name of application
   std::string d_object_name;

   // Hierarchy
   boost::shared_ptr<hier::PatchHierarchy> d_hierarchy;
   tbox::Dimension dim;
   
   // Solution vectors
   boost::shared_ptr< solv::SAMRAIVectorReal<double> > d_initial;
   boost::shared_ptr< solv::SAMRAIVectorReal<double> > d_x_tmp;
   boost::shared_ptr< solv::SAMRAIVectorReal<double> > d_x;
   boost::shared_ptr< solv::SAMRAIVectorReal<double> > d_x_r;
   boost::shared_ptr< solv::SAMRAIVectorReal<double> > d_x_ic;

   // A list of all registered vectors that need to be interpolated on a regrid
   std::vector<boost::shared_ptr< solv::SAMRAIVectorReal<double> > >  d_registeredVectors;

   // Auxillary vectors
   boost::shared_ptr< solv::SAMRAIVectorReal<double> > d_aux_scalar;
   boost::shared_ptr< solv::SAMRAIVectorReal<double> > d_aux_vector;
   boost::shared_ptr< solv::SAMRAIVectorReal<double> > d_aux_scalar_tmp;
   boost::shared_ptr< solv::SAMRAIVectorReal<double> > d_aux_vector_tmp;
   boost::shared_ptr< pdat::CellVariable<double> > d_f_src;
   boost::shared_ptr< pdat::CellVariable<double> > d_div_B;

   // Data used for storage for the gradient
   std::vector<int> d_x_grad_ids;       // storage for d_x
   std::vector<int> d_x0_grad_ids;      // storage for x (without ghosts)
   std::vector<int> d_auxs_grad_ids;    // storage for d_aux_scalar
   std::vector<int> d_auxv_grad_ids;    // storage for d_aux_vector

    // Other data
   int f_src_id;
   int div_B_id;

   int d_weight_id;
   
   bool d_RefineSchedulesGenerated;

   bool d_vectorsCloned;
   
   int *u0_id, *u_id, *auxs_id, *auxv_id, flux_id, flux_src_id;
   int *f_id, *u_tmp_id, *auxs_tmp_id, *auxv_tmp_id;

   hier::ComponentSelector d_problem_data;

   bool d_bIsInitialTime;

   double dt_exp;   // The current maximum timestep for explicit integration

   boost::shared_ptr<SAMRAI::appu::VisItDataWriter> d_VizWriter;

   input_CTX *input_data;

   // Level container info
   static const int MAX_LEVELS = 20;
   LevelContainer *level_container_array[MAX_LEVELS];

   // Data for applying the boundary conditions and the refine schedules
   std::string d_refine_op_str;
   boost::shared_ptr<pixie3dRefinePatchStrategy> d_refine_strategy; 
   std::vector<pixie3dRefinePatchStrategy::bcgrp_struct> d_BoundarySequenceGroups;
   boost::shared_ptr<xfer::TriangleRefineSchedule> *refineSchedule[MAX_LEVELS];

   // Data for the coarsen schedules
   std::string d_coarsen_op_str;
   boost::shared_ptr<xfer::CoarsenSchedule> coarsenSchedule[MAX_LEVELS];
   boost::shared_ptr<RefinementBoundaryInterpolation> d_coarseFineInterp;
    
   // Data for regrid
   std::vector<double> d_J_level;
   std::string d_regrid_op_str;

   // The names of the primary and auxillary variables
   std::string *depVarLabels;
   std::string *auxScalarLabels;
   std::string *auxVectorLabels;

   // Function to collect the data for a given id for all patches onto a single processor
   // Note: this is provided for writing debug files, requires global communication, and is not scalable
   std::vector<commPatchData> collectAllPatchData(boost::shared_ptr<hier::PatchLevel> level, int id, int root);

   // Get the boundary condition groups
   std::vector<pixie3dRefinePatchStrategy::bcgrp_struct> getBCgroup();
   
   // Function to check for overlapping boxes
   static bool overlappingBoxes( const boost::shared_ptr<hier::PatchLevel> level );

};
}
#endif
