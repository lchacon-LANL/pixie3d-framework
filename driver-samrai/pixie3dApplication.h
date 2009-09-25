//
// $Id: pixie3dApplication.h 2457 2006-04-18 14:21:27Z pernice $
// $Revision: 2457 $
// $Date: 2006-04-18 08:21:27 -0600 (Tue, 18 Apr 2006) $
//
// File:  pixie3dApplication.h
// Copyright:  (c) 2005 The Regents of the University of California
// Description:  Concrete version of ApplicationStrategy that provides
//               an interface between those interfaces and specific test
//               problems that satisfy the TestProblemStrategy interface.
//

#ifndef included_pixie3d_application
#define included_pixie3d_application

#include <string>
#include <vector>

#include "tbox/Array.h"
#include "tbox/Pointer.h"
#include "CellVariable.h"
#include "FaceVariable.h"
#include "ComponentSelector.h"
#include "PatchHierarchy.h"
#include "RefinePatchStrategy.h"
#include "RefineAlgorithm.h"
#include "RefineSchedule.h"
#include "CoarsenAlgorithm.h"
#include "CoarsenSchedule.h"
#include "VariableContext.h"
#include "FaceData.h"
#include "LevelContainer.h"

#include "ApplicationStrategy.h"
//#include "TagAndInitDataServer.h"
#include "pixie3dApplicationParameters.h"
// #include "pixie3dProblemStrategy.h"

#ifndef LACKS_NAMESPACE
using namespace SAMRAI;
#endif

/* User-defined application context */
#include "petscsnes.h"
#include "petscda.h"
typedef struct {
  PetscReal  tolgm;
  PetscReal  rtol;
  PetscReal  atol;
  PetscReal  damp;
  PetscReal  dt;
  PetscReal  tmax;
  PetscReal  mf_eps;
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
  PetscTruth asm_PC;
} input_CTX;


/** \class pixie3dApplication
 *
 * This is a concrete class that provides an interface between
 * the interface defined by ApplicationStrategy and the test
 * problem for pixie3d.
 */

// **** Remove xfer::RefinePatchStrategy<NDIM> ****
class pixie3dApplication : public ApplicationStrategy
{
public:
   // Constructor that takes a parameter list.  Calls initialize.
   pixie3dApplication( const pixie3dApplicationParameters* parameters );

   // Destructor.
   ~pixie3dApplication();

   // Initialize application using specified parameters.
   void initialize( const pixie3dApplicationParameters* parameters );
   
   // Set initial conditions on all levels
   void setInitialConditions( const double initial_time );

   // Set initial conditions on new level.
   void setInitialConditions( const double initial_time, tbox::Pointer< hier::PatchLevel<NDIM> > level );

   // Register location to write initial conditions.
   void setInitialConditions( tbox::Pointer< solv::SAMRAIVectorReal<NDIM,double> > ic );

   // Return ComponentSelector of data to allocate on a new level.
   hier::ComponentSelector getDataToAllocate();

   // Return ComponentSelector of data to time stamp on a new level.
   hier::ComponentSelector getDataToTimeStamp();

   // Set data values on new level.
   void setValuesOnNewLevel( tbox::Pointer< hier::PatchLevel<NDIM> > level );

   // Return a list of variables used by this application.
   tbox::Array< tbox::Pointer< hier::Variable<NDIM> > > getVariables();
   tbox::Pointer< solv::SAMRAIVectorReal<NDIM,double> > get_x() { return(d_x); }

   // Create the RefineSchedules needed to transfer data to a new level.
   tbox::Array< tbox::Pointer< xfer::RefineSchedule<NDIM> > > setupRegridRefineSchedules( 
      const tbox::Pointer< hier::BasePatchHierarchy<NDIM> > hierarchy, 
      const int level_number,  
      const tbox::Pointer< hier::BasePatchLevel<NDIM> > old_level );

   // Mark locations where additional refinement is desired.
   void tagCells( const tbox::Pointer< hier::PatchHierarchy<NDIM> > hierarchy,
                          const int level_number,
                          const double error_data_time,
                          const bool initial_time,
                          const int tag_index,
                          const bool uses_richardson_extrapolation_too );

   // Update data structures that change when the grid hierarchy changes.
   void resetHierarchyConfiguration( const tbox::Pointer< hier::PatchHierarchy<NDIM> > hierarchy,
                                             const int coarsest_level,
                                             const int finest_level );

   // Evaluate IVP forcing term.
   int evaluateFunction( tbox::Pointer< solv::SAMRAIVectorReal<NDIM,double> >  x,
                         tbox::Pointer< solv::SAMRAIVectorReal<NDIM,double> >  f );

   // Evaluate solution.
   void evaluateSolution( tbox::Pointer< solv::SAMRAIVectorReal<NDIM,double> >  y );

   // Get next timestep
   double getNextDt();

   // Set boundary conditions.
   //void setPhysicalBoundaryConditions( hier::Patch<NDIM>& patch, const double time, const hier::IntVector<NDIM>& ghost_width_to_fill);

   // Print identifying string.
   void printObjectName( ostream& os );

   // Allocate data 
   void allocateVectorData(SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM,double> > x, double time, bool flag=true ); 

private:
   
   void printVector( const tbox::Pointer< solv::SAMRAIVectorReal<NDIM,double> > vector);

   // Empty constructor.
   pixie3dApplication();

   void generateTransferSchedules( void );

   void coarsenVariables(void);

   void refineVariables(void);

   void synchronizeVariables(void);

   // Name of application
   string d_object_name;

   // Hierarchy
   tbox::Pointer< hier::PatchHierarchy<NDIM> > d_hierarchy;

   // Data Variables
   PetscReal time;
   double d_initial_time;
   tbox::Pointer< solv::SAMRAIVectorReal<NDIM,double> > d_equilibrium;
   tbox::Pointer< solv::SAMRAIVectorReal<NDIM,double> > d_initial;
   tbox::Pointer< solv::SAMRAIVectorReal<NDIM,double> > d_x_tmp;
   tbox::Pointer< solv::SAMRAIVectorReal<NDIM,double> > d_x;
   tbox::Pointer< solv::SAMRAIVectorReal<NDIM,double> > d_f;
   tbox::Pointer< solv::SAMRAIVectorReal<NDIM,double> > d_aux_scalar;
   tbox::Pointer< solv::SAMRAIVectorReal<NDIM,double> > d_aux_vector;
   tbox::Pointer< solv::SAMRAIVectorReal<NDIM,double> > d_aux_scalar_tmp;
   tbox::Pointer< solv::SAMRAIVectorReal<NDIM,double> > d_aux_vector_tmp;
   tbox::Pointer< pdat::CellVariable<NDIM,double> > d_f_src;
   int f_src_id;
   int d_number_solution_components;
   tbox::Array< tbox::Pointer< hier::Variable<NDIM> > > d_variable_list;
   void **level_container_array;

   int nbc_seq, *bc_seq;
   int u0_id[NVAR], u_id[NVAR], auxs_id[NAUXS], auxv_id[NAUXV];
   int f_id[NVAR], u_tmp_id[NVAR], auxs_tmp_id[NAUXS], auxv_tmp_id[NAUXV];

   // Misc. variables   
   hier::ComponentSelector d_application_data;
   tbox::Array< tbox::Pointer< xfer::RefineSchedule<NDIM> > > d_regrid_refine_scheds;
   tbox::Pointer<hier::VariableContext> d_application_ctx;
   input_CTX data;

   // Coarsen operator   
   //tbox::Pointer<xfer::CoarsenOperator<NDIM> > x_coarsen_op;

   // Boundary conditions
   //xfer::RefinePatchStrategy<NDIM>* bcfill;

   // Refinement operator
   //tbox::Pointer<xfer::RefineOperator<NDIM> > x_refine_op;

   string d_refine_op_str;
   xfer::RefineAlgorithm<NDIM> d_refine_algorithm;
   tbox::Array< tbox::Pointer<xfer::RefineSchedule<NDIM> > > d_refine_schedules;
   tbox::Array< tbox::Pointer<xfer::RefineSchedule<NDIM> > > d_level_schedules;
   xfer::RefinePatchStrategy<NDIM> *d_refine_strategy; 

   string d_coarsen_op_str;
   xfer::CoarsenAlgorithm<NDIM> d_cell_coarsen_alg;
   tbox::Array< tbox::Pointer<xfer::CoarsenSchedule<NDIM> > > d_cell_coarsen_schedules;

};
#endif
