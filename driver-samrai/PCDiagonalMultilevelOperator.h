#ifndef included_PCDiagonalMultilevelOperator
#define included_PCDiagonalMultilevelOperator


#include "SAMRAI/SAMRAI_config.h"
#include "SAMRAI/tbox/Array.h"
#include "SAMRAI/xfer/CoarsenSchedule.h"

#include <vector>

#include "operators/MultilevelOperatorParameters.h"
#include "operators/MultilevelOperator.h"

namespace SAMRAI {

namespace Pixie3d {

/**\class PCDiagonalMultilevelOperator
 * Class PCDiagonalMultilevelOperator is the multilevel operator that will be used within the preconditioner
 * associated with the density variable
*/

  class PCDiagonalMultilevelOperator: public SAMRSolvers::MultilevelOperator
{
public:
   /**
   * constructor
   * \param parameters
   *        A MultilevelOperatorParameters class used to provide arguments
   *        to initialize the CellDiffusionMultilevelOperator class
   */
   PCDiagonalMultilevelOperator(SAMRSolvers::MultilevelOperatorParameters *parameters);
  
   /**
   * destructor
   */
   ~PCDiagonalMultilevelOperator();

   /**
    * Compute forward apply operation on a level, with the default arguments for a and b
    * the residual, r=f-Au, should be calculated, else b*f+a*A*u is calculated.
    *  \param ln
    *            level number for which apply() should be done
    * \param f_id
    *            array of integer id's for right hand side
    * \param u_id
    *            array of integer id's for solution variables
    * \param r_id
    *            array of integer id's for residual
    * \param f_idx
    *            Array of integer id's for components of right hand side.
    *            These must correspond to the integer id's supplied in f_id
    * \param u_idx
    *            Array of integer id's for components of solution variables.
    *            These must correspond to the integer id's supplied in u_id
    * \param r_idx
    *            Array of integer id's for components of residual.
    *            These must correspond to the integer id's supplied in r_id
    * \param a
    *            Double scaling factor for A*u, default of -1.0
    * \param b
    *            Double scaling factor for rhs, default of 1.0
    */
   void apply(const int ln,
	      const int *f_id,
	      const int *u_id, 
	      const int *r_id,
	      const int *f_idx=NULL,
	      const int *u_idx=NULL,
	      const int *r_idx=NULL,
	      const double a = -1.0,
	      const double b = 1.0);

   /**
    * Compute forward apply operation on a range of levels, with the default arguments for a and b
    * the residual, r=f-Au, should be calculated, else b*f+a*A*u is calculated.
    *  \param coarse_ln
    *            lower bound on levels for which apply() should be done
    *  \param fine_ln
    *            upper bound on levels for which apply() should be done
    * \param f_id
    *            array of integer id's for right hand side
    * \param u_id
    *            array of integer id's for solution variables
    * \param r_id
    *            array of integer id's for residual
    * \param f_idx
    *            Array of integer id's for components of right hand side.
    *            These must correspond to the integer id's supplied in f_id
    * \param u_idx
    *            Array of integer id's for components of solution variables.
    *            These must correspond to the integer id's supplied in u_id
    * \param r_idx
    *            Array of integer id's for components of residual.
    *            These must correspond to the integer id's supplied in r_id
    * \param a
    *            Double scaling factor for A*u, default of -1.0
    * \param b
    *            Double scaling factor for rhs, default of 1.0
    */
   void apply(const int coarse_ln,
	      const int fine_ln,
	      const int *f_id,
	      const int *u_id, 
	      const int *r_id,
	      const int *f_idx=NULL,
	      const int *u_idx=NULL,
	      const int *r_idx=NULL,
	      const double a = -1.0,
	      const double b = 1.0);
   
   /**
   * Apply a boundary condition for a subset of variable, index pairs. Periodic boundaries
   * will be automatically detected and set.
   * \pre Initialization of the boundary condition types for the different variables
   *      should be done when the constructor is called.
   * \param var_id
   *         Variable descriptor index 
   * \param var_idx
   *         Component index
   * \param var_components
   *         We assume that the order of the descriptor indices supplied in var_id and
   *         the corresponding component indices supplied in var_idx need not correspond 
   *         with the order of the variables assumed within the discrete operator. For example
   *         the discrete equations may be ax+by=f, cx+dy=g. Within the operator, x is considered
   *         to be the first variable and y the second variable. However, var_id[0] may correspond
   *         y and var_id[1] to x. In which case var_components[0]=1 and var_components[1]=0
   *         var_components specifies the mapping from (var_id, var_idx) to the order assumed
   *         within the discrete operator.
   * \param number_of_variables
   *         For scalars the default is always true, for systems, in cases where it is desirable
   *         to call applyBoundaryCondition() for a subset of the variables, this specifies the
   *         number of variables for which we are calling applyBoundaryCondition().
   * \param reset_ghost_values
   *          boolean specifying whether interpolation from coarser levels and caching
   *          internally of the data is required before applying boundary conditions, default
   *          is false to minimize interpolations
   */
   void applyBoundaryCondition(const int ln,
			       const int *var_id,
			       const int *var_idx=NULL,
			       const int *var_components=NULL,
			       const int number_of_variables=-1,
			       const bool reset_ghost_values = false);
   
   /**
   * Interpolate source data on a coarser level to destination data on finer level
   * using scratch locations for intermediate calculations. This function is part
   * of the linear operator class because interpolations can require data that is
   * available only to the discrete linear operator, eg, boundary condition data.
   * \param flevel
   *        reference to pointer for fine level
   * \param clevel
   *        reference to pointer for coarse level
   * \param dst_id
   *        pointer to array of descriptor indices for destination data
   * \param src_id
   *        pointer to array of descriptor indices for source data
   * \param scratch_id
   *        pointer to array of descriptor indices for scratch data
   * \param refine_op
   *        name of refine operator to use
   */
   void interpolate(const tbox::Pointer<hier::PatchLevel> &flevel,
		    const tbox::Pointer<hier::PatchLevel> &clevel,
		    const int *dst_id,
		    const int *src_id,
		    const int *scratch_id,
		    std::string &refine_op);

   /**
    * Register class that performs interpolation at coarse/fine boundaries.
    */
   virtual void setRefinementBoundaryInterpolant(SAMRAI::RefinementBoundaryInterpolation *cf_interpolant){d_cf_interpolant=cf_interpolant;}

   /**
   * Returns a pointer to the cached SAMRAI::RefinementBoundaryInterpolation object
   */
   virtual SAMRAI::RefinementBoundaryInterpolation *getRefinementBoundaryInterpolant(void) { return d_cf_interpolant; }

   /**
   * return a pointer to a level linear operator object, can't use a Pointer here thanks to the Pointer semantics being messed up
   */
   SAMRSolvers::LevelOperator *getLevelOperator(const int ln);

   /**
   * This routine controls whether ghost cell values are reset either by interpolation on
   * coarse-fine boundaries or extrapolation on physical boundaries. In certain circumstances
   * data in the ghost cells may already have been reset and a communication step can be
   * avoided by setting this flag to false.
   * \param reset_ghost_values
   *        boolean flag, set to true to reset ghost values
   */
   void resetGhostValues(bool reset_ghost_values){ d_reset_ghost_values=reset_ghost_values; }

   /**
   * Resets the operator internally with new parameters if necessary
   * \param parameters
   *        MultilevelSolverParameters object that is NULL by default
   */
   void reset(SAMRSolvers::DiscreteOperatorParameters *parameters=NULL);

   void initializeBoundaryConditionStrategy(tbox::Pointer<tbox::Database> &db);


   tbox::Pointer< xfer::RefineSchedule > getRefineSchedule(const int ln,
							   const int var_id);

   /**
   * This routine allows a linear operator by default to create
   * new variables that have the same centering as the discrete operator
   * For operators that do not have all variables collocated the option is
   * provided to specify or override the default centering. This routine is
   * useful when a solver needs to create additional variables of its own
   * with the same centering of the linear operator. This routine must be provided
   * by all multilevel linear operators
   * \param name
   *        name for variable
   * \param context
   *        variable context to use
   * \param variable
   *        reference to variable (output)
   * \param nghosts
   *        number of ghosts - applicable only when a new variable is created
   * \param depth
   *        depth of variable
   * \param bOverride
   *        boolean to override default centering provided by operator
   * \param centering
   *        string specifying what centering to use if overriding default centering
   *        values are CELL, NODE, SIDE, FACE, EDGE, OUTERFACE, OUTERSIDE, OUTERNODE
   *        and the interlaced versions CCELL, CNODE etc provided by MLUtilities
   */
   int getVariableIndex(std::string &name, 
			tbox::Pointer<hier::VariableContext> &context,
			tbox::Pointer<hier::Variable> &variable,
			hier::IntVector nghosts,
			int depth = 1,
			bool bOverride = false,
			std::string centering = "");

   /**
    * Get information describing physical boundary conditions.  The
    * return value is a pointer to 2*d integers that describe the
    * physical boundary conditions (d is the spatial dimension).
    */
   const int* getPhysicalBoundaryConditions(void) const{ return d_bdry_types; }

   
protected:
   
   /**
   * Allocates space for the internally used face centered flux
   */
   void initializeInternalVariableData(void);

   void getFromInput(tbox::Pointer<tbox::Database> db);


   /**
   * Coarsen a souce term like variable from level ln+1 to level ln
   * currently this routine is identical to coarsenVariable but that may
   * change
   * \param ln
   *        level number
   * \param u_id
   *        descriptor index of solution to coarsen
   * \param src_id
   *        descriptor index of rhs to coarsen
   * \param b
   *        boolean to determine if b should be coarsened or not
   */
   void coarsenSolutionAndSourceTerm(const int ln, 
                                     const int u_id,
                                     const int src_id,
                                     const bool coarsen_rhs);

   /**
   * Compute fluxes.
   *  \param coarse_ln
   *            lower bound on levels for which setFlux() should be called
   *  \param fine_ln
   *            upper bound on levels for which setFlux() should be called
   * \param u_id
   *        descriptor index for cell centered variable to calculate flux from
   * \param u_idx
   *        component index for cell centered variable to calculate flux from
   */

   void setFlux(const int coarse_ln,
                const int fine_ln,
                const int *u_id,
                const int *u_idx=NULL);

   /**
   * Routine used to set up internal transfer schedules
   */
   void setupTransferSchedules(void);

 private:
   /**
    * Default constructor, should not be used in general
    */
   PCDiagonalMultilevelOperator();

   void initializeLevelOperators(SAMRSolvers::MultilevelOperatorParameters *parameters);

   int d_flux_id;

   bool d_coarsen_diffusive_fluxes;
   bool d_schedules_initialized;
   bool d_use_cf_interpolant;
   bool d_variable_order_interpolation;
   bool d_reset_ghost_values;

   std::string d_face_coarsen_op_str;
   std::string d_face_refine_op_str;
   std::string d_cell_coarsen_op_str;
   std::string d_cell_refine_op_str;
   std::string d_tangent_interp_scheme_str;
   std::string d_normal_interp_scheme_str;

   int *d_bdry_types;

   tbox::Pointer< pdat::FaceVariable<double> > d_flux;

   tbox::Array< tbox::Pointer<xfer::CoarsenSchedule > > d_flux_coarsen_schedule;
   tbox::Array< tbox::Pointer<xfer::CoarsenSchedule > > d_src_coarsen_schedule;

   tbox::Array< tbox::Pointer<xfer::RefineSchedule > > d_var_refine_schedule;
   tbox::Array< tbox::Pointer<xfer::RefineSchedule > > d_interpolate_schedule;

   SAMRAI::RefinementBoundaryInterpolation::InterpolationScheme d_tangent_interp_scheme;
   SAMRAI::RefinementBoundaryInterpolation::InterpolationScheme d_normal_interp_scheme;

   tbox::Array< tbox::Pointer<SAMRSolvers::LevelOperator> > d_level_operators;

};
  
}
}

#endif
