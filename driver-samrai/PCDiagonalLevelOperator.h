#ifndef included_PCDiagonalLevelOperator
#define included_PCDiagonalLevelOperator


#include "SAMRAI/SAMRAI_config.h"
#include "SAMRAI/hier/PatchLevel.h"
#include "SAMRAI/tbox/Pointer.h"
#include "SAMRAI/xfer/RefineSchedule.h"

#ifdef DEBUG_CHECK_ASSERTIONS
extern "C"{
#include "assert.h"
}
#endif

#include "SAMRAI/pdat/FaceVariable.h"
#include "patchdata/CCellVariable.h"
#include "PCDiagonalRefinePatchStrategy.h"

#include "operators/LevelOperatorParameters.h"
#include "operators/LevelOperator.h"

namespace SAMRAI {

namespace Pixie3d {

/**\class PCDiagonalLevelOperator
 * 
 * Class PCDiagonalLevelOperator represents linear operators on a cell-centered array on an SAMRAI patch level.
 * The discretization is the standard second order finite difference stencil. Instances
 * of this class are typically created internally by the PCDiagonalMultilevelOperator.
 * However, it may be used independent of that class.
 * 
 * The user must perform the following steps to use the PCDiagonalLevelOperator
 * 
 * -# Create a  PCDiagonalLevelOperator object using a 
 *    LevelOperatorParameters class as argument to the constructor.
 * -# LevelOperatorParameters object must contain a database object with the following keys and values.
 *    Those that have defaults do not need to be specified.    
 *      - alpha                       : double constant with default value 0.0
 *      - beta                        : double constant with default value 1.0
 *      - a_index                     : integer, descriptor index for cell centered array of 'a' values,\n
 *                                      default is -1, in which a is assumed to be all zeros
 *      - b_index                     : integer, descriptor index for face centered array of 'b' values,\n
 *                                      default is -1, in which b is assumed to be all ones
 *      - tangent_interp_scheme       : std::string, tangential interpolation to use at coarse-fine interfaces\n
 *                                      typical values are "LINEAR, ""QUADRATIC", "CUBIC" etc
 *      - normal_interp_scheme        : std::string, normal interpolation to use at coarse-fine interfaces\n
 *                                      typical values are "LINEAR, ""QUADRATIC", "CUBIC" etc
 *      - adjust_cf_coefficients      : boolean, set to TRUE only for simple grid configurations with no concave corners\n
 *                                      determines whether stencil coefficients at course fine interfaces should be\n
 *                                      adjusted to account for non uniform grid spacing at cf interfaces.
 *      - interpolate_ghost_values    : boolean
 *      - extrapolation_order         : integer, specifies order of extrapolation to use at physical boundaries, 0 does nothing
 *      - boundary_conditions         : integer array of size 2*NDIM specifies bc's on all sides.\n
 *                                      DIRICHLET=0, NEUMANN=1,MIXED=2,PERIODIC=3,ROBIN=4
 *      - use_cf_interpolant          : boolean, determines whether the SAMRAI::RefinementBoundaryInterpolation class is used or not
 *      - variable_order_interpolation: boolean, if TRUE lowers order of interpolation where it is hard to construct interpolants 
 *      - print_info_level            : integer, controls verbosity of deug information printed, higher is more
 *      - weight_id                   : integer, descriptor index for weighting in computing norms etc, default is -1
 */

class PCDiagonalLevelOperator: public SAMRSolvers::LevelOperator
{
public:
   PCDiagonalLevelOperator(SAMRSolvers::LevelOperatorParameters *parameters);

   ~PCDiagonalLevelOperator();

   /**
    * Compute forward apply operation, with the default arguments for a and b
    * the residual, r=f-Au, should be calculated, else b*f+a*A*u is calculated.
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
   void apply(const int *f_id,
              const int *u_id, 
              const int *r_id,
              const int *f_idx=NULL,
              const int *u_idx=NULL,
              const int *r_idx=NULL,
              const double a = -1.0,
              const double b = 1.0);


   /**
   * Apply a boundary condition for a specific variable and index. Periodic boundaries
   * will be automatically detected and set.
   * \pre setBoundaryConditionTypes must have been already called
   *      to initialize the boundary condition types for the different variables or initialization
   *      has been done when the constructor was called.
   * \param var_id
   *         variable descriptor index 
   * \param var_idx
   *         component index
   */
   void applyBoundaryCondition(const int *var_id,
                               const int *var_idx=NULL,
                               const int *var_components=NULL,
                               const int number_of_variables=-1);
   
   /**
   * Compute fluxes.
   * \pre This routine assumes data in ghost cells has been constant refined
   *      from coarser grids on c-f boundaries and extrapolated on physical
   *      boundaries. Using the cf_interpolant data is realigned before calculating
   *      fluxes.
   * \param u_id
   *        descriptor index for cell centered variable to calculate flux from
   * \param u_idx
   *        component index for cell centered variable to calculate flux from
   * \param flux_id
   *        descriptor index for face centered variable to store flux
   */   
   void setFlux(const int flux_id,
                const int *u_id,
                const int *u_idx=NULL);

   /**
   * Set extrapolation order to use on physical domain boundaries
   * \param extrapolation_order
   *        Order of extrapolation, currently should only be set to 1 or 2.
   * \param var_component
   *         variable component default set to 0
   */
   void setExtrapolationOrder(const int extrapolation_order);

   /**
   * Specify interpolation schemes in specified direction relative
   * to coarse/fine interface.
   * \param scheme
   *        interpolation scheme to use normal to a coarse fine interface
   * \param dir
   *        dir=0, is tangential and dir=1, is normal to interface
   * \param var_component
   *         variable component default set to 0
   */
   void setInterpolationScheme(SAMRAI::RefinementBoundaryInterpolation::InterpolationScheme scheme,
                               const int dir);
   
   /**
   * Return the number of non zero stencil elements for block i,j of the stencil.
   * Block i,j of the stencil contains non zero connections between the i-th and j-th variables.
   * \param i
   *        i is the first variable
   * \param j
   *        index of second variable
   */
   int getStencilSize(const int i=0,
		      const int j=0,
		      const int k=0){return 0;}
   
   /**
   * Returns a pointer to the stencil block data for block i,j,k of the stencil.
   * Block i,j,k of the stencil contains non zero connections between the i-th, j-th, and k-th variables.
   * Will return a NULL pointer if all elements of that stencil block are zero.
   * \param p
   *        patch number for which to get stencil data
   * \param i
   *        index of first variable
   * \param j
   *        index of second variable
   * \param k
   *        index of third variable
   */
   tbox::Pointer< hier::PatchData > getStencilBlock(tbox::Pointer<hier::Patch> p, 
						    const int i=0, 
						    const int j=0,
						    const int k=0);

   /**
   * Returns integer array of offsets for non zero elements of the i,j stencil block.
   * Block i,j of the stencil contains non zero connections between the i-th and j-th variables.
   * \param i
   *        index of first variable
   * \param j
   *        index of second variable   
   * \param k
   *        index of third variable   
   */
   std::vector<int> getStencilOffsets(const int i=0, const int j=0, const int k=0);

   /**
    * Compute b*f+a*A*u from fluxes, default values for a and b yield the residual
    */
   void apply(const int flux_id, 
              const int *f_id,
              const int *u_id,
              const int *r_id,
              const int *f_idx=NULL,
              const int *u_idx=NULL,
              const int *r_idx=NULL,
              const double a=-1, 
              const double b=1.0);

   /**
   * Returns a pointer to a LevelOperator. The purpose of this function
   * is to create a LevelOperator for a PatchLevel that is a subset in index space of
   * the PatchLevel on which the existing LevelOperator is constructed. It allows for
   * the construction of a LevelOperator for restricted PatchLevels used in algorithms
   * such as AFAC and AFACx.
   * \param parameters
   *        LevelOperatorParameters object that should contain a pointer to a valid PatchLevel object
   */
   SAMRSolvers::LevelOperator *constructOperator(SAMRSolvers::LevelOperatorParameters *parameters);

   /**
    * Get information describing physical boundary conditions.  The
    * return value is a pointer to 2*d integers that describe the
    * physical boundary conditions (d is the spatial dimension).
    */
   const int* getPhysicalBoundaryConditions(void) const{ return d_bdry_types; }

   int getVariableIndex(std::string &name, 
                        tbox::Pointer<hier::VariableContext> &context,
                        tbox::Pointer<hier::Variable > &var,
                        hier::IntVector nghosts,
                        int depth=1,
                        bool bOverride=false,
                        std::string centering="");

      
   void reset(SAMRSolvers::DiscreteOperatorParameters *params);
   
protected:
   
   void getFromInput(const tbox::Pointer<tbox::Database> &db);

   void initializeInternalVariableData(void);

   /**
   * Internal function used to populate entries of a database object which
   * will then be used in creating a  LevelOperator
   */
   void initializeDatabase(tbox::Pointer<tbox::Database> db);

private:
   PCDiagonalLevelOperator();

   bool d_interpolate_ghost_values;
   bool d_sibling_fill_cached;
   bool d_coefficients_changed;
   bool d_variable_order_interpolation;
   bool d_reset_ghost_cells;
   bool d_use_cf_interpolant;
   
   int d_extrapolation_order;           // extrapolation order to use for ghost cells on physical boundaries
   int *d_bdry_types;

   int d_flux_id;

   tbox::Pointer<xfer::RefineSchedule > d_sibling_fill_schedule;

   tbox::Pointer<pdat::FaceVariable<double> > d_flux;

   SAMRAI::RefinementBoundaryInterpolation::InterpolationScheme d_tangent_interp_scheme;
   SAMRAI::RefinementBoundaryInterpolation::InterpolationScheme d_normal_interp_scheme;

};
}
}


#endif
