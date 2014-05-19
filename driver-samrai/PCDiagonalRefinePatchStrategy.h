
#ifndef included_PCDiagonalRefinePatchStrategy
#define included_PCDiagonalRefinePatchStrategy

#include "SAMRAI/SAMRAI_config.h"
#include "SAMRAI/xfer/RefinePatchStrategy.h"
#include "SAMRAI/hier/PatchHierarchy.h"
#include "SAMRAI/hier/Patch.h"
#include "boost/shared_ptr.hpp"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/Box.h"
#include "SAMRAI/tbox/Dimension.h"

#include "boundaries/BoundaryConditionParameters.h"

namespace SAMRAI {

namespace Pixie3d {

/** \class PCDiagonalRefinePatchStrategy
 *
 * Class PCDiagonalRefinePatchStrategy provides a concrete
 * xfer::RefinePatchStrategy.  It
 * is used to fill ghost cells that represent physical boundary
 * conditions, which, for a preconditioner, are merely homogeneous
 * analogs of the actual problem being solved.  Other mechanisms are
 * needed to ensure proper residual updates, in particular the
 * computation of fluxes at physical boundaries.
 *
 * @see xfer::RefinePatchStrategy
 */

class PCDiagonalRefinePatchStrategy :
   public xfer::RefinePatchStrategy
{
public:
   /**
    * Blank constructor for PCDiagonalRefinePatchStrategy.
    */
   PCDiagonalRefinePatchStrategy( boost::shared_ptr<BoundaryConditionParameters> parameters);

   /**
    * Virtual destructor for PCDiagonalRefinePatchStrategy
    */

   virtual ~PCDiagonalRefinePatchStrategy();

   /**
    * Set boundary types.
    */
   void setBoundaryTypes(int* bdry_types);

   void setExtrapolationOrder(int order);

   /**
    * Re-do extrapolation of corner ghost cells after interpolation.
    */
   void extrapolateCornerGhostCells(
      boost::shared_ptr<hier::PatchHierarchy > hierarchy, 
      const int ln,
      const int var_id,
      const int var_idx=0);

   /**
    * Set solution ghost cell values along physical boundaries.
    *
    * Function is overloaded from xfer::RefinePatchStrategy.
    */
   void setPhysicalBoundaryConditions(hier::Patch& patch,
                                  const double time,
                                  const hier::IntVector& ghost_width_to_fill);

   /**
    * Set descriptor index for storage location where ghost cells 
    * need to be filled. 
    */
   void setIndexToFill(const int id, const int idx=0) 
   {
      d_data_id = id;
      d_data_idx=idx;
   }

   /**
    * Function for applying user-defined data processing prior
    * to refinement.  
    *
    * Function is overloaded from xfer::RefinePatchStrategy.  An empty
    * implementation is provided here.
    */
   void preprocessRefine(hier::Patch& fine,
                         const hier::Patch& coarse,
                         const hier::Box& fine_box,
                         const hier::IntVector& ratio)
   {
     (void) fine;
     (void) coarse;
     (void) fine_box;
     (void) ratio;
   } 

   /**
    * Function for applying user-defined data processing after
    * refinement.  
    * 
    * Function is overloaded from xfer::RefinePatchStrategy.  An empty
    * implementation is provided here.
    */
   void postprocessRefine(hier::Patch& fine,
                          const hier::Patch& coarse,
                          const hier::Box& fine_box,
                          const hier::IntVector& ratio);

   /**
    * Return maximum stencil width needed for user-defined 
    * data interpolation operations.  Default is to return 
    * zero, assuming no user-defined operations provided.
    */
   hier::IntVector getRefineOpStencilWidth(const tbox::Dimension &dim) const 
     {
       return(hier::IntVector::getZero(dim));
     }

   void getFromInput(const boost::shared_ptr<tbox::Database> &db);

   /**
   * Specify level of diagnostic information printed during iterations.
   * \param print_level
   *        zero prints none or minimal information, higher numbers provide increasingly
   *        verbose debugging information.
   */
   virtual void setDebugPrintInfoLevel(int print_level){d_debug_print_info_level=print_level;}

private:
   /**
    * Blank constructor for PCDiagonalRefinePatchStrategy.
    */
   PCDiagonalRefinePatchStrategy();

   bool contains( const SAMRAI::hier::Box & bigBox,
      const SAMRAI::hier::Box & smallBox,
      const SAMRAI::hier::IntVector & ghosts);

   int d_data_id;
   int d_data_idx;
   std::vector<int> d_bdry_types;
   int d_extrapolation_order;
   int d_debug_print_info_level;

   boost::shared_ptr<hier::PatchHierarchy> d_hierarchy;

};
}
}
#endif
