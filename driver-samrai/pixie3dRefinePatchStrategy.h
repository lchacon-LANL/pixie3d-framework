//
// $Id: pixie3dRefinePatchStrategy.h 1828 2005-08-23 19:51:19Z pernice $
// $Revision: 1828 $
// $Date: 2005-08-23 13:51:19 -0600 (Tue, 23 Aug 2005) $
//
// File:  pixie3dRefinePatchStrategy.h
// Copyright: (c) 2002 The Regents of the University of California
// Description:  Header file for applying boundary conditions to magnetic field variable
//

#ifndef included_pixie3dRefinePatchStrategy
#define included_pixie3dRefinePatchStrategy
#endif

#include "SAMRAI_config.h"

#include "Box.h"
#include "IntVector.h"
#include "Patch.h"
#include "RefinePatchStrategy.h"
#include "GridGeometry.h"

/*
*************************************************************************
*                                                                       *
* External declarations for FORTRAN 77 routines                         *
*                                                                       *
*************************************************************************
*/

namespace SAMRAI{

class pixie3dRefinePatchStrategy : public xfer::RefinePatchStrategy<NDIM>
{
public:
   /**
    * Constructor.
    */
   pixie3dRefinePatchStrategy(void);

   /**
    * Virtual destructor.
    */
   virtual ~pixie3dRefinePatchStrategy();

   /**
    * Set solution ghost cell values along physical boundaries.
    *
    * Function is overloaded from xfer_RefinePatchStrategyX.
    */
   void setPhysicalBoundaryConditions( hier::Patch<NDIM>& patch,
                                       const double time,
                                       const hier::IntVector<NDIM>& ghost_width_to_fill );

   /**
    * Function for applying user-defined data processing prior
    * to refinement.  
    *
    * Function is overloaded from xfer\_RefinePatchStrategyX.  An empty
    * implementation is provided here.
    */
   void preprocessRefine( hier::Patch<NDIM>& fine,
                          const hier::Patch<NDIM>& coarse,
                          const hier::Box<NDIM>& fine_box,
                          const hier::IntVector<NDIM>& ratio )
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
    * Function is overloaded from xfer\_RefinePatchStrategyX.  An empty
    * implementation is provided here.
    */
   void postprocessRefine( hier::Patch<NDIM>& fine,
                           const hier::Patch<NDIM>& coarse,
                           const hier::Box<NDIM>& fine_box,
                           const hier::IntVector<NDIM>& ratio ) 
   {
      (void) fine;
      (void) coarse;
      (void) fine_box;
      (void) ratio;
   } 

   /**
    * Return maximum stencil width needed for user-defined 
    * data interpolation operations.  Default is to return 
    * zero, assuming no user-defined operations provided.
    */
   hier::IntVector<NDIM> getRefineOpStencilWidth( void ) const 
     {
       return(hier::IntVector<NDIM>(0));
     }

   /**
    * Set descriptor index of data locations where boundary conditions
    * must be set.
    */
   void setRefineStrategyDataId( const int var_id ){ d_data_id = var_id;}

   
   void setGridGeometry(tbox::Pointer< hier::GridGeometry<NDIM > > grid_geometry){d_grid_geometry=grid_geometry;}

   void setPixie3dHierarchyData(void **hierarchy_data){d_level_container_array = hierarchy_data; }

   void setPixie3dDataIDs(bool copy, int *u, int *u_tmp, int *auxs, int *auxs_tmp, int *auxv, int *auxv_tmp );

   bool checkPhysicalBoundary( hier::Patch<NDIM>& patch);

private:

   int d_data_id;
   bool copy_data;
   int u_id[NVAR], auxs_id[NAUXS], auxv_id[NAUXV];
   int u_tmp_id[NVAR], auxs_tmp_id[NAUXS], auxv_tmp_id[NAUXV];


   tbox::Pointer< hier::GridGeometry<NDIM > > d_grid_geometry;

   void **d_level_container_array;

};
}