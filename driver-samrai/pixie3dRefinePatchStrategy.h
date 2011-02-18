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

#include "SAMRAI_config.h"

#include "Box.h"
#include "IntVector.h"
#include "Patch.h"
#include "RefinePatchStrategy.h"
#include "GridGeometry.h"
#include "PatchHierarchy.h"



namespace SAMRAI{

class pixie3dRefinePatchStrategy : public xfer::RefinePatchStrategy<NDIM>
{
public:


    // Constructor.
    pixie3dRefinePatchStrategy(void);

    // Virtual destructor.
    virtual ~pixie3dRefinePatchStrategy();


    // Set solution ghost cell values along physical boundaries.
    // Function is overloaded from xfer_RefinePatchStrategyX.
    void setPhysicalBoundaryConditions( hier::Patch<NDIM>& patch,
                                        const double time,
                                        const hier::IntVector<NDIM>& ghost_width_to_fill );


    /* Function for applying user-defined data processing prior to refinement.  
     * This function is called after setPhysicalBoundaryConditions.
     * Function is overloaded from xfer_RefinePatchStrategyX.  
     */
    void preprocessRefine( hier::Patch<NDIM>& fine,
                           const hier::Patch<NDIM>& coarse,
                           const hier::Box<NDIM>& fine_box,
                           const hier::IntVector<NDIM>& ratio );


    // Function for applying user-defined data processing after refinement.  
    // Function is overloaded from xfer\_RefinePatchStrategyX. 
    void postprocessRefine( hier::Patch<NDIM>& fine,
                           const hier::Patch<NDIM>& coarse,
                           const hier::Box<NDIM>& fine_box,
                           const hier::IntVector<NDIM>& ratio );


    /* Return maximum stencil width needed for user-defined 
     * data interpolation operations.  Default is to return 
     * zero, assuming no user-defined operations provided.
     */
    hier::IntVector<NDIM> getRefineOpStencilWidth( void ) const  { return(hier::IntVector<NDIM>(0)); }


    // Set the hierarchy
    void setHierarchy(tbox::Pointer< hier::PatchHierarchy<NDIM> > hierarchy) { d_hierarchy=hierarchy; }

    // Set the grid geometry
    void setGridGeometry(tbox::Pointer< hier::GridGeometry<NDIM > > grid_geometry) { d_grid_geometry=grid_geometry; }

    // Set the level contanier data
    void setPixie3dHierarchyData(void **hierarchy_data) { d_level_container_array = hierarchy_data; }

    // Set the ids of the data
    void setPixie3dDataIDs(bool copy, int nvar, int nauxs, int nauxv, int *u0, int *u, int *u_tmp, int *auxs, int *auxs_tmp, int *auxv, int *auxv_tmp );

    // Check the physical boundaries
    bool checkPhysicalBoundary( hier::Patch<NDIM>& patch);

    // Structures used to hold a single boundary condition group
    struct bcgrp_struct {
        int nbc_seq;        // Number of members of the current sequence
        int *bc_seq;        // The sequence ids (>0: primary variables, <0: auxillary variable)
        int *vector;        // Is the variable a vector (1) or scalar (0)
        int *fillBC;        // Do we need to fill BC (1), or are only the interiors needed (0)
        bcgrp_struct();
        bcgrp_struct(const int);
        bcgrp_struct(const bcgrp_struct&);
        bcgrp_struct& operator=(const bcgrp_struct&);
        ~bcgrp_struct();
    };

    // Function to set the boundary condition group that we are processing
    void setRefineStrategySequence( const bcgrp_struct bc_grp ) { d_bc_grp = bc_grp; }


private:
   
    // Private function to copy data on a patch from one variable to another
    void copy_data_patch( hier::Patch<NDIM>& patch, int src_id, int dst_id );
    
    // Information about the hierarchy
    tbox::Pointer< hier::PatchHierarchy<NDIM> > d_hierarchy;
    tbox::Pointer< hier::GridGeometry<NDIM > > d_grid_geometry;

    // The number of variables and their ids
    int d_nvar;
    int d_nauxs;
    int d_nauxv;
    bool copy_data;
    int *u0_id, *u_id, *auxs_id, *auxv_id;
    int *u_tmp_id, *auxs_tmp_id, *auxv_tmp_id;
   
    // The level contanier array
    void **d_level_container_array;

    // The boundry condition group
    bcgrp_struct d_bc_grp;

};

}
#endif

