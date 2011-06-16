//
// $Id: pixie3dRefinePatchStrategy.C 1828 2005-08-23 19:51:19Z pernice $
// $Revision: 1828 $
// $Date: 2005-08-23 13:51:19 -0600 (Tue, 23 Aug 2005) $
//
// File:  pixie3dRefinePatchStrategy.C
// Copyright: (c) 2002 The Regents of the University of California
// Description:  Implementation for applying boundary conditions to 
//               magnetic field.
//

#include "pixie3dRefinePatchStrategy.h"

#include "CartesianPatchGeometry.h"
#include "CellData.h"
#include "tbox/Pointer.h"
#include "LevelContainer.h"
#include "PatchContainer.h"
#include "GridGeometry.h"
#include "CartesianGridGeometry.h"

extern "C"{
#include <assert.h>
}

#include "fc_interface.h"
#include "fortran.h"
extern "C" {
#ifdef absoft
  extern void FORTRAN_NAME(APPLYBC) (void*, int*, int*);
#else
  extern void FORTRAN_NAME(applybc) (void*, int*, int*);
#endif
}
namespace SAMRAI{


/***********************************************************************
* Constructors/Deconstructor for bcgrp_struct.                         *
***********************************************************************/
// Empty constructor
pixie3dRefinePatchStrategy::bcgrp_struct::bcgrp_struct() {
    nbc_seq = -1;
    bc_seq = NULL;
    vector = NULL;
    fillBC = NULL;
}
// Constructor used to initialize the structure
pixie3dRefinePatchStrategy::bcgrp_struct::bcgrp_struct(const int N) {
    nbc_seq = N;
    bc_seq = new int[nbc_seq];
    vector = new int[nbc_seq];
    fillBC = new int[nbc_seq];
    for (int i=0; i<nbc_seq; i++) {
        bc_seq[i] = -1;
        vector[i] = -1;
        fillBC[i] = -1;
    }
}
// Copy constructor
pixie3dRefinePatchStrategy::bcgrp_struct::bcgrp_struct(const bcgrp_struct & rhs) {
    nbc_seq = rhs.nbc_seq;
    bc_seq = new int[nbc_seq];
    vector = new int[nbc_seq];
    fillBC = new int[nbc_seq];
    for (int i=0; i<nbc_seq; i++) {
        bc_seq[i] = rhs.bc_seq[i];
        vector[i] = rhs.vector[i];
        fillBC[i] = rhs.fillBC[i];
    }
}
// Assignment operator
pixie3dRefinePatchStrategy::bcgrp_struct& pixie3dRefinePatchStrategy::bcgrp_struct::operator=( const bcgrp_struct& rhs ) {
    if ( bc_seq!=NULL ) { delete [] bc_seq; bc_seq=NULL; }
    if ( vector!=NULL ) { delete [] vector; vector=NULL; }
    if ( fillBC!=NULL ) { delete [] fillBC; fillBC=NULL; }
    nbc_seq = rhs.nbc_seq;
    bc_seq = new int[nbc_seq];
    vector = new int[nbc_seq];
    fillBC = new int[nbc_seq];
    for (int i=0; i<nbc_seq; i++) {
        bc_seq[i] = rhs.bc_seq[i];
        vector[i] = rhs.vector[i];
        fillBC[i] = rhs.fillBC[i];
    }
    return *this;
}
// De-constructor
pixie3dRefinePatchStrategy::bcgrp_struct::~bcgrp_struct() {
    nbc_seq = 0;
    if ( bc_seq!=NULL ) { delete [] bc_seq; bc_seq=NULL; }
    if ( vector!=NULL ) { delete [] vector; vector=NULL; }
    if ( fillBC!=NULL ) { delete [] fillBC; fillBC=NULL; }
    bc_seq = NULL;
    vector = NULL;
    fillBC = NULL;
}

  
/***********************************************************************
*                                                                      *
* Constructor.                                                         *
*                                                                      *
***********************************************************************/
pixie3dRefinePatchStrategy::pixie3dRefinePatchStrategy(void)
{
    u0_id = NULL;
    u_id = NULL;
    u_tmp_id = NULL;
    auxs_id = NULL;
    auxs_tmp_id = NULL;
    auxv_id = NULL;
    auxv_tmp_id = NULL;
}


/***********************************************************************
*                                                                      *
* Destructor.                                                          *
*                                                                      *
***********************************************************************/
pixie3dRefinePatchStrategy::~pixie3dRefinePatchStrategy()
{
    delete [] u0_id;
    delete [] u_id;
    delete [] u_tmp_id;
    delete [] auxs_id;
    delete [] auxs_tmp_id;
    delete [] auxv_id;
    delete [] auxv_tmp_id;
}


/***********************************************************************
*                                                                      *
* Set ghost cell values at physical boundaries.                        *
*                                                                      *
***********************************************************************/
void 
pixie3dRefinePatchStrategy::setPhysicalBoundaryConditions( hier::Patch<NDIM>& patch,
                                                           const double time,
                                                           const hier::IntVector<NDIM>& ghost_width_to_fill)
{
    assert(d_nvar>=0);
    assert(d_nauxs>=0);
    assert(d_nauxv>=0);
    if ( ghost_width_to_fill==hier::IntVector<NDIM>(0) ) {
        // We don't need to fill the ghost cells, return
        // This can happen with temporary patches
        return;
    }
    /*if ( !checkPhysicalBoundary(patch) ) {
        // Patch does not touch physical boundary, no need to do anything, return
        return;
    }*/

    // Get some basic patch info
    const int ln = patch.getPatchLevelNumber();
    const int pn = patch.getPatchNumber();

    // Copy the data from the temporary patch to the pixie patch
    if ( copy_data ) {
        // Copy u
        for (int i=0; i<d_nvar; i++)
            copy_data_patch(patch,u_tmp_id[i],u_id[i]);
        // Copy auxillary scalars
        for (int i=0; i<d_nauxs; i++)
            copy_data_patch(patch,auxs_tmp_id[i],auxs_id[i]);
        // Copy auxillary vectors
        for (int i=0; i<d_nauxv; i++)
            copy_data_patch(patch,auxv_tmp_id[i],auxv_id[i]);
    }

    // Get the pixie patch (create a temporary patch if necessary)
    LevelContainer *level_container = NULL;
    if (patch.inHierarchy()) {
        // Patch is in the hierarchy, all is well
        level_container = (LevelContainer *) d_level_container_array[ln];
    } else {
//TBOX_ERROR("Not tested");
        // Patch is not in the hierarchy, create a temporary pixie patch
        level_container = new LevelContainer::LevelContainer(pn+1,d_hierarchy,
            d_nvar,u0_id,u_id,d_nauxs,auxs_id,d_nauxv,auxv_id);
        tbox::Pointer< hier::Patch<NDIM> > patch_ptr = tbox::Pointer< hier::Patch<NDIM> >::Pointer(&patch,false);
        level_container->CreatePatch(pn,patch_ptr);
    }
    assert(level_container!=NULL);
    void *pixie3d_data = level_container->getPtr(pn);
    assert(pixie3d_data!=NULL);

    // Apply the boundary conditions
    int nbc_seq = d_bc_grp.nbc_seq;
    int *bc_seq = new int[3*nbc_seq];
    for (int i=0; i<nbc_seq; i++) {
        bc_seq[i] = d_bc_grp.bc_seq[i];
        bc_seq[i+d_bc_grp.nbc_seq] = d_bc_grp.vector[i];
        bc_seq[i+2*d_bc_grp.nbc_seq] = d_bc_grp.fillBC[i];
    }
    #ifdef absoft
        FORTRAN_NAME(APPLYBC)(pixie3d_data,&d_bc_grp.nbc_seq,bc_seq);
    #else
        FORTRAN_NAME(applybc)(pixie3d_data,&d_bc_grp.nbc_seq,bc_seq);
    #endif
    delete [] bc_seq;

    // Delete the temporary pixie patch (if used)
    if (!patch.inHierarchy()) {
        delete level_container;
    }

    // Copy the data from the pixie patch to the temporary patch
    if ( copy_data ) {
        // Copy u         
        for (int i=0; i<d_nvar; i++)
            copy_data_patch(patch,u_id[i],u_tmp_id[i]);
        // Copy auxillary scalars
        for (int i=0; i<d_nauxs; i++)
            copy_data_patch(patch,auxs_id[i],auxs_tmp_id[i]);
        // Copy auxillary vectors
        for (int i=0; i<d_nauxv; i++)
            copy_data_patch(patch,auxv_id[i],auxv_tmp_id[i]);
    }
    
}


/***********************************************************************
*                                                                      *
* Preprocessor                                                         *
*                                                                      *
***********************************************************************/
void pixie3dRefinePatchStrategy::preprocessRefine( hier::Patch<NDIM>& fine,
                           const hier::Patch<NDIM>& coarse,
                           const hier::Box<NDIM>& fine_box,
                           const hier::IntVector<NDIM>& ratio ) 
{
    (void) fine;
    (void) coarse;
    (void) fine_box;
    (void) ratio;
} 


/***********************************************************************
*                                                                      *
* Postprocessor                                                        *
*                                                                      *
***********************************************************************/
void pixie3dRefinePatchStrategy::postprocessRefine( hier::Patch<NDIM>& patch,
					       const hier::Patch<NDIM>& coarse,
					       const hier::Box<NDIM>& fine_box,
					       const hier::IntVector<NDIM>& ratio )
{
    assert(d_nvar>=0);
    assert(d_nauxs>=0);
    assert(d_nauxv>=0);
    bool touches_boundary = this->checkPhysicalBoundary(patch);
    /*if ( touches_boundary ) {
        // Return for now if we called applybc, eventually we will call the post processor on all patches
        return;
    }*/


    // Perform a post processor step in which we fill all the trivial relationships
    // This needs to be done on all patches   
    if (patch.inHierarchy()) {

        // Patch is in the hierarchy, all is well
        const int ln = patch.getPatchLevelNumber();
        const int pn = patch.getPatchNumber();
	  
        LevelContainer *level_container = (LevelContainer *) d_level_container_array[ln];
        assert(level_container!=NULL);
        void *pixie3d_data = level_container->getPtr(pn);
        assert(pixie3d_data!=NULL);
 
        // Copy the data from the temporary patch to the pixie patch
        if ( copy_data ) {
            // Copy u         
            for (int i=0; i<d_nvar; i++)
                copy_data_patch(patch,u_tmp_id[i],u_id[i]);
            // Copy auxillary scalars
            for (int i=0; i<d_nauxs; i++)
                copy_data_patch(patch,auxs_tmp_id[i],auxs_id[i]);
            // Copy auxillary vectors
            for (int i=0; i<d_nauxv; i++)
                copy_data_patch(patch,auxv_tmp_id[i],auxv_id[i]);
        }
      
        // Apply the post-processor to fill the algebraic auxillary variables
        int *bc_seq = new int[3*d_bc_grp.nbc_seq];
        for (int i=0; i<d_bc_grp.nbc_seq; i++) {
            bc_seq[i] = d_bc_grp.bc_seq[i];
            bc_seq[i+d_bc_grp.nbc_seq] = d_bc_grp.vector[i];
            bc_seq[i+2*d_bc_grp.nbc_seq] = d_bc_grp.fillBC[i];
        }
        #ifdef absoft
            FORTRAN_NAME(APPLYBC)(pixie3d_data, &d_bc_grp.nbc_seq, bc_seq);
        #else
            FORTRAN_NAME(applybc)(pixie3d_data, &d_bc_grp.nbc_seq, bc_seq);
        #endif
        delete [] bc_seq;
	  
	    // Copy the data from the pixie patch to the temporary patch
        if ( copy_data ) {
            // Copy u         
            for (int i=0; i<d_nvar; i++)
                copy_data_patch(patch,u_id[i],u_tmp_id[i]);
            // Copy auxillary scalars
            for (int i=0; i<d_nauxs; i++)
                copy_data_patch(patch,auxs_id[i],auxs_tmp_id[i]);
            // Copy auxillary vectors
            for (int i=0; i<d_nauxv; i++)
                copy_data_patch(patch,auxv_id[i],auxv_tmp_id[i]);
        }
    } else {
        // Patch is not in the hierarchy, we should not need to apply the post processor
        abort();
    }
}


/***********************************************************************
*                                                                      *
* Set the id's for the variables (used for copying the data.            *
*                                                                      *
***********************************************************************/
void
pixie3dRefinePatchStrategy::setPixie3dDataIDs(bool copy, 
					      int nvar, int nauxs, int nauxv,
					      int *u0, int *u, int *u_tmp, int *auxs, int *auxs_tmp, int *auxv, int *auxv_tmp )
{
    // Copy some basic information
    d_nvar = nvar;
    d_nauxs = nauxs;
    d_nauxv = nauxv;
    copy_data = copy;

    // Deallocate existing arrays if they exist
    if ( u0_id != NULL )
        delete [] u0_id;
    if ( u_id != NULL )
        delete [] u_id;
    if ( u_tmp_id!=NULL )
        delete [] u_tmp_id;
    if ( auxs_id!=NULL )
        delete [] auxs_id;
    if ( auxs_tmp_id!=NULL )
        delete [] auxs_tmp_id;
    if ( auxv_id!=NULL )
        delete [] auxv_id;
    if ( auxv_tmp_id!=NULL )
        delete [] auxv_tmp_id;

    // Allocate and fill id arrays
    u0_id = new int[d_nvar];
    u_id = new int[d_nvar];
    u_tmp_id = new int[d_nvar];
    auxs_id = new int[d_nauxs];
    auxs_tmp_id = new int[d_nauxs];
    auxv_id = new int[d_nauxv];
    auxv_tmp_id = new int[d_nauxv];
    for (int i=0; i<d_nvar; i++) {
        u0_id[i] = u0[i];
        u_id[i] = u[i];
        u_tmp_id[i] = u_tmp[i];  
	}
    for (int i=0; i<d_nauxs; i++) {
        auxs_id[i] = auxs[i];
        auxs_tmp_id[i] = auxs_tmp[i];
	}
    for (int i=0; i<d_nauxv; i++) {
        auxv_id[i] = auxv[i];
        auxv_tmp_id[i] = auxv_tmp[i];
	}
}


/***********************************************************************
*                                                                      *
* Check if patch touches physical boundaries.                          *
* For temporary patches looks like SAMRAI has a bug
*                                                                      *
***********************************************************************/
bool pixie3dRefinePatchStrategy::checkPhysicalBoundary( hier::Patch<NDIM>& patch)
{
    tbox::Pointer< hier::PatchGeometry<NDIM> > patchGeom = patch.getPatchGeometry();
    const hier::Box<NDIM> box = patch.getBox();
    const hier::IntVector< NDIM > ratio = patch.getPatchGeometry()->getRatio();
    hier::BoxArray<NDIM> physicalDomain = d_grid_geometry->getPhysicalDomain();
    physicalDomain.refine(ratio);
    hier::Box<NDIM> physicalBox = physicalDomain[0];
    for (int i=0; i<NDIM; i++) {
        if (patchGeom->getTouchesRegularBoundary(i,0)) {
            const int patch_lower = box.lower(i);
            if ( patch_lower <= 0 )
                return true;
        }
        if (patchGeom->getTouchesRegularBoundary(i,1)) {
            const int patch_upper = box.upper(i);
            if ( patch_upper >= physicalBox.numberCells(i)-1 )
                return true;
        }
    }
    return(false);
}


/***********************************************************************
* Copy data from src_id to dst_id on the given patch                   *
***********************************************************************/
void pixie3dRefinePatchStrategy::copy_data_patch( hier::Patch<NDIM>& patch, int src_id, int dst_id ) {
    tbox::Pointer< pdat::CellData<NDIM,double> > src_data = patch.getPatchData(src_id);
    tbox::Pointer< pdat::CellData<NDIM,double> > dst_data = patch.getPatchData(dst_id);
    if ( dst_data.isNull() )
        return;
    assert(!src_data.isNull());
    dst_data->copy(*src_data);
}



}

