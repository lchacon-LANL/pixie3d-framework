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

#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/geom/CartesianGridGeometry.h"
#include "SAMRAI/pdat/CellData.h"
#include "boost/shared_ptr.hpp"
#include "SAMRAI/geom/GridGeometry.h"

#include "LevelContainer.h"
#include "PatchContainer.h"
#include "utilities/ProfilerApp.h"

extern "C"{
#include <assert.h>
}

#include "fc_interface.h"
#include "fortran.h"
extern "C" {
#ifdef absoft
  extern void FORTRAN_NAME(APPLYBC) (void*, int*, int*);
  extern void FORTRAN_NAME(POSTPROC)(void*, int*, int*);
  extern void FORTRAN_NAME(FIXDIVB)(void*);
#else
  extern void FORTRAN_NAME(applybc) (void*, int*, int*);
  extern void FORTRAN_NAME(postproc)(void*, int*, int*);
  extern void FORTRAN_NAME(fixdivb)(void*);
#endif
}
namespace SAMRAI{


struct NullDeleter{template<typename T> void operator()(T*){}};


/***********************************************************************
* Constructors/Deconstructor for bcgrp_struct.                         *
***********************************************************************/
// Empty constructor
pixie3dRefinePatchStrategy::bcgrp_struct::bcgrp_struct() {
    nbc_seq = 0;
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
// Number of ints needed to store the struct
size_t pixie3dRefinePatchStrategy::bcgrp_struct::size() {
    return 1+3*nbc_seq;
}
// Pack the data to an int buffer
void pixie3dRefinePatchStrategy::bcgrp_struct::pack(int *buffer) {
    buffer[0] = nbc_seq;
    for (int i=0; i<nbc_seq; i++) {
        buffer[1+3*i+0] = bc_seq[i];
        buffer[1+3*i+1] = vector[i];
        buffer[1+3*i+2] = fillBC[i];
    }
}
// Unpack the data from an int buffer
void pixie3dRefinePatchStrategy::bcgrp_struct::unpack(int *buffer) {
    if ( bc_seq!=NULL ) { delete [] bc_seq; bc_seq=NULL; }
    if ( vector!=NULL ) { delete [] vector; vector=NULL; }
    if ( fillBC!=NULL ) { delete [] fillBC; fillBC=NULL; }
    nbc_seq = buffer[0];
    bc_seq = new int[nbc_seq];
    vector = new int[nbc_seq];
    fillBC = new int[nbc_seq];
    for (int i=0; i<nbc_seq; i++) {
        bc_seq[i] = buffer[1+3*i+0];
        vector[i] = buffer[1+3*i+1];
        fillBC[i] = buffer[1+3*i+2];
    }
}


  
/***********************************************************************
*                                                                      *
* Constructor.                                                         *
*                                                                      *
***********************************************************************/
pixie3dRefinePatchStrategy::pixie3dRefinePatchStrategy(const tbox::Dimension &dim_in):
    RefinePatchStrategy(dim_in),
    dim(dim_in)
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
pixie3dRefinePatchStrategy::setPhysicalBoundaryConditions( hier::Patch& patch,
                                                           const double time,
                                                           const hier::IntVector& ghost_width_to_fill)
{
    //double start = MPI_Wtime();
    if ( ghost_width_to_fill==hier::IntVector(dim,0) ) {
        // We don't need to fill the ghost cells, return
        // This can happen with temporary patches
        return;
    }
    if ( !checkPhysicalBoundary(patch) ) {
        // Patch does not touch physical boundary, no need to do anything, return
        return;
    }

    // Apply the BC (and interior fill) on the patch
    boost::shared_ptr<hier::Patch> patch_ptr = boost::shared_ptr<hier::Patch>(&patch,NullDeleter());
    applyBC( patch_ptr );

    //double stop = MPI_Wtime();
    //printf("   call setPhysicalBoundaryConditions: %i, %f\n",patch.getPatchLevelNumber(),stop-start);
    
}


/***********************************************************************
*                                                                      *
* Preprocessor                                                         *
*                                                                      *
***********************************************************************/
void pixie3dRefinePatchStrategy::preprocessRefine( hier::Patch& patch,
                           const hier::Patch& coarse,
                           const hier::Box& fine_box,
                           const hier::IntVector& ratio ) 
{
} 


/***********************************************************************
*                                                                      *
* Postprocessor                                                        *
*                                                                      *
***********************************************************************/
void pixie3dRefinePatchStrategy::postprocessRefine( hier::Patch& patch,
					       const hier::Patch& coarse,
					       const hier::Box& fine_box,
					       const hier::IntVector& ratio )
{
}


/***********************************************************************
*                                                                      *
* Set the id's for the variables (used for copying the data.            *
*                                                                      *
***********************************************************************/
void
pixie3dRefinePatchStrategy::setPixie3dDataIDs( bool copy, 
    int nvar, int nauxs, int nauxv, int *u0, int *u, int *u_tmp, 
    int *auxs, int *auxs_tmp, int *auxv, int *auxv_tmp, int flux, int src )
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
    flux_id = flux;
    src_id = src;
}


/***********************************************************************
*                                                                      *
* Check if patch touches physical boundaries.                          *
*                                                                      *
***********************************************************************/
bool pixie3dRefinePatchStrategy::checkPhysicalBoundary( hier::Patch& patch)
{
    boost::shared_ptr<hier::PatchGeometry> patchGeom = patch.getPatchGeometry();
    return patchGeom->getTouchesRegularBoundary();
    /*const hier::Box box = patch.getBox();
    const hier::IntVector ratio = patch.getPatchGeometry()->getRatio();
    hier::BoxArray physicalDomain = d_grid_geometry->getPhysicalDomain();
    physicalDomain.refine(ratio);
    hier::Box physicalBox = physicalDomain[0];
    for (int i=0; i<dim->getValue(); i++) {
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
    return(false);*/
}


/***********************************************************************
* Copy data from src_id to dst_id on the given patch                   *
***********************************************************************/
void pixie3dRefinePatchStrategy::copy_data_patch( hier::Patch& patch, int src_id, int dst_id ) {
    boost::shared_ptr< pdat::CellData<double> > src_data = 
        boost::dynamic_pointer_cast<pdat::CellData<double> >( patch.getPatchData(src_id) );
    boost::shared_ptr< pdat::CellData<double> > dst_data = 
        boost::dynamic_pointer_cast<pdat::CellData<double> >( patch.getPatchData(dst_id) );
    if ( dst_data==NULL )
        return;
    assert(src_data!=NULL);
    dst_data->copy(*src_data);
}


/***********************************************************************
* Apply the boundary conditions to a patch                             *
* Note: this will also fill the trivial relations for the interiors    *
***********************************************************************/
void pixie3dRefinePatchStrategy::applyBC( boost::shared_ptr<hier::Patch> patch )
{

    // Some basic error checking
    assert(d_nvar>=0);
    assert(d_nauxs>=0);
    assert(d_nauxv>=0);

    // Copy the data from the temporary patch to the pixie patch
    if ( copy_data ) {
        // Copy u
        for (int i=0; i<d_nvar; i++)
            copy_data_patch(*patch,u_tmp_id[i],u_id[i]);
        // Copy auxillary scalars
        for (int i=0; i<d_nauxs; i++)
            copy_data_patch(*patch,auxs_tmp_id[i],auxs_id[i]);
        // Copy auxillary vectors
        for (int i=0; i<d_nauxv; i++)
            copy_data_patch(*patch,auxv_tmp_id[i],auxv_id[i]);
    }

    // Get the pixie patch (create a temporary patch if necessary)
    void *pixie3d_data = NULL;
    PatchContainer *pixie_patch = NULL;
    if ( patch->inHierarchy() ) {
        // Patch is in the hierarchy, all is well
        const int ln = patch->getPatchLevelNumber();
        LevelContainer *level_container = (LevelContainer *) d_level_container_array[ln];
        pixie3d_data = level_container->getPtr(patch);
    } else {
        // Patch is not in the hierarchy, create a temporary pixie patch
        pixie_patch = new PatchContainer(d_hierarchy,patch,d_nvar,u0_id,u_id,d_nauxs,auxs_id,d_nauxv,auxv_id,flux_id,src_id);        
        pixie3d_data = pixie_patch->getPtr();
    }
    assert(pixie3d_data!=NULL);

    // Correct div(B) = 0
    #ifdef absoft
        FORTRAN_NAME(FIXDIVB)(pixie3d_data);
    #else
        FORTRAN_NAME(fixdivb)(pixie3d_data);
    #endif

    // Apply the boundary conditions (this also fills interior auxillary variables)
    int nbc_seq = d_bc_grp.nbc_seq;
    int *bc_seq = new int[3*nbc_seq];
    for (int i=0; i<nbc_seq; i++) {
        bc_seq[i] = d_bc_grp.bc_seq[i];
        bc_seq[i+d_bc_grp.nbc_seq] = d_bc_grp.vector[i];
        bc_seq[i+2*d_bc_grp.nbc_seq] = d_bc_grp.fillBC[i];
    }
    PROFILE_START("Call applybc");
    #ifdef absoft
        FORTRAN_NAME(APPLYBC)(pixie3d_data,&d_bc_grp.nbc_seq,bc_seq);
    #else
        FORTRAN_NAME(applybc)(pixie3d_data,&d_bc_grp.nbc_seq,bc_seq);
    #endif
    PROFILE_STOP("Call applybc");
    PROFILE_START("Call postproc");
    #ifdef absoft
        FORTRAN_NAME(POSTPROC)(pixie3d_data,&d_bc_grp.nbc_seq,bc_seq);
    #else
        FORTRAN_NAME(postproc)(pixie3d_data,&d_bc_grp.nbc_seq,bc_seq);
    #endif
    PROFILE_STOP("Call postproc");
    delete [] bc_seq;

    // Delete the temporary pixie patch (if used)
    if ( pixie_patch != NULL ) {
        delete pixie_patch;
    }

    // Copy the data from the pixie patch to the temporary patch
    if ( copy_data ) {
        // Copy u         
        for (int i=0; i<d_nvar; i++)
            copy_data_patch(*patch,u_id[i],u_tmp_id[i]);
        // Copy auxillary scalars
        for (int i=0; i<d_nauxs; i++)
            copy_data_patch(*patch,auxs_id[i],auxs_tmp_id[i]);
        // Copy auxillary vectors
        for (int i=0; i<d_nauxv; i++)
            copy_data_patch(*patch,auxv_id[i],auxv_tmp_id[i]);
    }
    
}


}

