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
#include "GridGeometry.h"

extern "C"{
#include <assert.h>
}

#include "fc_interface.h"
#include "fortran.h"
extern "C" {
#ifdef absoft
extern void FORTRAN_NAME(APPLYBC) (int*, void*);
#else
extern void FORTRAN_NAME(applybc) (int*, void*);
#endif
}
namespace SAMRAI{
/*
************************************************************************
*                                                                      *
* Constructor.                                                         *
*                                                                      *
************************************************************************
*/
pixie3dRefinePatchStrategy::pixie3dRefinePatchStrategy(void)
{
  d_data_id = -1;
  
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
         for (int i=0; i<d_nvar; i++) {
            tbox::Pointer< pdat::CellData<NDIM,double> > tmp1 = patch.getPatchData(u_tmp_id[i]);
            tbox::Pointer< pdat::CellData<NDIM,double> > tmp2 = patch.getPatchData(u_id[i]);
            tmp2->copy(*tmp1);
         }
         // Copy auxillary scalars
         for (int i=0; i<d_nauxs; i++) {
            tbox::Pointer< pdat::CellData<NDIM,double> > tmp1 = patch.getPatchData(auxs_tmp_id[i]);
            tbox::Pointer< pdat::CellData<NDIM,double> > tmp2 = patch.getPatchData(auxs_id[i]);
            tmp2->copy(*tmp1);
         }
         // Copy auxillary vectors
         for (int i=0; i<d_nauxv; i++) {
            tbox::Pointer< pdat::CellData<NDIM,double> > tmp1 = patch.getPatchData(auxv_tmp_id[i]);
            tbox::Pointer< pdat::CellData<NDIM,double> > tmp2 = patch.getPatchData(auxv_id[i]);
            tmp2->copy(*tmp1);
         }
      }

      // f_apply_pixie3d_bc_(&d_data_id, pixie3d_data);
      #ifdef absoft
         FORTRAN_NAME(APPLYBC)(&d_data_id, pixie3d_data);
      #else
         FORTRAN_NAME(applybc)(&d_data_id, pixie3d_data);
      #endif

      // Copy the data from the pixie patch to the temporary patch
      if ( copy_data ) {
         // Copy u         
         for (int i=0; i<d_nvar; i++) {
            tbox::Pointer< pdat::CellData<NDIM,double> > tmp1 = patch.getPatchData(u_id[i]);
            tbox::Pointer< pdat::CellData<NDIM,double> > tmp2 = patch.getPatchData(u_tmp_id[i]);
            tmp2->copy(*tmp1);
         }
         // Copy auxillary scalars
         for (int i=0; i<d_nauxs; i++) {
            tbox::Pointer< pdat::CellData<NDIM,double> > tmp1 = patch.getPatchData(auxs_id[i]);
            tbox::Pointer< pdat::CellData<NDIM,double> > tmp2 = patch.getPatchData(auxs_tmp_id[i]);
            tmp2->copy(*tmp1);
         }
         // Copy auxillary vectors
         for (int i=0; i<d_nauxv; i++) {
            tbox::Pointer< pdat::CellData<NDIM,double> > tmp1 = patch.getPatchData(auxv_id[i]);
            tbox::Pointer< pdat::CellData<NDIM,double> > tmp2 = patch.getPatchData(auxv_tmp_id[i]);
            tmp2->copy(*tmp1);
         }
      }
   } else if( checkPhysicalBoundary(patch) ) {
      // Patch is not in the hierarchy and touches the boundary (this is not finished)
      const hier::IntVector< NDIM > ratio = patch.getPatchGeometry()->getRatio();
      hier::BoxArray<NDIM> physicalDomain = d_grid_geometry->getPhysicalDomain();
      physicalDomain.refine(ratio);

      // we assume the domain is a single box
      if(d_grid_geometry->getDomainIsSingleBox()) {
         hier::Box<NDIM> physicalBox = physicalDomain(0);
         const int nx = physicalBox.numberCells(0);
         const int ny = physicalBox.numberCells(1);
         const int nz = physicalBox.numberCells(2);
      } else {
         abort();
      }
      abort();
   } else {
      // Patch is not in the hierarchy and does not touch the boundary
      // No need to do anything
   }

}


/***********************************************************************
*                                                                      *
* Set the id's for the variables (used for coping the data.            *
*                                                                      *
***********************************************************************/
void
pixie3dRefinePatchStrategy::setPixie3dDataIDs(bool copy, int *u,
					      int nvar, int nauxs, int nauxv,
					      int *u_tmp, int *auxs, int *auxs_tmp, int *auxv, int *auxv_tmp )
{
  d_nvar = nvar;
  d_nauxs = nauxs;
  d_nauxv = nauxv;

  if(u_id==NULL)
    {
      u_id = new int[d_nvar];
    }
  
  if(u_tmp_id==NULL)
    {
      u_tmp_id = new int[d_nvar];
    }

  if(auxs_id==NULL)
    {
      auxs_id = new int[d_nauxs];
    }
  
  if(auxs_tmp_id==NULL)
    {
      auxs_tmp_id = new int[d_nauxs];
    }
  
  if(auxv_id==NULL)
    {
      auxv_id = new int[d_nauxs];
    }
  
  if(auxv_tmp_id==NULL)
    {
      auxv_tmp_id = new int[d_nauxs];
    }
  
  copy_data = copy;
  if ( copy )
    {
      for (int i=0; i<d_nvar; i++)
	{
	  u_id[i] = u[i];
	  u_tmp_id[i] = u_tmp[i];  
	}
      
      for (int i=0; i<d_nauxs; i++)
	{
	  auxs_id[i] = auxs[i];
	  auxs_tmp_id[i] = auxs_tmp[i];
	}
      
      for (int i=0; i<d_nauxv; i++)
	{
	  auxv_id[i] = auxv[i];
	  auxv_tmp_id[i] = auxv_tmp[i];
	}
    }
}


/***********************************************************************
*                                                                      *
* Check if patch touches physical boundaries.                          *
*                                                                      *
***********************************************************************/
bool 
pixie3dRefinePatchStrategy::checkPhysicalBoundary( hier::Patch<NDIM>& patch)
{
   tbox::Pointer< hier::PatchGeometry<NDIM> > patchGeom = patch.getPatchGeometry();
   const hier::Box<NDIM> box = patch.getBox();
   const hier::IntVector< NDIM > ratio = patch.getPatchGeometry()->getRatio();
   hier::BoxArray<NDIM> physicalDomain = d_grid_geometry->getPhysicalDomain();
   physicalDomain.refine(ratio);
   hier::Box<NDIM> physicalBox = physicalDomain(0);
   for (int i=0; i<NDIM; i++) {
      if (patchGeom->getTouchesRegularBoundary(i,0)) {
         const int patch_lower = box.lower(i);
         if ( patch_lower <= 0 )
            return true;
      }
      if (patchGeom->getTouchesRegularBoundary(i,1)) {
         const int patch_upper = box.upper(i);
         if ( patch_upper >= physicalBox.numberCells(i) )
            return true;
      }
   }
   return(false);
}


}
