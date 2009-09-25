#include "varrayContainer.h"
#include <iostream>

extern "C"{
   void f_create_var_array_(void **p_data, int&);
   void f_delete_var_array_(void *p_data);
   void fill_var_array_(void *p_data, int &, double *, int &, int &, int &, int &, int &, int &);
}

varrayContainer::varrayContainer()
{
   data = NULL;
}

varrayContainer::varrayContainer(tbox::Pointer< hier::Patch<NDIM> > patch, int n_var, int *u_id)
{
   int gcw = 1;
   int xs, xe, ys, ye, zs, ze;
   hier::IntVector<NDIM> gcwc;
   // Get the size of the patch         
   const hier::Index<NDIM> ifirst = patch->getBox().lower();
   const hier::Index<NDIM> ilast  = patch->getBox().upper();
   xs = ifirst(0)+1;
   ys = ifirst(1)+1;
   zs = ifirst(2)+1;
   xe = ilast(0)+1;
   ye = ilast(1)+1;
   ze = ilast(2)+1;
   double *u_ptr[n_var];
   tbox::Pointer< pdat::CellData<NDIM,double> > tmp;
   // Get the pointers to u
   for (int i=0; i<n_var; i++) {
      // Get the pointers to u
      tmp = patch->getPatchData(u_id[i]);
      u_ptr[i] = tmp->getPointer();
      // Check ghost cell width
      gcwc = tmp->getGhostCellWidth();
      for (int j=0; j<NDIM; j++)
         if ( gcw!=gcwc(j) ) { TBOX_ERROR("Ghost Cell width must be 1"); }
    }
   // create the patch object
   CreatePatchFortran(xs,ys,zs,xe,ye,ze,gcw,n_var,u_ptr);
}


void varrayContainer::CreatePatchFortran( int xs, int ys, int zs, int xe, int ye, int ze,
   int gcw, int n_var, double **u_ptr)
{
   int xsg, ysg, zsg, xeg, yeg, zeg;
   xsg = xs-gcw;
   ysg = ys-gcw;
   zsg = zs-gcw;
   xeg = xe+gcw;
   yeg = ye+gcw;
   zeg = ze+gcw;
   // Create var_array for u
   f_create_var_array_(&data,n_var);
   for(int i=0; i<n_var;i++) {
      fill_var_array_(data,i,u_ptr[i],xsg,ysg,zsg,xeg,yeg,zeg);
   }
}

varrayContainer::~varrayContainer()
{
   // Delete the var_array on the fortran side
   f_delete_var_array_(data);
   data = NULL;
}

