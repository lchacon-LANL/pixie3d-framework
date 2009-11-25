#include "LevelContainer.h"
#include <iostream>

extern "C"{
   void f_create_var_array_(void **p_data, int&);
   void f_create_aux_array_(void **p_data, int&, int&);
   void f_create_grid_mg_def_(void **p_data, int &, int &, int &, int &, int &, int &, int &, int &, int &, int &);
   void f_create_patch_data_(void **p_data, void *, void *, void *, void *);
   void f_delete_patch_data_(void *p_data);
   void fill_var_array_(void *p_data, int &, double *, int &, int &, int &, int &, int &, int &);
   void fill_aux_array_var_(void *p_data, int &, double *, int &, int &, int &, int &, int &, int &);
   void fill_aux_array_vec_(void *p_data, int &, double *, int &, int &, int &, int &, int &, int &);
}

LevelContainer::LevelContainer()
{
   data = NULL;
}

LevelContainer::LevelContainer(int n, int nx, int ny, int nz)
{
   N = n;
   // Create space for a patch object for each patch (n patches)
   data = new void *[N];
   u    = new void *[N];
   u0   = new void *[N];
   aux  = new void *[N];
   gparams = new void *[N];
   filled = new int [N];

   for (int i=0; i<N; i++)
     {
       filled[i] = 0;
     }

   // Set global number of points for the level
   nglx = nx;
   ngly = ny;
   nglz = nz;
}


void LevelContainer::CreatePatch(int patch_id, tbox::Pointer< hier::Patch<NDIM> > patch,
      int n_var, int *u0_id, int *u_id, int n_auxs, int *auxs_id, int n_auxv, int *auxv_id)
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
   double *u_ptr[n_var], *u0_ptr [n_var], *auxs_ptr[n_auxs], *auxv_ptr[n_auxv];
   tbox::Pointer< pdat::CellData<NDIM,double> > tmp;

   // Get the pointers to u_0, u_n
   for (int i=0; i<n_var; i++) 
     {
       // Get the pointers to u_0
       tmp = patch->getPatchData(u0_id[i]);
       u0_ptr[i] = tmp->getPointer();
       // Check ghost cell width
       gcwc = tmp->getGhostCellWidth();
       for (int j=0; j<NDIM; j++)
	 {
	   if ( gcw!=gcwc(j) ) 
	     { 
	       TBOX_ERROR("Ghost Cell width must be 1"); 
	     }
	 }

       // Get the pointers to u_n
       tmp = patch->getPatchData(u_id[i]);
       u_ptr[i] = tmp->getPointer();
       // Check ghost cell width
       gcwc = tmp->getGhostCellWidth();
      for (int j=0; j<NDIM; j++)
	{
	  if ( gcw!=gcwc(j) ) 
	    { 
	      TBOX_ERROR("Ghost Cell width must be 1"); 
	    }
	}
     }

   // Get the pointers to the auxillary variables
   for (int i=0; i<n_auxs; i++) 
     {
       // Get the pointers to scalar variables
       tmp = patch->getPatchData(auxs_id[i]);
       auxs_ptr[i] = tmp->getPointer();
       // Check ghost cell width
       gcwc = tmp->getGhostCellWidth();

       for (int j=0; j<NDIM; j++)
	 {
	   if ( gcw!=gcwc(j) ) 
	     { 
	       TBOX_ERROR("Ghost Cell width must be 1"); 
	     }
	 }
     }

   for (int i=0; i<n_auxv; i++) 
     {
       // Get the pointers to scalar variables
       tmp = patch->getPatchData(auxv_id[i]);
       auxv_ptr[i] = tmp->getPointer();
       // Check ghost cell width
       gcwc = tmp->getGhostCellWidth();
       for (int j=0; j<NDIM; j++)
	 {
	   if ( gcw!=gcwc(j) ) 
	     {
	       TBOX_ERROR("Ghost Cell width must be 1"); 
	     }
	 }
     }

   // create the patch object
   CreatePatchFortran(patch_id,xs,ys,zs,xe,ye,ze,gcw,n_var,u0_ptr,u_ptr,n_auxs,auxs_ptr,n_auxv,auxv_ptr);
}


void LevelContainer::CreatePatchFortran( int patch, int xs, int ys, int zs, int xe, int ye, int ze, int gcw, 
   int n_var, double **u0_ptr, double **u_ptr, int n_auxs, double **auxs_ptr, int n_auxv, double **auxv_ptr)
{
   int xsg, ysg, zsg, xeg, yeg, zeg;
   xsg = xs-gcw;
   ysg = ys-gcw;
   zsg = zs-gcw;
   xeg = xe+gcw;
   yeg = ye+gcw;
   zeg = ze+gcw;
   int xsgt, ysgt, zsgt, xegt, yegt, zegt;
   xsgt = xsg-xs+1;
   xegt = xeg-xs+1;
   ysgt = ysg-ys+1;
   yegt = yeg-ys+1;
   zsgt = zsg-zs+1;
   zegt = zeg-zs+1;
   
   // Create var_array for u, u0
   f_create_var_array_(&u[patch],n_var);
   f_create_var_array_(&u0[patch],n_var);

   for(int i=0; i<n_var;i++) 
     {
       fill_var_array_( u[patch],i, u_ptr[i],xsg,ysg,zsg,xeg,yeg,zeg);
       fill_var_array_(u0[patch],i,u0_ptr[i],xsg,ysg,zsg,xeg,yeg,zeg);
     }

   // Create aux_array auxillary variables
   f_create_aux_array_(&aux[patch],n_auxs,n_auxv);

   for(int i=0; i<n_auxs;i++)
     {
       fill_aux_array_var_(aux[patch],i,auxs_ptr[i],xsgt,ysgt,zsgt,xegt,yegt,zegt);
     }

   for(int i=0; i<n_auxv;i++)
     {
       fill_aux_array_vec_(aux[patch],i,auxv_ptr[i],xsgt,ysgt,zsgt,xegt,yegt,zegt);
     }

   // Create grid_mg_def for the patch information
   f_create_grid_mg_def_(&gparams[patch],nglx,ngly,nglz,xs,ys,zs,xe,ye,ze,gcw);

   // Create the patch object
   f_create_patch_data_(&data[patch],u0[patch],u[patch],aux[patch],gparams[patch]);

   filled[patch] = 1;
}


void * LevelContainer::getPtr(int patch)
{
   if ( patch<0 || patch>=N )
     {
       return NULL;
     }
   else if ( filled[patch] != 1 )
     {
       return NULL;
     }
   else
     {
       return data[patch];
     }
}


LevelContainer::~LevelContainer()
{
   for(int i=0; i<N; i++) 
     {
       if ( filled[i] == 1 ) 
	 {
	   // Delete the var_array on the fortran side
	   f_delete_patch_data_(data[i]);
	   // Set all pointers to null for safety
	   u[i] = NULL;
	   u0[i] = NULL;
	   aux[i] = NULL;
	   data[i] = NULL;
	   gparams[i] = NULL;
	 }
   }

   // Delete the internal variables
   delete[] u;
   delete[] u0;
   delete[] aux;
   delete[] data;
   delete[] gparams;
   delete[] filled;
   N = 0;
}

