#include "PatchContainer.h"
#include <iostream>
#include "CartesianGridGeometry.h"
#include "CartesianPatchGeometry.h"


extern "C"{
   void f_create_var_array_(void **p_data, int&);
   void f_create_aux_array_(void **p_data, int&, int&);
   void f_create_grid_mg_def_(void **p_data, const double *, const double *, int &, int &, int &, int &, int &, int &, int &, int &, int &, int &);
   void f_create_patch_data_(void **p_data, void *, void *, void *, void *);
   void creategridstructures_(void *p_data);
   void formequilibrium_(void *p_data);
   void initialize_u_n_(void *p_data);
   void forminitialcondition_(void *p_data);
   void setupvarinitseq_(void *p_data, int*);
   void f_delete_patch_data_(void *p_data);
   void fill_var_array_(void *p_data, int &, double *, int &, int &, int &, int &, int &, int &);
   void fill_aux_array_var_(void *p_data, int &, double *, int &, int &, int &, int &, int &, int &);
   void fill_aux_array_vec_(void *p_data, int &, double *, int &, int &, int &, int &, int &, int &);
}


// Empty constructor
PatchContainer::PatchContainer()
{
    u = NULL;
    u0 = NULL;
    aux = NULL;
    data = NULL;
    gparams = NULL;
    data = NULL;
}


// Constructor to create and initialize the patch container object
PatchContainer::PatchContainer(tbox::Pointer< hier::PatchHierarchy<NDIM> > d_hierarchy, tbox::Pointer< hier::Patch<NDIM> > &patch, 
    int n_var, int *u0_id, int *u_id, int n_auxs, int *auxs_id, int n_auxv, int *auxv_id)
{
    int gcw = 1;
    hier::IntVector<NDIM> gcwc;
    double *u_ptr[n_var], *u0_ptr [n_var], *auxs_ptr[n_auxs], *auxv_ptr[n_auxv];
    tbox::Pointer< pdat::CellData<NDIM,double> > tmp;
    // Get the dimensions of the domain
    tbox::Pointer<geom::CartesianGridGeometry<NDIM> > grid_geometry = d_hierarchy->getGridGeometry();
    const double *lower = grid_geometry->getXLower();
    const double *upper = grid_geometry->getXUpper();
    const SAMRAI::hier::BoxArray<NDIM> &physicalDomain = grid_geometry->getPhysicalDomain() ;
    int nbox[NDIM];
    for (int i=0; i<NDIM; i++)
        nbox[i] = physicalDomain[0].numberCells(i);
    // Get the size of the patch
    const hier::Index<NDIM> ifirst = patch->getBox().lower();
    const hier::Index<NDIM> ilast  = patch->getBox().upper();
    int xs = ifirst(0)+1;
    int ys = ifirst(1)+1;
    int zs = ifirst(2)+1;
    int xe = ilast(0)+1;
    int ye = ilast(1)+1;
    int ze = ilast(2)+1;
    tbox::Pointer<hier::PatchGeometry<NDIM> > patch_geometry = patch->getPatchGeometry();
    const hier::IntVector<NDIM> ratio = patch_geometry->getRatio();
    int ng[NDIM];
    for (int i=0; i<NDIM; i++)
        ng[i] = nbox[i]*ratio(i);

    // Get the pointers to u_0, u_n
    for (int i=0; i<n_var; i++)  {
        // Get the pointers to u_0
        tmp = patch->getPatchData(u0_id[i]);
        u0_ptr[i] = tmp->getPointer();
        // Check ghost cell width
        gcwc = tmp->getGhostCellWidth();
        for (int j=0; j<NDIM; j++) {
	        if ( gcw!=gcwc(j) ) 
                TBOX_ERROR("Ghost Cell width must be 1"); 
	    }
        // Get the pointers to u_n
        tmp = patch->getPatchData(u_id[i]);
        u_ptr[i] = tmp->getPointer();
        // Check ghost cell width
        gcwc = tmp->getGhostCellWidth();
        for (int j=0; j<NDIM; j++) {
            if ( gcw!=gcwc(j) ) 
                TBOX_ERROR("Ghost Cell width must be 1"); 
	    }
    }

   // Get the pointers to the auxillary variables
   for (int i=0; i<n_auxs; i++) {
       // Get the pointers to scalar variables
       tmp = patch->getPatchData(auxs_id[i]);
       auxs_ptr[i] = tmp->getPointer();
       // Check ghost cell width
       gcwc = tmp->getGhostCellWidth();
       for (int j=0; j<NDIM; j++) {
            if ( gcw!=gcwc(j) ) 
                TBOX_ERROR("Ghost Cell width must be 1"); 
        }
    }
    for (int i=0; i<n_auxv; i++) {
        // Get the pointers to vector variables
        tmp = patch->getPatchData(auxv_id[i]);
        auxv_ptr[i] = tmp->getPointer();
        // Check ghost cell width
        gcwc = tmp->getGhostCellWidth();
        for (int j=0; j<NDIM; j++) {
            if ( gcw!=gcwc(j) ) 
	       TBOX_ERROR("Ghost Cell width must be 1"); 
        }
   }

   // Create the patch object
   CreatePatchFortran(lower,upper,ng,xs,ys,zs,xe,ye,ze,gcw,n_var,u0_ptr,u_ptr,n_auxs,auxs_ptr,n_auxv,auxv_ptr);
   //f_delete_patch_data_(data);
   //CreatePatchFortran(lower,upper,ng,xs,ys,zs,xe,ye,ze,gcw,n_var,u0_ptr,u_ptr,n_auxs,auxs_ptr,n_auxv,auxv_ptr);

}


// Function to initialize the fortran side of the patch
void PatchContainer::CreatePatchFortran( const double *lowerCoordinates,const double *upperCoordinates, const int *ng, 
					 int xs, int ys, int zs, int xe, int ye, int ze, int gcw, 
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
   
    // Set global number of points for the level
    int nglx = ng[0];
    int ngly = ng[1];
    int nglz = ng[2];

    // Create var_array for u, u0
    f_create_var_array_(&u,n_var);
    f_create_var_array_(&u0,n_var);
    for(int i=0; i<n_var; i++) {
        fill_var_array_( u,i, u_ptr[i],xsgt,ysgt,zsgt,xegt,yegt,zegt);
        fill_var_array_(u0,i,u0_ptr[i],xsgt,ysgt,zsgt,xegt,yegt,zegt);
    }

    // Create aux_array auxillary variables
    f_create_aux_array_(&aux,n_auxs,n_auxv);
    for(int i=0; i<n_auxs; i++)
        fill_aux_array_var_(aux,i,auxs_ptr[i],xsgt,ysgt,zsgt,xegt,yegt,zegt);
    for(int i=0; i<n_auxv; i++)
        fill_aux_array_vec_(aux,i,auxv_ptr[i],xsgt,ysgt,zsgt,xegt,yegt,zegt);
   
    // Create grid_mg_def for the patch information
    f_create_grid_mg_def_(&gparams,lowerCoordinates, upperCoordinates, nglx,ngly,nglz,xs,ys,zs,xe,ye,ze,gcw);

    // Create the patch object
    f_create_patch_data_(&data,u0,u,aux,gparams);

    // Create the grid structures (moved from setInitialConditions)
    creategridstructures_(data);

    // Initilize the fortran data
    formequilibrium_(data);
    int it = 0;     // it = 0 for equilibrium setup, 1 for normal
    setupvarinitseq_(data,&it);

    // Initialize u_n
    initialize_u_n_(data);
}


// Deconstructor
PatchContainer::~PatchContainer()
{
    // Delete the var_array on the fortran side
    f_delete_patch_data_(data);
    // Set all pointers to null for safety
    u = NULL;
    u0 = NULL;
    aux = NULL;
    data = NULL;
    gparams = NULL;
}



