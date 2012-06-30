#include "PatchContainer.h"
#include <iostream>
#include "SAMRAI/geom/CartesianGridGeometry.h"
#include "SAMRAI/geom/CartesianPatchGeometry.h"

extern "C" {
#include "assert.h"
}

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


// Function to get which boundaries of a patch touch a periodic boundary
static void touchesPeriodicBoundary( SAMRAI::tbox::Pointer<SAMRAI::hier::Patch >& patch, 
    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel>& level, bool* periodic )
{
    const short int dim = level->getDim().getValue();
    const SAMRAI::tbox::Array<SAMRAI::hier::BoxContainer> domain_array = level->getPhysicalDomainArray();
    if ( domain_array.size() != 1 ) 
        TBOX_ERROR("Only 1 domain box is supported");
    const SAMRAI::hier::Index ifirst_global = domain_array[0].getBoundingBox().lower();
    const SAMRAI::hier::Index ilast_global  = domain_array[0].getBoundingBox().upper();
    const SAMRAI::hier::IntVector ratio = level->getRatioToLevelZero();
    const SAMRAI::hier::IntVector shift = level->getGridGeometry()->getPeriodicShift(ratio);
    const SAMRAI::hier::Index ifirst = patch->getBox().lower();
    const SAMRAI::hier::Index ilast  = patch->getBox().upper();
    for (int i=0; i<dim; i++) {
        periodic[2*i+0] = false;
        periodic[2*i+1] = false;
        if ( shift[i]!=0 ) {
            if ( ifirst(i)==0 )
                periodic[2*i+0] = true;
            if ( ilast(i)==ilast_global(i) )
                periodic[2*i+1] = true;
        } 
    }
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
    tmp_mem = NULL;
}


// Constructor to create and initialize the patch container object
PatchContainer::PatchContainer(tbox::Pointer<hier::PatchHierarchy> d_hierarchy, tbox::Pointer<hier::Patch> &patch, 
    int n_var, int *u0_id, int *u_id, int n_auxs, int *auxs_id, int n_auxv, int *auxv_id)
{
    tbox::Dimension dim(patch->getDim());
    tmp_mem = NULL;
    assert(!patch.isNull());
    int gcw = 1;
    hier::IntVector gcwc(dim,0);
    double *u_ptr[n_var], *u0_ptr[n_var], *auxs_ptr[n_auxs], *auxv_ptr[n_auxv];
    tbox::Pointer< pdat::CellData<double> > tmp;
    // Get the dimensions of the domain
    tbox::Pointer<geom::CartesianGridGeometry> grid_geometry = d_hierarchy->getGridGeometry();
    const double *lowerX = grid_geometry->getXLower();
    const double *upperX = grid_geometry->getXUpper();
    double lower[10], upper[10];
    for (int i=0; i<dim.getValue(); i++) {
        lower[i] = lowerX[i];
        upper[i] = upperX[i];
    }
    if ( grid_geometry->getNumberBlocks() != 1 )
        TBOX_ERROR("Multiblock domains are not supported");
    const SAMRAI::tbox::Array<SAMRAI::hier::BoxContainer> domain_array = d_hierarchy->getPatchLevel(0)->getPhysicalDomainArray();
    if ( domain_array.size() != 1 ) 
        TBOX_ERROR("Only 1 domain box is supported");
    const SAMRAI::hier::Box physicalDomain = domain_array[0].getBoundingBox();
    // Get the size of the patch
    const hier::Index ifirst = patch->getBox().lower();
    const hier::Index ilast  = patch->getBox().upper();
    int xs = ifirst(0)+1;
    int ys = ifirst(1)+1;
    int zs = ifirst(2)+1;
    int xe = ilast(0)+1;
    int ye = ilast(1)+1;
    int ze = ilast(2)+1;
    tbox::Pointer<hier::PatchGeometry> patch_geometry = patch->getPatchGeometry();
    const hier::IntVector ratio = patch_geometry->getRatio();
    int nbox[10];
    for (int i=0; i<10; i++)
        nbox[i] = 0;
    for (int i=0; i<dim.getValue(); i++)
        nbox[i] = physicalDomain.numberCells(i)*ratio(i);

    // Check that the patch size is >= 3x3
    // There is an error with 2x2 patches when creating the equlibrium.
    // Grid has to support second-order extrapolation to boundary, which requires at least 3 points in domain
    for (int i=0; i<dim.getValue(); i++) {
        int n = ilast(i)-ifirst(i)+1;
        if ( n<3 && n!=nbox[i] )
            TBOX_ERROR("Patches must be at least 3x3");
    }

    // Get the pointers to u_0, u_n
    for (int i=0; i<n_var; i++)  {
        // Get the pointers to u_0
        tmp = patch->getPatchData(u0_id[i]);
        if ( tmp.isNull() ) {
            // u0[i] is missing
            if ( patch->inHierarchy() )
                TBOX_ERROR("u0 is missing in patch"); 
            u0_ptr[i] = NULL;
        } else {
            u0_ptr[i] = tmp->getPointer();
            // Check ghost cell width
            gcwc = tmp->getGhostCellWidth();
            for (int j=0; j<dim.getValue(); j++) {
    	        if ( gcw!=gcwc(j) ) 
                    TBOX_ERROR("Ghost Cell width must be 1"); 
    	    }
        }
        // Get the pointers to u_n
        tmp = patch->getPatchData(u_id[i]);
        if ( tmp.isNull() ) {
            // u[i] is missing
            if ( patch->inHierarchy() )
                TBOX_ERROR("u is missing in patch"); 
            u_ptr[i] = NULL;
        } else {
            u_ptr[i] = tmp->getPointer();
            // Check ghost cell width
            gcwc = tmp->getGhostCellWidth();
            for (int j=0; j<dim.getValue(); j++) {
                if ( gcw!=gcwc(j) ) 
                    TBOX_ERROR("Ghost Cell width must be 1"); 
	        }
        }
    }

    // Get the pointers to the auxillary variables
    for (int i=0; i<n_auxs; i++) {
        // Get the pointers to scalar variables
        tmp = patch->getPatchData(auxs_id[i]);
        if ( tmp.isNull() ) {
            // auxs[i] is missing
            if ( patch->inHierarchy() )
                TBOX_ERROR("auxs is missing in patch"); 
            auxs_ptr[i] = NULL;
        } else {
            auxs_ptr[i] = tmp->getPointer();
            // Check ghost cell width
            gcwc = tmp->getGhostCellWidth();
            for (int j=0; j<dim.getValue(); j++) {
                if ( gcw!=gcwc(j) ) 
                    TBOX_ERROR("Ghost Cell width must be 1"); 
            }
        }
    }
    for (int i=0; i<n_auxv; i++) {
        // Get the pointers to vector variables
        tmp = patch->getPatchData(auxv_id[i]);
        if ( tmp.isNull() ) {
            // auxs[i] is missing
            if ( patch->inHierarchy() )
                TBOX_ERROR("auxv is missing in patch"); 
            auxv_ptr[i] = NULL;
        } else {
            auxv_ptr[i] = tmp->getPointer();
            // Check ghost cell width
            gcwc = tmp->getGhostCellWidth();
            for (int j=0; j<dim.getValue(); j++) {
                if ( gcw!=gcwc(j) ) 
	                TBOX_ERROR("Ghost Cell width must be 1"); 
            }
        }
    }

    // Check the patch range
    if ( xs>=1 && ys>=1 && zs>=1 && xe<=nbox[0] && ye<=nbox[1] && ze<=nbox[2] ) {
        // The patch is within the logical domain
    } else if ( patch->inHierarchy() ) {
        // The patch is in the hierarchy, no ranges should be outside the logical domain
        TBOX_ERROR("Error, patch is outside domain"); 
    } else {
        // A temporary patch is outside the logical domain
        if ( (xe-xs+1) > nbox[0] ) {
            // Temporary patch has more cells than domain, extend domain
            nbox[0] *= 2;
            upper[0] += upper[0]-lower[0];
        }
        if ( (ye-ys+1) > nbox[1] ) {
            // Temporary patch has more cells than domain, extend domain
            nbox[1] *= 2;
            upper[1] += upper[1]-lower[1];
        }
        if ( (ze-zs+1) > nbox[2] ) {
            // Temporary patch has more cells than domain, extend domain
            nbox[2] *= 2;
            upper[2] += upper[2]-lower[2];
        }
        // Shift the domain as necessary to cover the patch
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel> level = d_hierarchy->getPatchLevel(patch->getPatchLevelNumber());
        bool is_periodic[6];
        touchesPeriodicBoundary( patch, level, is_periodic );
        tbox::Pointer<hier::PatchGeometry> PatchGeom = patch->getPatchGeometry();
        int shift[10];
        for (int i=0; i<10; i++)
            shift[i] = 0;
        if ( xs < 1 )
            shift[0] = 1-xs;
    	if ( xe > nbox[0] )
            shift[0] = nbox[0]-xe;
        if ( ys < 1 )
            shift[1] = 1-ys;
    	if ( ye > nbox[1] )
            shift[1] = nbox[1]-ye;
        if ( zs < 1 )
            shift[2] = 1-zs;
    	if ( ze > nbox[2] )
            shift[2] = nbox[2]-ze;
        for (int i=0; i<dim.getValue(); i++) {
            if ( shift[i]==0 )
                continue;
            if ( !is_periodic[2*i+0] || !is_periodic[2*i+1] )
                TBOX_ERROR("Error, patch is outside physical domain"); 
            if ( i==0 ) {
                xs += shift[i];
                xe += shift[i];
            } else if ( i==1 ) {
                ys += shift[i];
                ye += shift[i];
            } else if ( i==2 ) {
                zs += shift[i];
                ze += shift[i];
            }
            lower[i] -= (upper[i]-lower[i])*((double) shift[i])/((double) nbox[i]);
            upper[i] -= (upper[i]-lower[i])*((double) shift[i])/((double) nbox[i]);
        }
        // Check that the shifted patch is within the domain (1:nbox)
        if ( xs<1 || ys<1 || zs<1 || xe>nbox[0] || ye>nbox[1] || ze>nbox[2] )
            TBOX_ERROR("Unable to create a valid temporary patch"); 
    }

    // Create the patch object
    CreatePatchFortran(lower,upper,nbox,xs,ys,zs,xe,ye,ze,gcw,n_var,u0_ptr,u_ptr,n_auxs,auxs_ptr,n_auxv,auxv_ptr);

}


// Function to initialize the fortran side of the patch
void PatchContainer::CreatePatchFortran( const double *lowerCoordinates,const double *upperCoordinates, const int *ng, 
    int xs, int ys, int zs, int xe, int ye, int ze, int gcw, int n_var, double **u0_ptr, double **u_ptr, int n_auxs, 
    double **auxs_ptr, int n_auxv, double **auxv_ptr)
{
    if ( xs<1 || ys<1 || zs<1 || xe>ng[0] || ye>ng[1] || ze>ng[2] )
        TBOX_ERROR("Bad index"); 
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

    // Allocate temporary memory of u0 if necessary
    bool use_tmp_mem = false;
    tmp_mem = NULL;
    for(int i=0; i<n_var; i++) {
        if ( u0_ptr[i]==NULL )
            use_tmp_mem = true;
    }
    if ( use_tmp_mem == true ) {
        int N = (xeg-xsg+1)*(yeg-ysg+1)*(zeg-zsg+1);
        tmp_mem = new double[n_var*N];
        for(int i=0; i<n_var; i++) {
            if ( u0_ptr[i]==NULL )
                u0_ptr[i] = &tmp_mem[i*N];
        }
    }

    // Create var_array for u, u0
    f_create_var_array_(&u,n_var);
    f_create_var_array_(&u0,n_var);
    for (int i=0; i<n_var; i++) {
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
    // Delete the temporary memory
    if ( tmp_mem != NULL )
        delete [] tmp_mem;
    // Set all pointers to null for safety
    u = NULL;
    u0 = NULL;
    aux = NULL;
    data = NULL;
    gparams = NULL;
    tmp_mem = NULL;
}



