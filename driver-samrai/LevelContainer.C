#include "LevelContainer.h"
#include <iostream>
#include "PatchContainer.h"
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


#if NDIM != 3
#error Not programed for dimension other than 3
#endif


LevelContainer::LevelContainer()
{
   data = NULL;
   patch_ptr = NULL;
}

LevelContainer::LevelContainer(int n, tbox::Pointer< hier::PatchHierarchy<NDIM> > hierarchy,
      int n_var_in, int *u0_id_in, int *u_id_in, int n_auxs_in, int *auxs_id_in, int n_auxv_in, int *auxv_id_in )
{
    N = n;
    d_hierarchy = hierarchy;
    // Create space for a patch object for each patch (n patches)
    data = new PatchContainer::PatchContainer *[N];
    patch_ptr = new tbox::Pointer< hier::Patch<NDIM> > [N];
    for (int i=0; i<N; i++) {
        data[i] = NULL;
        patch_ptr[i].setNull();;
    }
    n_var = n_var_in;
    u_id = new int[n_var];
    u0_id = new int[n_var];
    for (int i=0; i<n_var; i++) {
        u_id[i] = u_id_in[i];
        u0_id[i] = u0_id_in[i];
    }
    n_auxs = n_auxs_in;
    auxs_id = new int[n_auxs];
    for (int i=0; i<n_auxs; i++) 
        auxs_id[i] = auxs_id_in[i];
    n_auxv = n_auxv_in;
    auxv_id = new int[n_auxv];
    for (int i=0; i<n_auxv; i++) 
        auxv_id[i] = auxv_id_in[i];
}


void LevelContainer::CreatePatch(int patch_id, tbox::Pointer< hier::Patch<NDIM> > &patch )
{
    patch_ptr[patch_id] = patch;
    data[patch_id] = new PatchContainer::PatchContainer(d_hierarchy,patch_ptr[patch_id],n_var,u0_id,u_id,n_auxs,auxs_id,n_auxv,auxv_id);
    //delete data[patch_id];
    //data[patch_id] = NULL;
    //data[patch_id] = new PatchContainer::PatchContainer(d_hierarchy,patch_ptr[patch_id],n_var,u0_id,u_id,n_auxs,auxs_id,n_auxv,auxv_id);
}



void * LevelContainer::getPtr(int patch)
{
    if ( patch<0 || patch>=N )
        return NULL;
    if ( data[patch] == NULL )
       return NULL;
    //delete data[patch];
    //data[patch] = new PatchContainer::PatchContainer(d_hierarchy,patch_ptr[patch],n_var,u0_id,u_id,n_auxs,auxs_id,n_auxv,auxv_id);  
    void *data_ptr = data[patch]->getPtr();
    return data_ptr;
}


LevelContainer::~LevelContainer()
{
    // Delete the patch data
    if ( data != NULL ) {
        for(int i=0; i<N; i++) {
            if ( data[i] != NULL )
                delete data[i];
        }
        delete [] data;
        delete [] patch_ptr;
        data = NULL;
    }
    // Delete the internal variables
    delete [] u0_id;
    delete [] u_id;
    delete [] auxs_id;
    delete [] auxv_id;
    N = 0;
}


