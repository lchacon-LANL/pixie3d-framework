#include <iostream>
#include "LevelContainer.h"
#include "PatchContainer.h"
#include "SAMRAI/geom/CartesianGridGeometry.h"
#include "SAMRAI/geom/CartesianPatchGeometry.h"


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



/********************************************************************
* Constructor                                                       *
********************************************************************/
LevelContainer::LevelContainer( tbox::Pointer<hier::PatchHierarchy> hierarchy, 
    tbox::Pointer< hier::PatchLevel> level,
    int n_var_in, int *u0_id_in, int *u_id_in, int n_auxs_in, 
    int *auxs_id_in, int n_auxv_in, int *auxv_id_in )
{
    // Copy the inputs
    d_hierarchy = hierarchy;
    d_level = level;
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
    // Create the patch objects
    data.resize(0);
    data.reserve(level->getNumberOfPatches());
    int index = 0;
    for (hier::PatchLevel::Iterator p(level); p; p++) {
        tbox::Pointer<hier::Patch> patch = *p;
        int local_id = patch->getLocalId().getValue();
        patch_map.insert( std::pair<int,int>(local_id,index) );
        data[index] = new PatchContainer(d_hierarchy,patch,n_var,u0_id,u_id,n_auxs,auxs_id,n_auxv,auxv_id);
        index++;
    }
}



void *LevelContainer::getPtr(tbox::Pointer<hier::Patch> patch)
{
    if ( patch.isNull() )
        return NULL;
    int local_id = patch->getLocalId().getValue();
    int index = patch_map.find(local_id)->second;
    if ( data[index] == NULL )
       return NULL;
    //delete data[index];
    //data[index] = new PatchContainer::PatchContainer(d_hierarchy,patch,n_var,u0_id,u_id,n_auxs,auxs_id,n_auxv,auxv_id);  
    void *data_ptr = data[index]->getPtr();
    return data_ptr;
}


LevelContainer::~LevelContainer()
{
    // Delete the patch data
    for (size_t i=0; i<data.size(); i++) {
        if ( data[i] != NULL )
            delete data[i];
    }
    data.clear();
    patch_map.clear();
    // Delete the internal variables
    delete [] u0_id;
    delete [] u_id;
    delete [] auxs_id;
    delete [] auxv_id;
}


