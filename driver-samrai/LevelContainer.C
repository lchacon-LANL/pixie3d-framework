#include <iostream>
#include "LevelContainer.h"
#include "PatchContainer.h"
#include "SAMRAI/geom/CartesianGridGeometry.h"
#include "SAMRAI/geom/CartesianPatchGeometry.h"





/********************************************************************
* Constructor                                                       *
********************************************************************/
LevelContainer::LevelContainer( boost::shared_ptr<hier::PatchHierarchy> hierarchy, 
    boost::shared_ptr< hier::PatchLevel> level,
    int n_var_in, int *u0_id_in, int *u_id_in, int n_auxs_in, 
    int *auxs_id_in, int n_auxv_in, int *auxv_id_in, int flux_id_in, int src_id_in )
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
    flux_id = flux_id_in;
    src_id = src_id_in;
    // Create the patch objects
    data.resize(0);
    data.reserve(level->getNumberOfPatches());
    int index = 0;
    for (hier::PatchLevel::Iterator p=level->begin(); p!=level->end(); p++) {
        boost::shared_ptr<hier::Patch> patch = *p;
        int local_id = patch->getLocalId().getValue();
        patch_map.insert( std::pair<int,int>(local_id,index) );
        data[index] = new PatchContainer(d_hierarchy,patch,n_var,u0_id,u_id,n_auxs,auxs_id,n_auxv,auxv_id,flux_id,src_id);
        index++;
    }
}



void *LevelContainer::getPtr(boost::shared_ptr<hier::Patch> patch)
{
    if ( patch.get()==NULL )
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


