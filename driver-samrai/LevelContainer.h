#ifndef included_levelcontainer
#define included_levelcontainer

#include <vector>
#include "tbox/Array.h"
#include "tbox/Pointer.h"
#include "CellData.h"
#include "CellVariable.h"
#include "ComponentSelector.h"
#include "PatchHierarchy.h"
#include "VariableContext.h"
#include "PatchContainer.h"

#ifndef LACKS_NAMESPACE
using namespace SAMRAI;
#endif



class LevelContainer{
public:
    LevelContainer( const int n, tbox::Pointer< hier::PatchHierarchy<NDIM> > hierarchy,
        int n_var, int *u0_id, int *u_id, int n_auxs, int *auxs_id, int n_auxv, int *auxv_id);
    ~LevelContainer();
    void CreatePatch(int patch_id, tbox::Pointer< hier::Patch<NDIM> > &patch );
    void *getPtr(int patch);
    void *getPtr2(int patch);

private:
    LevelContainer();
    int N;
    tbox::Pointer< hier::PatchHierarchy<NDIM> > d_hierarchy;
    int n_var, *u0_id, *u_id, n_auxs, *auxs_id, n_auxv, *auxv_id;
    PatchContainer **data;
    tbox::Pointer< hier::Patch<NDIM> > *patch_ptr;
};

#endif


