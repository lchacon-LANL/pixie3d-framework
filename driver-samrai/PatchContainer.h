#ifndef included_patchcontainer
#define included_patchcontainer

#include <vector>
#include "SAMRAI/tbox/Array.h"
#include "SAMRAI/tbox/Pointer.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/pdat/CellVariable.h"
#include "SAMRAI/hier/ComponentSelector.h"
#include "SAMRAI/hier/PatchHierarchy.h"
#include "SAMRAI/hier/VariableContext.h"

#ifndef LACKS_NAMESPACE
using namespace SAMRAI;
#endif

class PatchContainer{
public:
    PatchContainer( tbox::Pointer<hier::PatchHierarchy> hierarchy, tbox::Pointer<hier::Patch> &patch, 
        int n_var, int *u0_id, int *u_id, int n_auxs, int *auxs_id, int n_auxv, int *auxv_id );  
    ~PatchContainer();
    void *getPtr() { return data; }

private:
    PatchContainer();
    void CreatePatchFortran( const double *lower, const double *upper, const int *ng, 
       int xs, int ys, int zs, int xe, int ye, int ze, int gcw, int n_var, double **u0_ptr, double **u_ptr, 
       int n_auxs, double **auxs_ptr, int n_auxv, double **auxv_ptr );
    void *data;
    void *u;
    void *u0;
    void *aux;
    void *gparams;
    double *tmp_mem;
};

#endif
