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

#ifndef LACKS_NAMESPACE
using namespace SAMRAI;
#endif

class LevelContainer{
public:
   LevelContainer(const int n, int nx, int ny, int nz);
   ~LevelContainer();
   void CreatePatch(int patch_id, tbox::Pointer< hier::Patch<NDIM> > patch,
      int n_var, int *u0_id, int *u_id, int n_auxs, int *auxs_id, int n_auxv, int *auxv_id);  
   void *getPtr(int patch);

private:
   LevelContainer();
   void CreatePatchFortran( int patch, int xs, int ys, int zs, int xe, int ye, 
      int ze, int gcw, int n_var, double **u0_ptr, double **u_ptr, 
      int n_auxs, double **auxs_ptr, int n_auxv, double **auxv_ptr);
   void **data;
   void **u;
   void **u0;
   void **aux;
   void **gparams;
   int *filled;
   int N, nglx, ngly, nglz;
};

#endif
