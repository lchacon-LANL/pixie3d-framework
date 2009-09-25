#ifndef included_varraycontainer
#define included_varraycontainer

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

class varrayContainer{
public:
   varrayContainer(tbox::Pointer< hier::Patch<NDIM> > patch, int n_var, int *u_id);
   ~varrayContainer();
   
   void *getPtr(){return data;}

private:
   varrayContainer();
   void CreatePatchFortran( int xs, int ys, int zs, int xe, int ye, 
      int ze, int gcw, int n_var, double **u_ptr );
   void *data;
};

#endif
