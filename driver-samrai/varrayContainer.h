#ifndef included_varraycontainer
#define included_varraycontainer

#include <vector>
#include "SAMRAI/tbox/Pointer.h"
#include "SAMRAI/hier/Patch.h"

#ifndef LACKS_NAMESPACE
using namespace SAMRAI;
#endif

class varrayContainer{
public:
   varrayContainer(tbox::Pointer<hier::Patch> patch, int n_var, int *u_id);
   ~varrayContainer();
   
   void *getPtr(){return data;}

private:
   varrayContainer();
   void CreatePatchFortran( int xs, int ys, int zs, int xe, int ye, 
      int ze, int gcw, int n_var, double **u_ptr );
   void *data;
};

#endif