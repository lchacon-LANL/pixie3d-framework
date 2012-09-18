#ifndef included_ImplicitPixie3dApplicationParameters
#define included_ImplicitPixie3dApplicationParameters

#include "pixie3dApplicationParameters.h"

namespace SAMRAI{
namespace Pixie3d{

class ImplicitPixie3dApplicationParameters: public pixie3dApplicationParameters
{
 public:
  ImplicitPixie3dApplicationParameters();
  ImplicitPixie3dApplicationParameters(tbox::Pointer<tbox::Database> database);
  ~ImplicitPixie3dApplicationParameters();    
};  

}
}
#endif
