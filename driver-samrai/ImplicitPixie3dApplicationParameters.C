#include "ImplicitPixie3dApplicationParameters.h"

namespace SAMRAI{
namespace Pixie3d{

ImplicitPixie3dApplicationParameters::ImplicitPixie3dApplicationParameters():pixie3dApplicationParameters()
{
}
  
ImplicitPixie3dApplicationParameters::ImplicitPixie3dApplicationParameters(tbox::Pointer<tbox::Database> database):pixie3dApplicationParameters(database)
{
}
  
ImplicitPixie3dApplicationParameters::~ImplicitPixie3dApplicationParameters()
{
}
  
}
}
