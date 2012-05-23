#include "Pixie3dPreconditionerParameters.h"

namespace SAMRAI {
namespace SAMRSolvers {

Pixie3dPreconditionerParameters::Pixie3dPreconditionerParameters()
{
   d_db.setNull();
   d_hierarchy.setNull();
   d_cf_interpolant          = NULL;
}

Pixie3dPreconditionerParameters::Pixie3dPreconditionerParameters(const tbox::Pointer<tbox::Database> &db)
   :d_db(db)
{
   d_hierarchy.setNull();
   d_cf_interpolant          = NULL;
}

Pixie3dPreconditionerParameters::~Pixie3dPreconditionerParameters()
{
   d_cf_interpolant          = NULL;
}

}
}
