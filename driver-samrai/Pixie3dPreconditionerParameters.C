#include "Pixie3dPreconditionerParameters.h"

namespace SAMRAI {
namespace Pixie3d {

Pixie3dPreconditionerParameters::Pixie3dPreconditionerParameters()
{
   d_db.reset();
   d_hierarchy.reset();
   d_cf_interpolant.reset();
}

Pixie3dPreconditionerParameters::Pixie3dPreconditionerParameters(const boost::shared_ptr<tbox::Database> &db)
   :d_db(db)
{
   d_hierarchy.reset();
   d_cf_interpolant.reset();
}

Pixie3dPreconditionerParameters::~Pixie3dPreconditionerParameters()
{
   d_cf_interpolant.reset();
}

}
}
