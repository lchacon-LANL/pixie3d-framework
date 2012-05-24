#ifndef included_Pixie3dPreconditionerParameters
#define included_Pixie3dPreconditionerParameters


#include "SAMRAI/tbox/Database.h"
#include "SAMRAI/tbox/Pointer.h"
#include "SAMRAI/hier/PatchHierarchy.h"

#include "preconditioner_base/PreconditionerParameters.h"
#include "source/ParameterBase.h"
#include "interpolation/RefinementBoundaryInterpolation.h"

namespace SAMRAI{
namespace Pixie3d{

class Pixie3dPreconditionerParameters: public SAMRSolvers::PreconditionerParameters
{
public:
   Pixie3dPreconditionerParameters();

   Pixie3dPreconditionerParameters( const tbox::Pointer<tbox::Database>& database);

   ~Pixie3dPreconditionerParameters();

   /**
   *  Database object which needs to be initialized specific to the preconditioner.
   *  Documentation for parameters required by each preconditioner can be found in the
   *  documentation for the preconditioner.
   */
   tbox::Pointer<tbox::Database> d_db;

   /**
   * Pointer to the patch hierarchy object used by the preconditioner.
   */
   tbox::Pointer<hier::PatchHierarchy> d_hierarchy;

   /**
   * Pointer to the RefinementBoundaryInterpolation object used by the preconditioner to initialize and set coarse-fine boundary values.
   */
   SAMRAI::RefinementBoundaryInterpolation *d_cf_interpolant;      // object storing refinement boundary geometry and ghost value info.
   
protected:

private:

};

}
}
#endif
