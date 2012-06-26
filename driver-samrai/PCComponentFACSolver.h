#ifndef included_PCComponentFACSolver
#define included_PCComponentFACSolver

#include "solvers/CellFACPreconditioner.h"

namespace SAMRAI{
namespace Pixie3d{

class PCComponentFACSolver: public SAMRSolvers::CellFACPreconditioner
{

public:  

  PCComponentFACSolver( SAMRSolvers::MultilevelSolverParameters *parameters );

  ~PCComponentFACSolver();

protected:
  /**
   * initialization of the level smoothers/solvers with pixie specific smoothers
   */

  void initializeLevelSolvers( tbox::Pointer<tbox::Database> db );
private:  

};
  
}
}


#endif
