#ifndef included_PCComponentFACSolver
#define included_PCComponentFACSolver

#include "solvers/CellFACPreconditioner.h"

namespace SAMRAI{
namespace Pixie3d{

class PCComponentFACSolver: public SAMRSolvers::CellFACPreconditioner
{

public:  

  PCComponentFACSolver( boost::shared_ptr<SAMRSolvers::MultilevelSolverParameters> parameters );

  ~PCComponentFACSolver();

protected:
  /**
   * initialization of the level smoothers/solvers with pixie specific smoothers
   */

  void initializeLevelSolvers( boost::shared_ptr<tbox::Database> db );
private:  

};
  
}
}


#endif
