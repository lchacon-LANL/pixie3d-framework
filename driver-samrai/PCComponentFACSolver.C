#include "PCComponentFACSolver.h"

namespace SAMRAI{
namespace Pixie3d{

PCComponentFACSolver ::PCComponentFACSolver( SAMRSolvers::MultilevelSolverParameters *parameters ):SAMRSolvers::CellFACPreconditioner(parameters)
{
  initializeLevelSolvers(parameters->d_db);
  
  if(d_operator!=NULL)
    {
      SAMRSolvers::MultilevelOperator *mlOperator = dynamic_cast<SAMRSolvers::MultilevelOperator *>(d_operator);
#ifdef DEBUG_CHECK_ASSERTIONS
      TBOX_ASSERT(mlOperator!=NULL);
#endif
      registerOperator(mlOperator);
    }
}

PCComponentFACSolver ::~PCComponentFACSolver()
{
}

void
PCComponentFACSolver::initializeLevelSolvers( tbox::Pointer<tbox::Database> db )
{
  tbox::pout << "Reached here" << std::endl;
}
  
}
}


  
