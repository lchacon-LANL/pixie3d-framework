#include "PCComponentFACSolver.h"

namespace SAMRAI{
namespace Pixie3d{

PCComponentFACSolver ::PCComponentFACSolver( boost::shared_ptr<SAMRSolvers::MultilevelSolverParameters> parameters ):
    SAMRSolvers::CellFACPreconditioner(parameters)
{
  initializeLevelSolvers(parameters->d_db);
  
  if(d_operator!=NULL)
    {
      boost::shared_ptr<SAMRSolvers::MultilevelOperator> mlOperator = 
        boost::dynamic_pointer_cast<SAMRSolvers::MultilevelOperator>(d_operator);
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
PCComponentFACSolver::initializeLevelSolvers( boost::shared_ptr<tbox::Database> db )
{
  tbox::pout << "Reached here" << std::endl;
}
  
}
}


  
