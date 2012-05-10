#ifndef included_Pixie3dPreconditioner
#define included_Pixie3dPreconditioner

#include "SAMRAI/pdat/FaceVariable.h"
#include "SAMRAI/xfer/CoarsenSchedule.h"
#include "SAMRAI/xfer/CoarsenAlgorithm.h"
#include "SAMRAI/tbox/List.h"

#include "preconditioner_base/PreconditionerStrategy.h"
#include "interpolation/RefinementBoundaryInterpolation.h"
#include "solvers/multilevel/MultilevelSolverFactory.h"
#include "solvers/MultilevelSolver.h"
#include "operators/MultilevelOperatorFactory.h"

#include "Pixie3dPreconditionerParameters.h"

#include <vector>

typedef SAMRAI::tbox::List< SAMRAI::tbox::Pointer<SAMRAI::xfer::CoarsenSchedule > > crsList;

namespace SAMRAI{

namespace SAMRSolvers{

class Pixie3dPreconditioner: public PreconditionerStrategy
{
public:
   Pixie3dPreconditioner(Pixie3dPreconditionerParameters *parameters);

   ~Pixie3dPreconditioner();
   
   int setupPreconditioner( PreconditionerParameters *parameters );
   
   int applyPreconditioner( tbox::Pointer< solv::SAMRAIVectorReal<double> > r,
                            tbox::Pointer< solv::SAMRAIVectorReal<double> > z );

   /**
    * Functions to read data from input and restart databases. If the
    * boolean flag is true, all data members are read from restart.
    * They can later be overwritten from values in the input file.
    * When the flag is false, all data values are set from those given
    * in input.
    *
    * If assertion checking is enabled, an unrecoverable exception
    * results if the database pointer is null.
   */
   void getFromInput( tbox::Pointer<tbox::Database> &db,
                      bool is_from_restart = false);
   
   void setRefinementBoundaryInterpolant(RefinementBoundaryInterpolation *cf_interpolant);

protected:

private:

   // private constuctor to prevent it being called
   Pixie3dPreconditioner();

   void preprocessPCApply( int r_id );

   void postprocessPCApply( int z_id );

   void interpolateVariable( const int src_id, 
                             const int dest_id, 
                             RefinementBoundaryInterpolation::InterpolationScheme tangential_interp_scheme,
                             RefinementBoundaryInterpolation::InterpolationScheme normal_interp_scheme );

   void initializeSolvers(tbox::Pointer<tbox::Database> &db);

   void initializeOperators(tbox::Pointer<tbox::Database> &db);

   void coarsenVariable( const int var_id,
                         std::string coarsen_op_str);
   /*
    * Retrieve list of potentially usable sibling fill schedules for a level.
    */
   crsList* getCoarsenSchedules(int ln) { return( &d_coarsen_fill_schedules[ln] ); }

   /*
    * Retrieve list of potentially usable sibling fill schedules for a level.
    */
    tbox::Pointer<xfer::CoarsenSchedule > getCoarsenSchedule(int ln, xfer::CoarsenAlgorithm &);

   bool  d_preconditioner_print_flag;

   double d_dt;

   /*
    * Cached list of fill schedules, to avoid creating one for every
    * sibling fill.
    */
   std::vector< crsList > d_coarsen_fill_schedules;

   tbox::Pointer<hier::PatchHierarchy > d_hierarchy;

   RefinementBoundaryInterpolation *d_cf_interpolant;

};

}

}
#endif
