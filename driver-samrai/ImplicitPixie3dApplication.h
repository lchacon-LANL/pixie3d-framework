#ifndef included_ImplicitPixie3dApplication
#define included_ImplicitPixie3dApplication

#include "SAMRAI/algs/ImplicitEquationStrategy.h"
#include "SAMRAI/tbox/Serializable.h"
#include "SAMRAI/solv/SNESAbstractFunctions.h"
#include "pixie3dApplication.h"
#include "ImplicitPixie3dApplicationParameters.h"
#include "Pixie3dPreconditioner.h"

namespace SAMRAI{
namespace Pixie3d{
    
class ImplicitPixie3dApplication:
    public SAMRAI::pixie3dApplication,
    public SAMRAI::algs::ImplicitEquationStrategy, 
    public SAMRAI::solv::SNESAbstractFunctions,
    public SAMRAI::tbox::Serializable
{
public: 

  /*!
   * DEfault constructor
   */
  ImplicitPixie3dApplication();

  /*!
   * Constructor that takes a parameter list.  Calls initialize.
   */
  ImplicitPixie3dApplication(ImplicitPixie3dApplicationParameters *parameters);

  /*!
   * Destructor
   */
  ~ImplicitPixie3dApplication();  
  
  /**
   * Initialize application using specified parameters.
   */
  void initialize(ImplicitPixie3dApplicationParameters *parameters);
  
  /*
   * Interface functions overloaded from solv::SNESAbstractFunctions.
   */
  int evaluateNonlinearFunction(Vec xcur, Vec fcur);
  
  /*
   * Interface functions overloaded from solv::SNESAbstractFunctions.
   */
  int setupPreconditioner(Vec x);

  /*
   * Interface functions overloaded from solv::SNESAbstractFunctions.
   */
   int applyPreconditioner(Vec r, Vec z);
  
  /*
   * Interface functions overloaded from solv::SNESAbstractFunctions.
   */
  int evaluateJacobian(Vec x)
  {
    (void) x;
    return 0;
  }
  
  /*
   * Interface functions overloaded from solv::SNESAbstractFunctions.
   */
  int jacobianTimesVector(Vec v, Vec Jv)
  {
    (void) v;
    (void) Jv;
    return 0;
  }

  /*!
   * User-supplied nonlinear function evaluation.  Returns 0 if successful.
   * Arguments:
   *
   * @param xcur (IN) the current iterate for the nonlinear system
   * @param fcur (OUT) current function value
   *
   * IMPORTANT: This function must not modify xcur.
   *
   * Function overloaded from SAMRAI::solv::SNESAbstractFunctions.
   */
  int evaluateNonlinearFunction(tbox::Pointer< solv::SAMRAIVectorReal<double> > x, tbox::Pointer< solv::SAMRAIVectorReal<double> > f);  
  
  /*!
   * User-supplied preconditioner setup function.  The setup
   * function is called to provide matrix data for the subsequent
   * call(s) to applyPreconditioner().  The integer return value
   * is a flag indicating success if 0 is returned, and failure
   * otherwise.  Together setupPreconditioner() and applyPreconditioner()
   * form a right preconditoner for the PETSc Krylov solver.
   * Returns 0 if successful.
   *
   *
   * Function overloaded from SAMRAI::solv::SNESAbstractFunctions.
   */
  int setupPreconditioner(tbox::Pointer< solv::SAMRAIVectorReal<double> > x);
  
   /*!
    * User-supplied preconditioner solve function.  This function must
    * solve \f$M z = r\f$, where \f$M\f$ is the right preconditioner
    * matrix formed by setupPreconditioner(). The integer return value
    * is a flag indicating success if 0 is returned, and failure otherwise.
    * Arguments:
    *
    * @param r (IN) right-hand side of preconditioning system
    * @param z (OUT) result of applying preconditioner to right-hand side
    *
    * IMPORTANT: This function must not modify r.
    *
    * Function overloaded from SAMRAI::solv::SNESAbstractFunctions.
    */
  int applyPreconditioner(tbox::Pointer< solv::SAMRAIVectorReal<double> > r,  tbox::Pointer< solv::SAMRAIVectorReal<double> > z);
  
  /**
   * Set the nonlinear solution vector so that the new solution data is
   * solved for when the nonlinear solver advances the solution.
   *
   * When assertion checking is active, passing in a null pointer for
   * the solution vector will result in an unrecoverable exception.
   *
   * Function overloaded from algs::ImplicitEquationStrategy.
   *
   */
  void setupSolutionVector(tbox::Pointer< solv::SAMRAIVectorReal<double> > solution);
  
  /**
   * Return time increment for advancing the solution at the first timestep.
   *
   * Function overloaded from algs::ImplicitEquationStrategy.
   */
  double getInitialDt();
  
  /**
   * Return the next time increment through which to advance the solution.
   * The good_solution is the value returned by a call to checkNewSolution(),
   * which determines whether the computed solution is acceptable or not.
   * The integer solver_retcode is the return code generated by the
   * nonlinear solver.   This value must be interpreted in a manner
    * consistant with the solver in use.
    *
    * Function overloaded from algs::ImplicitEquationStrategy.
    */
  double getNextDt(const bool good_solution, const int solver_retcode);
  
  /**
   * Set the initial guess for the time advanced solution at the start
   * of the nonlinear iteration.  The boolean argument first_step
   * indicates whether we are at the first step on the current hierarchy
   * configuration.  This is true when the hierarchy is constructed
   * initially and after regridding.  In these cases, setting the initial
   * iterate using extrapolation, for example, may not be possible.
   *
   * Function overloaded from algs::ImplicitEquationStrategy.
   */
  void setInitialGuess(const bool first_step, const double current_time,
		       const double current_dt, const double old_dt);
  
  /**
   * Check the computed solution and return true if it is acceptable;
   * otherwise return false.  The integer solver_retcode is the return
   * code generated by the nonlinear solver.  This value must be
   * interpreted in a manner consistent with the solver in use.
   *
   * Function overloaded from algs::ImplicitEquationStrategy.
   */
  bool checkNewSolution(const int solver_retcode);
  
  /**
   * Update solution storage and dependent quantities after computing an
   * acceptable time advanced solution.   The new_time value is the new
   * solution time.
   *
   * Function overloaded from algs::ImplicitEquationStrategy.
   */
   void updateSolution(const double new_time);
   
   /**
    * Write data members to given data base for restart.
    *
    * When assertion checking is enabled, passing in a null pointer for 
    * the database will result in an unrecoverable exception.
    *
    * Overloaded from tbox::Serializable.
    */
   void putToDatabase(tbox::Pointer<tbox::Database> db);

   /**
    * Reset cached information that depends on the hierarchy configuration.  
    *
    * Function overloaded from mesh::StandardTagAndInitStrategy.
    */
   void resetHierarchyConfiguration(
           const tbox::Pointer<hier::PatchHierarchy> hierarchy,
           const int coarsest_level,
           const int finest_level );

   /**
    * create the preconditioner
    */
   void createPreconditioner( void );

   /**
    * Destroy the preconditioner.
    */
   void destroyPreconditioner(void);

 private:
    
   Pixie3dPreconditionerParameters *
   createPreconditionerParameters( tbox::Pointer<tbox::Database> &db );

   /*
    * The nonlinear solution process requires a solution vector; we cache
    * a pointer to it here.
    */
   tbox::Pointer< solv::SAMRAIVectorReal<double> > d_newSolutionVector;

   /*
    * We are hard coding the BDF2 solution scheme for now as the time integration
    * scheme. This vector will store the current computed solution
    */
   tbox::Pointer< solv::SAMRAIVectorReal<double> > d_currentSolutionVector;

   /*
    * We are hard coding the BDF2 solution scheme for now as the time integration
    * scheme. This vector will store the solution at the previous time level
    */
   tbox::Pointer< solv::SAMRAIVectorReal<double> > d_previousSolutionVector;

   /*
    * We are hard coding the BDF2 solution scheme for now as the time integration
    * scheme. This vector will be used to store intermediate quantities
    */
   tbox::Pointer< solv::SAMRAIVectorReal<double> > d_scratchVector;

   /**
    * Pointer to pixie preconditioner
    */
   tbox::Pointer< Pixie3dPreconditioner > d_preconditioner;
   
   /**
    * Pointer to database with preconditioner parameters
    */
   tbox::Pointer<tbox::Database> d_pc_db;

   /*
    * Current solution time and time increment used in the solution process.
    * New time is current time + current dt.
    */
   bool d_first_step;
   bool d_first_regrid;

   /**
    * bool flag to ensure current, previous, and new vectors are only cloned once
    */
   bool d_vectorsCloned;
   
   double d_current_time;
   double d_new_time;

   /**
    * The initial time step
    */
   double d_initial_dt;

   /**
    * The current time step
    */
   double d_current_dt;

   /**
    * The previous time step
    */
   double d_old_dt;

   /**
    * The maximum user specified time step
    */
   double d_max_timestep;
   
};

}
} 
#endif
