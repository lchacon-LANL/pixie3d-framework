static char help[] = "Options:
	-snes_mf : use matrix-free Newton methods (default)\n\
	-user_precond: use user-defined preconditioner (-snes_mf)\n\
        -nmax <nmax>: max. iteration # of time steps\n";

/*
 *    Concepts: SNES matrix-free/finite-difference Jacobian methods
 *    Concepts: SNES user-provided preconditioner;
 *    Processors: n (parallel)
 *    
 */

#include "petscsnes.h"
#include "petscda.h"
#include "fortran.h"
#include <stdio.h>
#include <string.h>

#include "hdf5.h"
#include <assert.h>
#include <math.h>

#define POWFLOP 5 /* assume a pow() takes five flops */
#define DIFF_NORM 1.0e-10

typedef struct {
  PetscScalar var[NVAR];
} Field;

/* User-defined application context */
typedef struct {
  PetscReal tolgm;
  PetscReal rtol;
  PetscReal atol;
  PetscReal damp;
  PetscReal dt;
  PetscReal tmax;
  int       nxd;
  int       nyd;
  int       nzd;
  int       numtime;
  int       maxitnwt;
  int       maxksp;
  int       maxitgm;
  int       method;
  int       global;
  int       iguess;
  int       bcs[6];
} input_CTX;

typedef struct {
  DA	      da;
  input_CTX   indata;
  int         nwits;
  int         gmits;
  int         ierr;
  PetscReal   dt;
  PetscReal   time;
  Vec         x0;
  Vec         xold;
  Vec         fold;
  Vec         fsrc;
  PetscTruth  hdf5c;
} AppCtx;


#ifdef absoft
extern void FORTRAN_NAME(EVALUATENONLINEARRESIDUAL) 
                      (Field*, Field*, int*, int*, int*, int*, int*, int*,
		       int*, int*, int*, int*, int*, int*);
extern void FORTRAN_NAME(FORMEQUILIBRIUM) (Field*, int*, int*, int*, int*, int*, int*);
extern void FORTRAN_NAME(FORMINITIALCONDITION) (Field*,int*,int*,int*,int*,int*,int*,PetscScalar*);

extern void FORTRAN_NAME(WRITEOUTPUTDATA) (Field*,int*,int*,int*,int*,int*,int*,int*,int*);
extern void FORTRAN_NAME(CORRECTTIMESTEP) (PetscScalar*,PetscScalar*,PetscScalar*,int*,PetscScalar*);
extern void FORTRAN_NAME(FORTRANDESTROY) ();
extern void FORTRAN_NAME(READINPUTFILE) (input_CTX*);

#else
extern void FORTRAN_NAME(evaluatenonlinearresidual) 
                      (Field*, Field*, int*, int*, int*, int*, int*, int*,
		       int*, int*, int*, int*, int*, int*);
extern void FORTRAN_NAME(formequilibrium) (Field*, int*, int*, int*, int*, int*, int*);
extern void FORTRAN_NAME(forminitialcondition) (Field*,int*,int*,int*,int*,int*,int*,PetscScalar*);

extern void FORTRAN_NAME(writeoutputdata) (Field*,int*,int*,int*,int*,int*,int*,int*,int*);
extern void FORTRAN_NAME(correcttimestep) (PetscScalar*,PetscScalar*,PetscScalar*,int*,PetscScalar*);
extern void FORTRAN_NAME(fortrandestroy) ();
extern void FORTRAN_NAME(readinputfile) (input_CTX*);

#endif

/* User-defined routines */


extern int ReadInputNInitialize(input_CTX*, DAPeriodicType*);
extern int FormFunction        (SNES, Vec, Vec, void*);
extern int FormEquilibrium     (SNES, Vec, void*);
extern int FormInitialCondition(SNES, Vec, void*);
extern int correctTimeStep     (SNES, Vec, void*);
extern int Monitor             (SNES, int, double, void*);
extern int MatrixFreePreconditioner(void*,Vec,Vec);

#if !defined(lahey)
int main(int argc, char **argv)
#else
int MAIN__(int argc, char **argv)
#endif
{
  SNES     snes;	        /* SNES context */
  KSP      ksp;            /* KSP context  */
  PC       pc;             /* PC context   */
  Vec	   x, r;	        /* Vec x: solution, r: residual */
  AppCtx   user;
  
  int     ierr;
  int     nxd = 32;         /* x-direction grid size, i.e., nx */
  int     nyd = 32;         /* y-direction grid size, i.e., ny */
  int     nzd = 1;          /* z-direction grid size, i.e., nz */
  int     numtime;
  int     i, steps=0;
  int     np,npx,npy,npz;           /* num. of processors used */
  int     my_rank;        /* id of my processor */
  PetscReal atol;
  PetscReal rtol;
  PetscReal tolgm;
  PetscReal time,tmax;
  int       maxitnwt;
  int       maxitgm;
 
  DAPeriodicType     BC=DA_NONPERIODIC;

  Mat                J;     /* Jacobian matrix for FD approximation */
  PetscTruth         matrix_free;
  PetscTruth         user_precond;

  char *datasetname[] = {"Rho", "Vx", "Vy", "Vz", "Bx", "By", "Bz", "Temp"};
  
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Begin program
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  */	
  ierr = PetscInitialize(&argc, &argv, (char *)0, help);CHKERRQ(ierr);

  MPI_Comm_rank(PETSC_COMM_WORLD, &my_rank);
  MPI_Comm_size(PETSC_COMM_WORLD, &np);


  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Set problem parameters
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  */

  user.ierr = 0;

  ierr = ReadInputNInitialize(&user.indata, &BC);CHKERRQ(ierr);
  
  nxd = user.indata.nxd;
  nyd = user.indata.nyd;
  nzd = user.indata.nzd;
  
  atol 	   = PETSC_DEFAULT;
  rtol 	   = user.indata.rtol;
  tolgm    = user.indata.tolgm;
  numtime  = user.indata.numtime;
  if (numtime < 0) numtime = 10000000;
  tmax     = user.indata.tmax;
  if (tmax == 0.) tmax = 1e30;
  maxitnwt = user.indata.maxitnwt;
  if (maxitnwt == 0) maxitnwt = (int) PetscMax((1.5*log(rtol)/log(tolgm)),10.);

  user.dt = user.indata.dt;

  /* Set runtime options */
  
  ierr = PetscOptionsGetInt (PETSC_NULL,"-nmax",&numtime,PETSC_NULL);CHKERRQ(ierr);
  
  ierr = PetscOptionsGetReal(PETSC_NULL,"-tmax",&tmax   ,PETSC_NULL);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Set processor distribution
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  */

  npx = PETSC_DECIDE;
  npy = PETSC_DECIDE;
  npz = PETSC_DECIDE;
  if (np >= 2) {
    if (np%2 == 0) {
      /*      npx = 2;*/
    }
  }

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Create SNES context 
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  */

  ierr = SNESCreate(PETSC_COMM_WORLD, &snes);CHKERRQ(ierr);

  ierr = DACreate3d(PETSC_COMM_WORLD,BC,DA_STENCIL_BOX,nxd,nyd,nzd,npx,npy,npz,\
		    NVAR,1,PETSC_NULL,PETSC_NULL,PETSC_NULL,&user.da);CHKERRQ(ierr);
  
  ierr = DACreateGlobalVector(user.da,&x);CHKERRQ(ierr);
  ierr = VecDuplicate(x,&r);CHKERRQ(ierr);
  ierr = VecDuplicate(x,&user.x0  );CHKERRQ(ierr);
  ierr = VecDuplicate(x,&user.xold);CHKERRQ(ierr);
  ierr = VecDuplicate(x,&user.fold);CHKERRQ(ierr);
  ierr = VecDuplicate(x,&user.fsrc);CHKERRQ(ierr);
  
  ierr = SNESSetFunction(snes,r,FormFunction,(void*)&user);CHKERRQ(ierr);
  
  ierr = SNESSetMonitor(snes,Monitor,(void*)&user,PETSC_NULL);CHKERRQ(ierr);

  
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Customize nonlinear solver; set runtime options
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  
  /* 
     Set linear solver defaults for this problem. By extracting the
     SLES, KSP, and PC contexts from the SNES context, we can then
     directly call any SLES, KSP, and PC routines to set various options.
  */

  /* SNESSetTolerances(snes,atol,rtol,stol,maxit,maxf); */
  
  ierr = SNESSetTolerances(snes,atol,rtol,PETSC_DEFAULT,maxitnwt,PETSC_DEFAULT);CHKERRQ(ierr);

  /* KSPSetTolerances(ksp,rtol,atol,dtol,maxits); */
  
  /* Set preconditioner for matrix-free method */
  ierr = PetscOptionsHasName(PETSC_NULL,"-snes_mf",&matrix_free);CHKERRQ(ierr);
  if (matrix_free) {
    ierr = SNESGetKSP (snes, &ksp); CHKERRQ(ierr);
    ierr = KSPGetPC   (ksp, &pc); CHKERRQ(ierr);
    ierr = KSPSetTolerances(ksp,tolgm,PETSC_DEFAULT,PETSC_DEFAULT,maxitgm);CHKERRQ(ierr);
    ierr = PetscOptionsHasName(PETSC_NULL,"-user_precond",&user_precond);CHKERRQ(ierr);
    if (user_precond) {
      ierr = PCSetType(pc,PCSHELL);CHKERRQ(ierr);
      ierr = PCShellSetApply(pc,MatrixFreePreconditioner,PETSC_NULL);CHKERRQ(ierr);
      ierr = PCShellSetName(pc,"user-defined preconditioner");CHKERRQ(ierr);
    } else {
      ierr = PCSetType(pc,PCNONE); CHKERRQ(ierr);
    }
  } else {
    SETERRQ(1, "You should use \"matrix-free method\", i.e., -snes_mf");
  }

  /* Set SNES run time options */
  
  ierr = SNESSetFromOptions(snes);CHKERRQ(ierr);
  
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Initialize calculation 
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  ierr = FormEquilibrium(snes, x, &user);CHKERRQ(ierr);

  ierr = VecCopy(x, user.x0  );CHKERRQ(ierr);

  ierr = VecCopy(x, user.xold);CHKERRQ(ierr);

  ierr = EvaluateFunction(snes,x,user.fsrc,&user);CHKERRQ(ierr);

  ierr = FormInitialCondition(snes, x, &user);CHKERRQ(ierr);
  /*ierr = FormInitialCondition(snes, user.fsrc, &user);CHKERRQ(ierr);*/

  time = user.time+user.dt;

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Time stepping  
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  
  for (steps=1; (steps<=numtime)&&(time<=tmax); steps++,time+=user.dt) {
    
    ierr = ProcessOldSolution(snes,x,&user);CHKERRQ(ierr);

    ierr = VecCopy(x, user.xold);CHKERRQ(ierr);

    ierr = SNESSolve(snes, x); CHKERRQ(ierr);

    ierr = SNESGetIterationNumber(snes,&user.nwits);CHKERRQ(ierr);

    ierr = SNESGetNumberLinearIterations(snes,&user.gmits);CHKERRQ(ierr);

    /*
    PetscPrintf(PETSC_COMM_WORLD,"Time = %g, Tmax = %g \n",time,tmax);
    PetscPrintf(PETSC_COMM_WORLD,"Steps = %d Number of Newton iterations = %d\n"\
      ,steps, user.nwits);
    PetscPrintf(PETSC_COMM_WORLD,"Steps = %d Number of Krylov iterations = %d\n"\
      ,steps, user.gmits);
    */
  }
  
  ierr = ProcessOldSolution(snes,x,&user);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Deallocate memory and finish
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  
#ifdef absoft
  FORTRAN_NAME(FORTRANDESTROY) ();
#else
  FORTRAN_NAME(fortrandestroy) ();
#endif

  ierr = VecDestroy(x);CHKERRQ(ierr);
  ierr = VecDestroy(r);CHKERRQ(ierr);
  ierr = VecDestroy(user.x0  );CHKERRQ(ierr);
  ierr = VecDestroy(user.xold);CHKERRQ(ierr);
  ierr = VecDestroy(user.fold);CHKERRQ(ierr);
  ierr = VecDestroy(user.fsrc);CHKERRQ(ierr);
  ierr = SNESDestroy(snes);CHKERRQ(ierr);
  ierr = DADestroy(user.da);CHKERRQ(ierr);
  PetscFinalize();
  
  PetscFunctionReturn(0);
}

/* -------------------- Read input file routine ------------- */
#undef __FUNCT__
#define __FUNCT__ "ReadInputNInitialize"
int ReadInputNInitialize(input_CTX* data, DAPeriodicType* BC)
{

  char       dummy[256], dummy_ch, *bc_def = "'def'";
  int        dummy_int, i, ierr;
  PetscReal  dummy_float;
  FILE       *fp;

  PetscFunctionBegin;

#ifdef absoft
  FORTRAN_NAME(READINPUTFILE) (data);
#else
  FORTRAN_NAME(readinputfile) (data);
#endif

  if (data->bcs[0] && data->bcs[2] && data->bcs[4]) {
    *BC=DA_XYZPERIODIC;
  } else if (data->bcs[0] && data->bcs[2] && !data->bcs[4]) {
    *BC=DA_XYPERIODIC;
  } else if (data->bcs[0] && !data->bcs[2] && data->bcs[4]) {
    *BC=DA_XZPERIODIC;
  } else if (data->bcs[0] && !data->bcs[2] && !data->bcs[4]) {
    *BC=DA_XPERIODIC;
  } else if (!data->bcs[0] && data->bcs[2] && data->bcs[4]) {
    *BC=DA_YZPERIODIC;
  } else if (!data->bcs[0] && data->bcs[2] && !data->bcs[4]) {
    *BC=DA_YPERIODIC;
  } else if (!data->bcs[0] && !data->bcs[2] && data->bcs[4]) {
    *BC=DA_ZPERIODIC;
  } else {
    *BC=DA_NONPERIODIC;
  }

  ierr = PetscOptionsGetInt(PETSC_NULL,"-nx",&data->nxd,PETSC_NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetInt(PETSC_NULL,"-ny",&data->nyd,PETSC_NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetInt(PETSC_NULL,"-nz",&data->nzd,PETSC_NULL);CHKERRQ(ierr);
  
  PetscFunctionReturn(0);
}


/* -------------------- User-defined Monitor() routine ------------- */
#undef __FUNCT__
#define __FUNCT__ "Monitor"
int Monitor(SNES snes, int its, double fnorm, void *ctx)
{

  PetscFunctionBegin;
  
  /* EMPTY */
  
  PetscFunctionReturn(0);
}

/* --------------------  Form equilibrium ------------------------- */
#undef __FUNCT__
#define __FUNCT__ "FormEquilibrium"
int FormEquilibrium(SNES snes,Vec X,void *ptr)
{
  AppCtx  *user = (AppCtx*)ptr;
  Field	  ***xvec;
  int	  ierr,i,j,k,xs,ys,zs,xm,ym,zm,ze,ye,xe,mx,my,mz;
  Vec     localX;

  PetscFunctionBegin;

  /*
   * Initialize counters
   */

  user->nwits = 0;
  user->gmits = 0;
    
  /*
   * Get pointers to vector data
   */
  ierr = DAVecGetArray(user->da, X, (void**)&xvec);CHKERRQ(ierr);
  
  /*
   * Get local grid boundaries
   */
  ierr = DAGetCorners(user->da, &xs, &ys, &zs, &xm, &ym, &zm);CHKERRQ(ierr);
  
  zs = zs + 1;
  ys = ys + 1;
  xs = xs + 1;
  
  ze = zs + zm - 1; 
  ye = ys + ym - 1; 
  xe = xs + xm - 1;

  /*
   * Compute function over the locally owned part of the grid
   */

#ifdef absoft
  FORTRAN_NAME(FORMEQUILIBRIUM)(&(xvec[zs-1][ys-1][xs-1])\
                                           ,&xs,&xe,&ys,&ye,&zs,&ze);
#else
  FORTRAN_NAME(formequilibrium)(&(xvec[zs-1][ys-1][xs-1])\
                                           ,&xs,&xe,&ys,&ye,&zs,&ze);
#endif


  /* Restore vectors */
  ierr = DAVecRestoreArray(user->da, X, (void**)&xvec);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

/* --------------------  Form initial condition ----------------- */
#undef __FUNCT__
#define __FUNCT__ "FormInitialCondition"
int FormInitialCondition(SNES snes,Vec X,void *ptr)
{
  AppCtx  *user = (AppCtx*)ptr;
  Field	  ***xvec;
  int	  ierr,i,j,k,xs,ys,zs,xm,ym,zm,mx,my,mz;
  int     ze,ye,xe;
  int     xs_g,ys_g,zs_g,ze_g,ye_g,xe_g,xm_g,ym_g,zm_g;
  Vec     localX;

  PetscFunctionBegin;

  /*
   * Initialize counters
   */

  user->nwits = 0;
  user->gmits = 0;

  /*
   * Scatter ghost points to local vector,using the 2-step process
   * DAGlobalToLocalBegin(),DAGlobalToLocalEnd().
   * By placing code between these two statements, computations can be
   * done while messages are in transition.
   */
  
  ierr = DAGetLocalVector(user->da,&localX);CHKERRQ(ierr);
  
  ierr = DAGlobalToLocalBegin(user->da,X,INSERT_VALUES,localX);CHKERRQ(ierr);
  ierr = DAGlobalToLocalEnd  (user->da,X,INSERT_VALUES,localX);CHKERRQ(ierr);
  
  /*
   * Get pointers to vector data
   */
  ierr = DAVecGetArray(user->da, localX, (void**)&xvec);CHKERRQ(ierr);
  
  /*
   * Get local grid boundaries
   */
  ierr = DAGetCorners(user->da, &xs, &ys, &zs, &xm, &ym, &zm);CHKERRQ(ierr);
  
  zs = zs + 1;
  ys = ys + 1;
  xs = xs + 1;
  
  ze = zs + zm - 1; 
  ye = ys + ym - 1; 
  xe = xs + xm - 1;
  
  /*
   * Obtain ghost cell boundaries
   */
  ierr = DAGetGhostCorners(user->da,&xs_g,&ys_g,&zs_g,&xm_g,&ym_g,&zm_g);CHKERRQ(ierr);

  zs_g = zs_g + 1;
  ys_g = ys_g + 1;
  xs_g = xs_g + 1;

  ze_g = zs_g + zm_g - 1;
  ye_g = ys_g + ym_g - 1;
  xe_g = xs_g + xm_g - 1;

  /*
   * Compute initial condition over the locally owned part of the grid, including 
   * ghost cells (which are needed for restarting the calculation)
   */

#ifdef absoft
  FORTRAN_NAME(FORMINITIALCONDITION)(&(xvec[zs_g-1][ys_g-1][xs_g-1])\
                                           ,&xs_g,&xe_g,&ys_g,&ye_g,&zs_g,&ze_g,&user->time);
#else
  FORTRAN_NAME(forminitialcondition)(&(xvec[zs_g-1][ys_g-1][xs_g-1])\
                                           ,&xs_g,&xe_g,&ys_g,&ye_g,&zs_g,&ze_g,&user->time);
#endif

  /* Restore vectors */
  ierr = DAVecRestoreArray(user->da, localX, (void**)&xvec);CHKERRQ(ierr);
  ierr = DARestoreLocalVector(user->da,&localX);CHKERRQ(ierr);
  ierr = DALocalToGlobal(user->da,localX,INSERT_VALUES,X);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}


/* --------------------  Process old time solution ----------------- */
#undef __FUNCT__
#define __FUNCT__ "ProcessOldSolution"
int ProcessOldSolution(SNES snes,Vec X,void *ptr)
{
  AppCtx        *user = (AppCtx*)ptr;
  Field         ***xvec;
  PetscScalar   hx, hy, hz, xp, yp;
  int           ierr,i,j,k,xs,ys,zs,xm,ym,zm,mx,my,mz;
  int           ze,ye,xe;
  int           xs_g,ys_g,zs_g,ze_g,ye_g,xe_g,xm_g,ym_g,zm_g;
  Vec           localX;

  PetscFunctionBegin;

  /*
   * Form old time fluxes
   */

  ierr = EvaluateFunction(snes, X, user->fold, (void*)user);CHKERRQ(ierr);

  /*
   * Scatter ghost points to local vector,using the 2-step process
   * DAGlobalToLocalBegin(),DAGlobalToLocalEnd().
   * By placing code between these two statements, computations can be
   * done while messages are in transition.
   */
  
  ierr = DAGetLocalVector(user->da,&localX);CHKERRQ(ierr);
  
  ierr = DAGlobalToLocalBegin(user->da,X,INSERT_VALUES,localX);CHKERRQ(ierr);
  ierr = DAGlobalToLocalEnd  (user->da,X,INSERT_VALUES,localX);CHKERRQ(ierr);
  
  /*
   * Get pointers to vector data
   */
  ierr = DAVecGetArray(user->da, localX, (void**)&xvec);CHKERRQ(ierr);
  
  /*
   * Get local grid boundaries
   */
  ierr = DAGetCorners(user->da, &xs, &ys, &zs, &xm, &ym, &zm);CHKERRQ(ierr);
  
  zs = zs + 1;
  ys = ys + 1;
  xs = xs + 1;
  
  ze = zs + zm - 1; 
  ye = ys + ym - 1; 
  xe = xs + xm - 1;
  
  /*
   * Obtain ghost cell boundaries
   */
  ierr = DAGetGhostCorners(user->da,&xs_g,&ys_g,&zs_g,&xm_g,&ym_g,&zm_g);CHKERRQ(ierr);

  zs_g = zs_g + 1;
  ys_g = ys_g + 1;
  xs_g = xs_g + 1;

  ze_g = zs_g + zm_g - 1;
  ye_g = ys_g + ym_g - 1;
  xe_g = xs_g + xm_g - 1;

  /*
   * Write output data
   */

#ifdef absoft
  FORTRAN_NAME(WRITEOUTPUTDATA)(&(xvec[zs_g-1][ys_g-1][xs_g-1])\
				   ,&xs_g,&xe_g,&ys_g,&ye_g,&zs_g,&ze_g\
				   ,&user->gmits,&user->nwits);
#else
  FORTRAN_NAME(writeoutputdata)(&(xvec[zs_g-1][ys_g-1][xs_g-1])\
				   ,&xs_g,&xe_g,&ys_g,&ye_g,&zs_g,&ze_g\
				   ,&user->gmits,&user->nwits);
#endif

  /*
   * Correct time step
   */

  ierr = correctTimeStep(snes,X,user);CHKERRQ(ierr);

  /* 
   * Restore vectors 
   */
  ierr = DAVecRestoreArray(user->da, localX, (void**)&xvec);CHKERRQ(ierr);
  ierr = DARestoreLocalVector(user->da,&localX);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

/* --------------------  Correct time step ----------------- */
#undef __FUNCT__
#define __FUNCT__ "correctTimeStep"
int correctTimeStep(SNES snes,Vec X,void *ptr)
{
  AppCtx        *user = (AppCtx*)ptr;
  int           ierr,ieq;
  PetscScalar   dnp[NVAR],dn[NVAR],dnh[NVAR];
  Vec           dummy;
  PetscScalar   one=1.0, mone=-1.0, half=0.5;

  PetscFunctionBegin;

  ierr = VecDuplicate(X,&dummy); CHKERRQ(ierr);

  /*
   * Evaluate dnp (norm of perturbation from equilibrium of time level n+1)
   */

  ierr = VecCopy(X,dummy); CHKERRQ(ierr);

  ierr = VecAXPY(&mone,user->x0,dummy);CHKERRQ(ierr);

  for (ieq=0;ieq < NVAR;ieq++) {
    ierr = VecStrideNorm(dummy,ieq,NORM_2,&dnp[ieq]);CHKERRQ(ierr);
  }

  /*
   * Evaluate dn (norm of perturbation from equilibrium of time level n)
   */

  ierr = VecCopy(user->xold,dummy); CHKERRQ(ierr);

  ierr = VecAXPY(&mone,user->x0,dummy);CHKERRQ(ierr);

  for (ieq=0;ieq < NVAR;ieq++) {
    ierr = VecStrideNorm(dummy,ieq,NORM_2,&dn[ieq]);CHKERRQ(ierr);
  }

  /*
   * Evaluate dnh (norm of perturbation from equilibrium of time level n+1/2)
   */

  ierr = VecCopy(X,dummy); CHKERRQ(ierr);

  ierr = VecAXPY (&one       ,user->xold,dummy);CHKERRQ(ierr);
  ierr = VecAXPBY(&mone,&half,user->x0  ,dummy);CHKERRQ(ierr);

  for (ieq=0;ieq < NVAR;ieq++) {
    ierr = VecStrideNorm(dummy,ieq,NORM_2,&dnh[ieq]);CHKERRQ(ierr);
  }

  /*
   * Correct time step (using previous norm quantities)
   */

#ifdef absoft
  FORTRAN_NAME(CORRECTTIMESTEP)(&(dn[0]),&(dnh[0]),&(dnp[0]),&user->ierr,&user->dt);
#else
  FORTRAN_NAME(correcttimestep)(&(dn[0]),&(dnh[0]),&(dnp[0]),&user->ierr,&user->dt);
#endif

  /*
   * Deallocate memory
   */

  ierr = VecDestroy(dummy);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

/* --------------------  Evaluate Function G(x,t) --------------------- */
#undef __FUNCT__
#define __FUNCT__ "FormFunction"
int FormFunction(SNES snes,Vec X,Vec F,void* ptr)
{
  AppCtx  *user = (AppCtx*)ptr;
  
  int     ierr,i,j,k,l;
  int     xs,ys,zs,xm,ym,zm;
  Field   ***x,***f,***xold,***fold,***fsrc;
  Vec     localX,localXold;
  double  dt, theta;
  
  int     ye,xe,ze;
  int     ivar,my_rank, xs_g,ys_g,zs_g,ze_g,ye_g,xe_g,xm_g,ym_g,zm_g;

  double cnf[NVAR],one_over_dt[NVAR];

  PetscFunctionBegin;

  ierr = EvaluateFunction(snes, X, F, (void*)user);CHKERRQ(ierr);

  /*
   * Scatter ghost points to local vector,using the 2-step process
   * DAGlobalToLocalBegin(),DAGlobalToLocalEnd().
   * By placing code between these two statements, computations can be
   * done while messages are in transition.
   */
  
  ierr = DAGetLocalVector(user->da,&localX);CHKERRQ(ierr);
  ierr = VecDuplicate(localX,&localXold);CHKERRQ(ierr)

  ierr = DAGlobalToLocalBegin(user->da,X,INSERT_VALUES,localX);CHKERRQ(ierr);
  ierr = DAGlobalToLocalEnd  (user->da,X,INSERT_VALUES,localX);CHKERRQ(ierr);
  ierr = DAGlobalToLocalBegin(user->da,user->xold,INSERT_VALUES,localXold);CHKERRQ(ierr);
  ierr = DAGlobalToLocalEnd  (user->da,user->xold,INSERT_VALUES,localXold);CHKERRQ(ierr);

  /* 
   *Get pointers to vector data
   */

  ierr = DAVecGetArray(user->da,localX,(void**)&x);CHKERRQ(ierr);
  ierr = DAVecGetArray(user->da,F     ,(void**)&f);CHKERRQ(ierr);
  ierr = DAVecGetArray(user->da,localXold,(void**)&xold);CHKERRQ(ierr);
  ierr = DAVecGetArray(user->da,user->fold,(void**)&fold);CHKERRQ(ierr);
  ierr = DAVecGetArray(user->da,user->fsrc,(void**)&fsrc);CHKERRQ(ierr);

  /*
   * Get local grid boundaries
   */
  ierr = DAGetCorners(user->da, &xs, &ys, &zs, &xm, &ym, &zm);CHKERRQ(ierr);
	
  /*
   * Evaluate Time Step...
   */

#ifdef absoft
  FORTRAN_NAME(DEFINETSPARAMETERSROUTINE) (&(cnf[0]),&(one_over_dt[0]));
#else
  FORTRAN_NAME(definetsparametersroutine) (&(cnf[0]),&(one_over_dt[0]));
#endif

  for (k=zs; k<zs+zm; k++) {
    for (j=ys; j<ys+ym; j++) {
      for (i=xs; i<xs+xm; i++) {
	for (l=0; l<NVAR; l++) {
	  f[k][j][i].var[l] = (x[k][j][i].var[l] - xold[k][j][i].var[l])*one_over_dt[l]
	                     + (1.0 - cnf[l])*f   [k][j][i].var[l] 
	                     +        cnf[l] *fold[k][j][i].var[l]
                             -                fsrc[k][j][i].var[l];
	}
      }
    }
  }

  /*
   * Restore vectors
   */
  ierr = DAVecRestoreArray(user->da,localX,(void**)&x);CHKERRQ(ierr);
  ierr = DAVecRestoreArray(user->da,F     ,(void**)&f);CHKERRQ(ierr);

  ierr = DAVecRestoreArray(user->da,localXold,(void**)&xold);CHKERRQ(ierr);
  ierr = DAVecRestoreArray(user->da,user->fold,(void**)&fold);CHKERRQ(ierr);
  ierr = DAVecRestoreArray(user->da,user->fsrc,(void**)&fsrc);CHKERRQ(ierr);


  ierr = DARestoreLocalVector(user->da,&localX);CHKERRQ(ierr);
  ierr = DARestoreLocalVector(user->da,&localXold);CHKERRQ(ierr);

  ierr = PetscLogFlops((22 + 4*POWFLOP)*zm*ym*xm);CHKERRQ(ierr);
  PetscFunctionReturn(0);

}


/* --------------------  Evaluate Function F(x,t) --------------------- */
#undef __FUNCT__
#define __FUNCT__ "EvaluateFunction"
int EvaluateFunction(SNES snes,Vec X,Vec F,void* ptr)
{
  AppCtx  *user = (AppCtx*)ptr;
  
  int     ierr,i,j,k,l;
  int     xs,ys,zs,xm,ym,zm,mx,my,mz;
  Field   ***x,***f;
  Vec     localX;
  double  dt, theta;
  
  int     ye,xe,ze;
  int     xs_g,ys_g,zs_g,ze_g,ye_g,xe_g,xm_g,ym_g,zm_g;


  PetscFunctionBegin;

  /*
   * Scatter ghost points to local vector,using the 2-step process
   * DAGlobalToLocalBegin(),DAGlobalToLocalEnd().
   * By placing code between these two statements, computations can be
   * done while messages are in transition.
   */
  
  ierr = DAGetLocalVector(user->da,&localX);CHKERRQ(ierr);
  /*ierr = DAGetInfo(user->da,PETSC_NULL,&mx,&my,&mz,0,0,0,0,0,0,0);CHKERRQ(ierr);*/
  
  ierr = DAGlobalToLocalBegin(user->da,X,INSERT_VALUES,localX);CHKERRQ(ierr);
  ierr = DAGlobalToLocalEnd  (user->da,X,INSERT_VALUES,localX);CHKERRQ(ierr);
  
  /* 
   *Get pointers to vector data
   */

  ierr = DAVecGetArray(user->da,localX,(void**)&x);CHKERRQ(ierr);
  ierr = DAVecGetArray(user->da,F     ,(void**)&f);CHKERRQ(ierr);
  
  /*
   * Get local grid boundaries
   */
  ierr = DAGetCorners(user->da, &xs, &ys, &zs, &xm, &ym, &zm);CHKERRQ(ierr);
  
  zs = zs + 1;
  ys = ys + 1;
  xs = xs + 1;
  
  ze = zs + zm - 1; 
  ye = ys + ym - 1; 
  xe = xs + xm - 1;

  /*
   * Obtain ghost cell boundaries
   */
  ierr = DAGetGhostCorners(user->da,&xs_g,&ys_g,&zs_g,&xm_g,&ym_g,&zm_g);CHKERRQ(ierr);

  zs_g = zs_g + 1;
  ys_g = ys_g + 1;
  xs_g = xs_g + 1;

  ze_g = zs_g + zm_g - 1;
  ye_g = ys_g + ym_g - 1;
  xe_g = xs_g + xm_g - 1;
	
  /*
   * Evaluate nonlinear residual
   */

#ifdef absoft
  FORTRAN_NAME(EVALUATENONLINEARRESIDUAL) (&(x[zs_g-1][ys_g-1][xs_g-1]),&(f[zs-1][ys-1][xs-1])\
	  ,&xs,&xe,&ys,&ye,&zs,&ze,&xs_g,&xe_g,&ys_g,&ye_g,&zs_g,&ze_g);
#else
  FORTRAN_NAME(evaluatenonlinearresidual) (&(x[zs_g-1][ys_g-1][xs_g-1]),&(f[zs-1][ys-1][xs-1])\
	  ,&xs,&xe,&ys,&ye,&zs,&ze,&xs_g,&xe_g,&ys_g,&ye_g,&zs_g,&ze_g);
#endif

  /*
   * Restore vectors
   */
  ierr = DAVecRestoreArray(user->da,localX,(void**)&x);CHKERRQ(ierr);
  ierr = DAVecRestoreArray(user->da,F     ,(void**)&f);CHKERRQ(ierr);
  ierr = DARestoreLocalVector(user->da,&localX);CHKERRQ(ierr);
  ierr = PetscLogFlops((22 + 4*POWFLOP)*zm*ym*xm);CHKERRQ(ierr);
  PetscFunctionReturn(0);

}

#undef __FUNCT__
#define __FUNCT__ "MatrixFreePreconditioer"
/* for future usage... See ex6.c of SNES */
/* ------------------------------------------------------------------- */
/*
   Input Parameters:
   ctx - optional user-defined context, as set by PCShellSetApply()
   x - input vector

   Output Parameter:
   y - preconditioned vector
*/
int MatrixFreePreconditioner(void *ctx,Vec x,Vec y)
{

  PetscFunctionBegin;
  
  VecCopy(x,y);
  
  PetscFunctionReturn(0);
}
