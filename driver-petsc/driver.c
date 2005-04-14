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
  int       ilevel;
  int       nxd;
  int       nyd;
  int       nzd;
  int       npx;
  int       npy;
  int       npz;
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
  int         ierr;
  PetscReal   time;
  PetscReal   ksp_res0;
  int         ksp_its;
  PetscReal   snes_res0;
  int         snes_its;
  int         aspc_its;
  Vec         x0;
  Vec         xold;
  Vec         xk;
  Vec         dzk;
  Vec         rk;
  Mat         J;
  PetscTruth  hdf5c;
  PetscTruth  petsc_PC;
} AppCtx;


#ifdef absoft
extern void FORTRAN_NAME(EVALUATENONLINEARRESIDUAL) 
                      (Field*, Field*, int*, int*, int*, int*, int*, int*,
		       int*, int*, int*, int*, int*, int*);
extern void FORTRAN_NAME(FORMEQUILIBRIUM) (Field*, int*, int*, int*, int*, int*, int*);
extern void FORTRAN_NAME(FORMINITIALCONDITION) (Field*,int*,int*,int*,int*,int*,int*,PetscScalar*);
extern void FORTRAN_NAME(PROCESSOLDSOLUTION) (Field*,int*,int*,int*,int*,int*,int*,int*,int*);
extern void FORTRAN_NAME(CORRECTTIMESTEP) (PetscScalar*,PetscScalar*,PetscScalar*,int*,PetscScalar*);
extern void FORTRAN_NAME(FORTRANDESTROY) ();
extern void FORTRAN_NAME(READINPUTFILE) (input_CTX*);
extern void FORTRAN_NAME(SETUPSHELLPC) (Field*, int*, int*, int*, int*, int*, int*);
extern void FORTRAN_NAME(APPLYSHELLPC) (Field*, Field*, int*, int*, int*, int*, int*, int*);
#else
extern void FORTRAN_NAME(evaluatenonlinearresidual) 
                      (Field*, Field*, int*, int*, int*, int*, int*, int*,
		       int*, int*, int*, int*, int*, int*);
extern void FORTRAN_NAME(formequilibrium) (Field*, int*, int*, int*, int*, int*, int*);
extern void FORTRAN_NAME(forminitialcondition) (Field*,int*,int*,int*,int*,int*,int*,PetscScalar*);
extern void FORTRAN_NAME(processoldsolution) (Field*,int*,int*,int*,int*,int*,int*,int*,int*);
extern void FORTRAN_NAME(correcttimestep) (PetscScalar*,PetscScalar*,PetscScalar*,int*,PetscScalar*);
extern void FORTRAN_NAME(fortrandestroy) ();
extern void FORTRAN_NAME(readinputfile) (input_CTX*);
extern void FORTRAN_NAME(setupshellpc) (Field*, int*, int*, int*, int*, int*, int*);
extern void FORTRAN_NAME(applyshellpc) (Field*, Field*, int*, int*, int*, int*, int*, int*);
#endif

/* User-defined routines */

extern int ReadInputNInitialize(input_CTX*, DAPeriodicType*);
extern int FormFunction        (SNES, Vec, Vec, void*);
extern int EvaluateFunction    (SNES, Vec, Vec, void*);
extern int FormEquilibrium     (SNES, Vec, void*);
extern int FormInitialCondition(SNES, Vec, void*);
extern int correctTimeStep     (SNES, Vec, void*);
extern int MySNESMonitor       (SNES, int, PetscReal, void*);
extern int MyKSPMonitor        (KSP , int, PetscReal, void*);
extern int ApplyPC             (void*,Vec,Vec);
extern int ApplyASPC           (void*,Vec,Vec);
extern int SetupPreconditioner (void*);

#if !defined(lahey)
int main(int argc, char **argv)
#else
int MAIN__(int argc, char **argv)
#endif
{
  SNES      snes;	        /* SNES context */
  KSP       ksp;                /* KSP context  */
  PC        pc;                 /* PC context   */
  Vec	    x, r;	        /* Vec x: solution, r: residual */
  AppCtx    user;
  
  int       ierr;
  int       nxd = 32;         	/* x-direction grid size, i.e., nx */
  int       nyd = 32;         	/* y-direction grid size, i.e., ny */
  int       nzd = 1;          	/* z-direction grid size, i.e., nz */
  int       numtime;
  int       i, steps=0;
  int       np,npx,npy,npz;     /* num. of processors used */
  int       my_rank;            /* id of my processor */
  PetscReal atol;
  PetscReal rtol;
  PetscReal tolgm;
  PetscReal time,tmax;
  int       maxitnwt;
  int       maxitgm;
 
  PetscScalar   zero=0.0;

  DAPeriodicType     BC=DA_NONPERIODIC;

  PetscTruth         matrix_free,user_precond;

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Begin program
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  */	
  ierr = PetscInitialize(&argc, &argv, (char *)0, help);CHKERRQ(ierr);
  ierr = PetscInitializeFortran();CHKERRQ(ierr);

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

  npx = user.indata.npx;
  npy = user.indata.npy;
  npz = user.indata.npz;

  if (npx == 0) npx = PETSC_DECIDE;
  if (npy == 0) npy = PETSC_DECIDE;
  if (npz == 0) npz = PETSC_DECIDE;

  atol 	   = PETSC_DEFAULT;
  rtol 	   = user.indata.rtol;
  tolgm    = user.indata.tolgm;

  numtime  = user.indata.numtime;
  if (numtime < 0) numtime = 10000000;

  tmax     = user.indata.tmax;
  if (tmax == 0.) tmax = 1e30;

  maxitgm  = user.indata.maxksp;

  maxitnwt = user.indata.maxitnwt;
  if (maxitnwt == 0) maxitnwt = (int) PetscMax((1.5*log(rtol)/log(tolgm)),10.);

  user.aspc_its = 0;   /* Initialize Additive Schwartz method */

  user.petsc_PC = 0;   /* Do not use PETSC PC */

  /* Set runtime options */
  
  ierr = PetscOptionsGetInt(PETSC_NULL,"-nmax",&numtime,PETSC_NULL)	        ;CHKERRQ(ierr);

  ierr = PetscOptionsGetInt(PETSC_NULL,"-ilevel",&user.indata.ilevel,PETSC_NULL);CHKERRQ(ierr);
  
  ierr = PetscOptionsGetReal(PETSC_NULL,"-tmax",&tmax   ,PETSC_NULL)	        ;CHKERRQ(ierr);

  ierr = PetscOptionsHasName(PETSC_NULL,"-snes_mf",&matrix_free)                ;CHKERRQ(ierr);

  ierr = PetscOptionsGetInt(PETSC_NULL,"-aspc_its",&user.aspc_its,PETSC_NULL)   ;CHKERRQ(ierr);

  ierr = PetscOptionsGetLogical(PETSC_NULL,"-petsc_PC",&user.petsc_PC,PETSC_NULL);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Create DA context 
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  */

  ierr = SNESCreate(PETSC_COMM_WORLD, &snes)                                    ;CHKERRQ(ierr);

  /* Create DA */
  ierr = DACreate3d(PETSC_COMM_WORLD,BC,DA_STENCIL_BOX,nxd,nyd,nzd,npx,npy,npz,\
		    NVAR,1,PETSC_NULL,PETSC_NULL,PETSC_NULL,&user.da)	       	;CHKERRQ(ierr);
  
  ierr = DACreateGlobalVector(user.da,&x)		      	     	       	;CHKERRQ(ierr);

  /* Initialize vectors */
  ierr = VecDuplicate(x,&r)        	 		      	     	       	;CHKERRQ(ierr);
  ierr = VecDuplicate(x,&user.x0  )	 		      	     	       	;CHKERRQ(ierr);
  ierr = VecDuplicate(x,&user.xold)	 		      	     	       	;CHKERRQ(ierr);
  ierr = VecDuplicate(x,&user.dzk )	 		      	     	       	;CHKERRQ(ierr);
  ierr = VecDuplicate(x,&user.rk  )	 		      	     	       	;CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Create SNES context 
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  */

  /* Initialize SNES */
  ierr = SNESSetFunction(snes,r,EvaluateFunction,(void*)&user)	     	       	;CHKERRQ(ierr);

  /* Set SNES monitor */
  ierr = SNESSetMonitor(snes,MySNESMonitor,(void*)&user,PETSC_NULL)  	       	;CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Customize nonlinear and linear solvers; set runtime options
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  
  /* Customize SNES */
  
  ierr = SNESSetTolerances(snes,atol,rtol,PETSC_DEFAULT,maxitnwt,PETSC_DEFAULT)	;CHKERRQ(ierr);
  
  if (matrix_free) {

    /* Customize KSP */

    ierr = SNESGetKSP (snes, &ksp)                                             	;CHKERRQ(ierr);
    ierr = KSPSetType (ksp,KSPFGMRES)                                          	;CHKERRQ(ierr);
    ierr = KSPSetTolerances(ksp,tolgm,PETSC_DEFAULT,PETSC_DEFAULT,maxitgm)     	;CHKERRQ(ierr);
    ierr = KSPSetPreconditionerSide(ksp,PC_RIGHT);                             	;CHKERRQ(ierr);
    ierr = KSPSetMonitor(ksp,MyKSPMonitor,(void*)&user,PETSC_NULL)             	;CHKERRQ(ierr);
    if (user.indata.iguess == 1) ierr = KSPSetInitialGuessKnoll(ksp,PETSC_TRUE)	;CHKERRQ(ierr);

    /* Customize PC */

    ierr = KSPGetPC(ksp,&pc)                                                   	;CHKERRQ(ierr);
    if (!user.petsc_PC) {
      ierr = PCSetType(pc,PCSHELL)                                             	;CHKERRQ(ierr);
      if (user.aspc_its == 0) { 
	ierr = PCShellSetApply(pc,ApplyPC  ,(void*)&user)                      	;CHKERRQ(ierr);
      } else { 
	ierr = PCShellSetApply(pc,ApplyASPC,(void*)&user)                      	;CHKERRQ(ierr);
      }
      ierr = PCShellSetName(pc,"Physics-based preconditioner")                 	;CHKERRQ(ierr);
      ierr = PCShellSetSetUp(pc,SetupPreconditioner)                           	;CHKERRQ(ierr);
    } else {
      ierr = PCSetType(pc,PCNONE)                                              	;CHKERRQ(ierr);
    }

  } else {

    SETERRQ(1, "You should use \"matrix-free method\", i.e., -snes_mf");

  }

  /* Set SNES runtime options */
  
  ierr = SNESSetFromOptions(snes)                                              	;CHKERRQ(ierr);
  
  ierr = SNESSetUp(snes,x)                                                      ;CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Initialize calculation 
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  ierr = FormEquilibrium(snes, x, &user)     				       	;CHKERRQ(ierr);

  ierr = VecCopy(x, user.x0  )               				       	;CHKERRQ(ierr);

  ierr = VecCopy(x, user.xold)               				       	;CHKERRQ(ierr);

  ierr = FormInitialCondition(snes, x, &user)				       	;CHKERRQ(ierr);

  time = user.time+user.indata.dt;

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Time stepping  
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  
  for (steps=1; (steps<=numtime)&&(time<1.00001*tmax); steps++,time+=user.indata.dt) {

    /* Process nonlinear solution */
    ierr = ProcessOldSolution(snes,x,&user)           	    			;CHKERRQ(ierr);

    /* Save current time step */
    ierr = VecCopy(x, user.xold)                      	    			;CHKERRQ(ierr);

    /* Get initial Newton state info for PC */
    ierr = SNESGetSolution(snes,&user.xk)                                       ;CHKERRQ(ierr);
    ierr = SNESGetJacobian(snes,&user.J,PETSC_NULL,PETSC_NULL,PETSC_NULL)       ;CHKERRQ(ierr);
    
    /* Solve nonlinear problem */
    ierr = SNESSolve(snes, x)                         	    			;CHKERRQ(ierr);

    ierr = SNESGetIterationNumber(snes,&user.snes_its)	    			;CHKERRQ(ierr);

  }
  
  ierr = ProcessOldSolution(snes,x,&user)                                       ;CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Deallocate memory and finish
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  
#ifdef absoft
  FORTRAN_NAME(FORTRANDESTROY) ();
#else
  FORTRAN_NAME(fortrandestroy) ();
#endif

  ierr = VecDestroy(x)	      						       	;CHKERRQ(ierr);
  ierr = VecDestroy(r)	      						       	;CHKERRQ(ierr);
  ierr = VecDestroy(user.x0  )						       	;CHKERRQ(ierr);
  ierr = VecDestroy(user.xold)						       	;CHKERRQ(ierr);
  ierr = VecDestroy(user.dzk )						       	;CHKERRQ(ierr);
  ierr = VecDestroy(user.rk  )						       	;CHKERRQ(ierr);
  ierr = SNESDestroy(snes)    						       	;CHKERRQ(ierr);
  ierr = DADestroy(user.da)   						       	;CHKERRQ(ierr);

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
  
  PetscFunctionReturn(0);
}


/* -------------------- User-defined Monitor() routine ------------- */
#undef __FUNCT__
#define __FUNCT__ "MySNESMonitor"
int MySNESMonitor(SNES snes, int its, double fnorm, void *ctx)
{
  AppCtx  *user = (AppCtx*)ctx;

  int     ierr;
  PetscReal  rel_res;

  PetscFunctionBegin;

  /* Get nonlinear solver info relevant for PC */
  ierr = SNESGetSolution(snes,&user->xk)                                ; CHKERRQ(ierr);
  ierr = SNESGetJacobian(snes,&user->J,PETSC_NULL,PETSC_NULL,PETSC_NULL); CHKERRQ(ierr);

  /* Monitor */
  ierr = SNESGetNumberLinearIterations(snes,&user->ksp_its); CHKERRQ(ierr);
  if (user->indata.ilevel > 0) {
    if (its == 0) user->snes_res0 = fnorm;
    if (user->indata.ilevel > 1) PetscPrintf(PETSC_COMM_WORLD,"\n");
    rel_res = fnorm/user->snes_res0;
    PetscPrintf(PETSC_COMM_WORLD\
               ,"SNES its = %d ; res norm = %4.2e; res/rold norm = %4.2e; KSP its = %d \n"\
               ,its,fnorm,rel_res,user->ksp_its);
    if (user->indata.ilevel > 1) PetscPrintf(PETSC_COMM_WORLD,"\n");
  }

  PetscFunctionReturn(0);
}

/* -------------------- User-defined Monitor() routine ------------- */
#undef __FUNCT__
#define __FUNCT__ "MyKSPMonitor"
int MyKSPMonitor(KSP ksp, int its, double fnorm, void *ctx)
{
  AppCtx  *user = (AppCtx*)ctx;

  PetscReal  rel_res;

  PetscFunctionBegin;

  if (user->indata.ilevel > 1) {
    if (its == 0) user->ksp_res0 = fnorm;
    rel_res = fnorm/user->ksp_res0;
    PetscPrintf(PETSC_COMM_WORLD," KSP its = %d ; res/rold norm = %4.2e \n",its,rel_res);
  }
  
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

  PetscFunctionBegin;

  /*
   * Initialize counters
   */

  user->snes_its = 0;
  user->ksp_its = 0;
    
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

  /*ierr = EvaluateFunction(snes, X, user->fold, (void*)user);CHKERRQ(ierr);*/

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
   * Postprocess data
   */

#ifdef absoft
  FORTRAN_NAME(PROCESSOLDSOLUTION)(&(xvec[zs_g-1][ys_g-1][xs_g-1])\
				   ,&xs_g,&xe_g,&ys_g,&ye_g,&zs_g,&ze_g\
				   ,&user->ksp_its,&user->snes_its);
#else
  FORTRAN_NAME(processoldsolution)(&(xvec[zs_g-1][ys_g-1][xs_g-1])\
				   ,&xs_g,&xe_g,&ys_g,&ye_g,&zs_g,&ze_g\
				   ,&user->ksp_its,&user->snes_its);
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
  FORTRAN_NAME(CORRECTTIMESTEP)(&(dn[0]),&(dnh[0]),&(dnp[0]),&user->ierr,&user->indata.dt);
#else
  FORTRAN_NAME(correcttimestep)(&(dn[0]),&(dnh[0]),&(dnp[0]),&user->ierr,&user->indata.dt);
#endif

  /*
   * Deallocate memory
   */

  ierr = VecDestroy(dummy);CHKERRQ(ierr);

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
#define __FUNCT__ "ApplyASPC"
/* See ex6.c of SNES */
/* ------------------------------------------------------------------- */
/*
   Additive Schwartz preconditioner, using physics-based PC for diagonal
   inversion. 

   Input Parameters:
   ctx - optional user-defined context, as set by PCShellSetApply()
   y - input vector

   Output Parameter:
   z - preconditioned vector
*/
int ApplyASPC(void *ctx,Vec y,Vec z)
{
  AppCtx  *user = (AppCtx*)ctx;

  int	  ierr,k;
  PetscScalar  mone=-1.,omega=1.;

  PetscFunctionBegin;

  for (k=1;k<=user->aspc_its;++k) {

    /* Calculate residual rk=y-J*zk*/
    if (k == 1) {
      ierr = VecCopy(y,user->rk)            	        	     ;CHKERRQ(ierr);
    } else {
      ierr = MatMult(user->J ,z,user->rk)               	     ;CHKERRQ(ierr);
      /*ierr = VecView(user->rk,PETSC_VIEWER_STDOUT_WORLD)	     ;CHKERRQ(ierr);*/
      ierr = VecAYPX(&mone   ,y,user->rk)               	     ;CHKERRQ(ierr);
    }

    /* Calculate update dzk=P^-1(rk) */
    ierr = ApplyPC(user,user->rk,user->dzk)	   		     ;CHKERRQ(ierr);

    /* Correct solution: z = z + omega*dzk */
    ierr = VecAXPY(&omega,user->dzk,z)      	   		     ;CHKERRQ(ierr);

  }

  PetscFunctionReturn(0);

}


#undef __FUNCT__
#define __FUNCT__ "ApplyPC"
/* See ex6.c of SNES */
/* ------------------------------------------------------------------- */
/*
   Input Parameters:
   ctx - optional user-defined context, as set by PCShellSetApply()
   y - input vector

   Output Parameter:
   z - preconditioned vector

   Current does additive NKSchwartz
*/
int ApplyPC(void *ctx,Vec y,Vec x)
{
  AppCtx  *user = (AppCtx*)ctx;

  Field	  ***xvec,***yvec;
  int	  ierr,xs,ys,zs,xm,ym,zm,ze,ye,xe,mx,my,mz;

  PetscFunctionBegin;

  /*
   * Get pointers to vector data
   */
  ierr = DAVecGetArray(user->da, x, (void**)&xvec);CHKERRQ(ierr);
  ierr = DAVecGetArray(user->da, y, (void**)&yvec);CHKERRQ(ierr);
  
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
   * Compute initial condition over the locally owned part of the grid, including 
   * ghost cells (which are needed for restarting the calculation)
   */

#ifdef absoft
  FORTRAN_NAME(APPLYSHELLPC)(&(yvec[zs-1][ys-1][xs-1])\
			    ,&(xvec[zs-1][ys-1][xs-1])\
                            ,&xs,&xe,&ys,&ye,&zs,&ze);
#else
  FORTRAN_NAME(applyshellpc)(&(yvec[zs-1][ys-1][xs-1])\
			    ,&(xvec[zs-1][ys-1][xs-1])\
                            ,&xs,&xe,&ys,&ye,&zs,&ze);
#endif

  /* Restore vectors */
  ierr = DAVecRestoreArray(user->da, x, (void**)&xvec);CHKERRQ(ierr);
  ierr = DAVecRestoreArray(user->da, y, (void**)&yvec);CHKERRQ(ierr);

  PetscFunctionReturn(0);

}

#undef __FUNCT__
#define __FUNCT__ "SetupPreconditioner"
/* ------------------------------------------------------------------- */
/*
   Input Parameters:
   ctx - optional user-defined context, as set by PCShellSetSetUP()
*/
int SetupPreconditioner(void *ctx)
{
  AppCtx  *user = (AppCtx*)ctx;

  Field	  ***xvec;
  int	  ierr,i,j,k,xs,ys,zs,xm,ym,zm,ze,ye,xe,mx,my,mz;
  int     xs_g,ys_g,zs_g,ze_g,ye_g,xe_g,xm_g,ym_g,zm_g;
  Vec     localX;

  PetscFunctionBegin;

  /*
   * Scatter ghost points to local vector,using the 2-step process
   * DAGlobalToLocalBegin(),DAGlobalToLocalEnd().
   * By placing code between these two statements, computations can be
   * done while messages are in transition.
   */
  
  ierr = DAGetLocalVector(user->da,&localX);CHKERRQ(ierr);
  
  ierr = DAGlobalToLocalBegin(user->da,user->xk,INSERT_VALUES,localX);CHKERRQ(ierr);
  ierr = DAGlobalToLocalEnd  (user->da,user->xk,INSERT_VALUES,localX);CHKERRQ(ierr);
  
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
  FORTRAN_NAME(SETUPSHELLPC)(&(xvec[zs_g-1][ys_g-1][xs_g-1])\
                              ,&xs_g,&xe_g,&ys_g,&ye_g,&zs_g,&ze_g);
#else
  FORTRAN_NAME(setupshellpc)(&(xvec[zs_g-1][ys_g-1][xs_g-1])\
                              ,&xs_g,&xe_g,&ys_g,&ye_g,&zs_g,&ze_g);
#endif

  /* Restore vectors */
  ierr = DAVecRestoreArray(user->da, localX, (void**)&xvec);     CHKERRQ(ierr);
  ierr = DARestoreLocalVector(user->da,&localX);                 CHKERRQ(ierr);
  ierr = DALocalToGlobal(user->da,localX,INSERT_VALUES,user->xk);CHKERRQ(ierr);


  PetscFunctionReturn(0);
}
