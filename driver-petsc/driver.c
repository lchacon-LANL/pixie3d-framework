static char help[] = "Options:
	-snes_mf : use matrix-free Newton methods (default)\n\
	-user_precond: use user-defined preconditioner (-snes_mf)\n\
	-nx <nx> -ny <ny> -nz <nz> : grid size of each direction\n\
        -nmax <nmax>: max. iteration # of time steps\n";

/*T
 *    Concepts: SNES^matrix-free/finite-difference Jacobian methods
 *    Concepts: SNES^user-provided preconditioner;
 *    Processors: n (parallel)
 *    
T*/

#include "petscsnes.h"
#include "petscda.h"
#include "fortran.h"
#include <stdio.h>
#include <string.h>

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
  PetscReal dt0;
  int       nxd;
  int       nyd;
  int       nzd;
  int       maxitnwt;
  int       maxksp;
  int       maxitgm;
  int       method;
  int       global;
  int       iguess;
  int       bcs[6];
} input_CTX;

#ifdef absoft
extern void FORTRAN_NAME(EVALUATENONLINEARRESIDUAL) 
                      (Field*, Field*, int*, int*, int*, int*, int*, int*,
		       int*, int*, int*, int*, int*, int*);
extern void FORTRAN_NAME(INITIALIZECALCULATION) (Field*, int*, int*, int*, int*, int*, int*);
extern void FORTRAN_NAME(PROCESSOLDSOLUTION) (Field*,int*,int*,int*,int*,int*,int*,int*,int*);
extern void FORTRAN_NAME(FORTRANDESTROY) ();
extern void FORTRAN_NAME(READINPUTFILE) (input_CTX*);
#else
extern void FORTRAN_NAME(evaluatenonlinearresidual) 
                      (Field*, Field*, int*, int*, int*, int*, int*, int*,
		       int*, int*, int*, int*, int*, int*);
extern void FORTRAN_NAME(initializecalculation) (Field*, int*, int*, int*, int*, int*, int*);
extern void FORTRAN_NAME(processoldsolution) (Field*,int*,int*,int*,int*,int*,int*,int*,int*);
extern void FORTRAN_NAME(fortrandestroy) ();
extern void FORTRAN_NAME(readinputfile) (input_CTX*);
#endif


typedef struct {
  DA	    da;
  input_CTX indata;
  int       nmax,nwits,gmits;
} AppCtx;

/* User-defined routines */
extern int ReadInputNInitialize(input_CTX*, char*);
extern int FormFunction    (SNES, Vec, Vec, void*);
extern int FormInitialGuess(SNES, Vec, void*);
extern int Monitor         (SNES, int, double, void*);
extern int MatrixFreePreconditioner(void*,Vec,Vec);

#if !defined(lahey)
int main(int argc, char **argv)
#else
int MAIN__(int argc, char **argv)
#endif
{
  SNES     snes;	        /* SNES context */
  SLES     sles;           /* SLES context */
  KSP      ksp;            /* KSP context  */
  PC       pc;             /* PC context   */
  Vec	   x, r;	        /* Vec x: solution, r: residual */
  AppCtx   user;
  
  int     ierr;
  int     nxd = 32;         /* x-drirection grid size, i.e., nx */
  int     nyd = 32;         /* y-drirection grid size, i.e., ny */
  int     nzd = 1;          /* z-drirection grid size, i.e., nz */
  int     i, steps;
  int     size;           /* num. of processors used */
  int     my_rank;        /* id of my processor */

  PetscReal atol,rtol,tolgm;
  int       maxitnwt,maxitgm;
  char      BC[256];

  Mat                J;     /* Jacobian matrix for FD approximation */
  PetscTruth         matrix_free,user_precond;

  
  
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Begin program
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  */	
  ierr = PetscInitialize(&argc, &argv, (char *)0, help);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Set problem parameters
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  */
  ierr = ReadInputNInitialize(&user.indata, BC);CHKERRQ(ierr);
  
  ierr = PetscOptionsGetInt(PETSC_NULL,"-nmax",&user.nmax,PETSC_NULL);CHKERRQ(ierr);
  
  nxd = user.indata.nxd;
  nyd = user.indata.nyd;
  nzd = user.indata.nzd;
  
  atol = PETSC_DEFAULT;
  rtol = user.indata.rtol;
  tolgm = user.indata.tolgm;
  maxitnwt = user.indata.maxitnwt;
  if (maxitnwt == 0)
    maxitnwt = (int) PetscMax((1.5*log(rtol)/log(tolgm)),10.);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Create SNES context 
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  */

  ierr = SNESCreate(PETSC_COMM_WORLD, &snes);CHKERRQ(ierr);
  
  ierr = DACreate3d(PETSC_COMM_WORLD,*BC,DA_STENCIL_BOX, nxd,nyd,nzd,\
		    PETSC_DECIDE,PETSC_DECIDE,PETSC_DECIDE,\
		    NVAR,1,PETSC_NULL,PETSC_NULL,PETSC_NULL,&user.da);CHKERRQ(ierr);

  ierr = DACreateGlobalVector(user.da,&x);CHKERRQ(ierr);
  ierr = VecDuplicate(x,&r);CHKERRQ(ierr);
  
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
    ierr = SNESGetSLES(snes,&sles); CHKERRQ(ierr);
    ierr = SLESGetKSP (sles,&ksp ); CHKERRQ(ierr);
    ierr = SLESGetPC  (sles,&pc  ); CHKERRQ(ierr);
    ierr = KSPSetTolerances(ksp,tolgm,PETSC_DEFAULT,PETSC_DEFAULT,maxitgm);CHKERRQ(ierr);
    ierr = PetscOptionsHasName(PETSC_NULL,"-user_precond",&user_precond);CHKERRQ(ierr);
    if (user_precond) {
      /*
      SETERRQ(1, "user-defined preconditoner is not defined yet!!!");
      */
      ierr = PCSetType(pc,PCSHELL);CHKERRQ(ierr);
      ierr = PCShellSetApply(pc,MatrixFreePreconditioner,PETSC_NULL);CHKERRQ(ierr);
      ierr = PCShellSetName(pc,"user-defined preconditioner");CHKERRQ(ierr);
    } else {
      ierr = PCSetType(pc,PCNONE); CHKERRQ(ierr);
    }
  } else {
    SETERRQ(1, "You shoud use \"matrix-free method\", i.e., -snes_mf");
  }

  /* Set SNES run time options */
  
  ierr = SNESSetFromOptions(snes);CHKERRQ(ierr);
  
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Initialize calculation 
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  
  ierr = FormInitialGuess(snes, x, &user);CHKERRQ(ierr);
  
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Time stepping 
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  
  /* for (steps = 1; (steps < 100) && (norm_2 > DIFF_NORM); steps++) { */
  for (steps = 0; steps < user.nmax; steps++) {
    
    ierr = ProcessOldSolution(snes,x,&user);CHKERRQ(ierr);
    
    ierr = SNESSolve(snes,x,&user.nwits);CHKERRQ(ierr);
    
    ierr = SNESGetNumberLinearIterations(snes,&user.gmits);CHKERRQ(ierr);
    
    /*PetscPrintf(PETSC_COMM_WORLD,"Steps = %d Number of Newton iterations = %d\n"\
      ,steps, user.nwits);
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

  VecDestroy(x);           
  VecDestroy(r);
  SNESDestroy(snes);
  DADestroy(user.da);
  PetscFinalize();
  
  PetscFunctionReturn(0);
}

/* -------------------- Read input file routine ------------- */
#undef __FUNCT__
#define __FUNCT__ "ReadInputNInitialize"
int ReadInputNInitialize(input_CTX* data, char* BC)
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
    strcpy(BC, "DA_XYZPERIODIC");
  } else if (data->bcs[0] && data->bcs[2] && !data->bcs[4]) {
    strcpy(BC, "DA_XYPERIODIC");
  } else if (data->bcs[0] && !data->bcs[2] && data->bcs[4]) {
    strcpy(BC, "DA_XZPERIODIC");
  } else if (data->bcs[0] && !data->bcs[2] && !data->bcs[4]) {
    strcpy(BC, "DA_XPERIODIC");
  } else if (!data->bcs[0] && data->bcs[2] && data->bcs[4]) {
    strcpy(BC, "DA_YZPERIODIC");
  } else if (!data->bcs[0] && data->bcs[2] && !data->bcs[4]) {
    strcpy(BC, "DA_YPERIODIC");
  } else if (!data->bcs[0] && !data->bcs[2] && data->bcs[4]) {
    strcpy(BC, "DA_ZPERIODIC");
  } else {
    strcpy(BC, "DA_NONPERIODIC");
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

/* --------------------  Form initial approximation ----------------- */
#undef __FUNCT__
#define __FUNCT__ "FormInitialGuess"
int FormInitialGuess(SNES snes,Vec X,void *ptr)
{
  AppCtx	*user = (AppCtx*)ptr;
  Field	***xvec;
  int	ierr,i,j,k,xs,ys,zs,xm,ym,zm,mx,my,mz;
  int     ze,ye,xe;
  int     xs_g,ys_g,zs_g,ze_g,ye_g,xe_g,xm_g,ym_g,zm_g;
  

  PetscFunctionBegin;

  /*
   * Initialize counters
   */

  user->nwits = 0;
  user->gmits = 0;
  
  /*
   * Get pointers to vector data
   */
  DAVecGetArray(user->da, X, (void**)&xvec);
  
  /*
   * Get local grid boundaries
   */
  DAGetCorners(user->da, &xs, &ys, &zs, &xm, &ym, &zm);
  
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
   * Compute function over the locally owned part of the grid
   */

#ifdef absoft
  FORTRAN_NAME(INITIALIZECALCULATION)(&(xvec[zs_g-1][ys_g-1][xs_g-1])\
                                           ,&xs_g,&xe_g,&ys_g,&ye_g,&zs_g,&ze_g);
#else
  FORTRAN_NAME(initializecalculation)(&(xvec[zs_g-1][ys_g-1][xs_g-1])\
                                           ,&xs_g,&xe_g,&ys_g,&ye_g,&zs_g,&ze_g);
#endif

  /* Restore vectors */
  DAVecRestoreArray(user->da, X, (void**)&xvec);

  PetscFunctionReturn(0);
}

/* --------------------  Process old time solution ----------------- */
#undef __FUNCT__
#define __FUNCT__ "ProcessOldSolution"
int ProcessOldSolution(SNES snes,Vec X,void *ptr)
{
  AppCtx	*user = (AppCtx*)ptr;
  Field	***xvec;
  PetscScalar	hx, hy, hz, xp, yp;
  int	ierr,i,j,k,xs,ys,zs,xm,ym,zm,mx,my,mz;
  int     ze,ye,xe;
  int     xs_g,ys_g,zs_g,ze_g,ye_g,xe_g,xm_g,ym_g,zm_g;
  
  PetscFunctionBegin;
  /*
   * Get pointers to vector data
   */
  DAVecGetArray(user->da, X, (void**)&xvec);
  
  /*
   * Get local grid boundaries
   */
  DAGetCorners(user->da, &xs, &ys, &zs, &xm, &ym, &zm);
  
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
	
#ifdef absoft
  FORTRAN_NAME(PROCESSOLDSOLUTION)(&(xvec[zs_g-1][ys_g-1][xs_g-1])\
				   ,&xs_g,&xe_g,&ys_g,&ye_g,&zs_g,&ze_g\
				   ,&user->gmits,&user->nwits);
#else
  FORTRAN_NAME(processoldsolution)(&(xvec[zs_g-1][ys_g-1][xs_g-1])\
				   ,&xs_g,&xe_g,&ys_g,&ye_g,&zs_g,&ze_g\
				   ,&user->gmits,&user->nwits);
#endif

  /* Restore vectors */
  DAVecRestoreArray(user->da, X, (void**)&xvec);
  
  PetscFunctionReturn(0);
}

/* --------------------  Evaluate Function G(x,t) --------------------- */
#undef __FUNCT__
#define __FUNCT__ "FormFunction"
int FormFunction(SNES snes,Vec X,Vec F,void* ptr)
{
  AppCtx  *user = (AppCtx*)ptr;
  
  int     ierr,i,j,k,xs,ys,zs,xm,ym,zm,mx,my,mz;
  Field   ***x,***f;
  Vec     localX;
  double  dt, theta;
  
  int     ye,xe,ze;
  int     ivar,my_rank, xs_g,ys_g,zs_g,ze_g,ye_g,xe_g,xm_g,ym_g,zm_g;

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
