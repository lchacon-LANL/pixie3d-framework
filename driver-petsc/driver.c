static char help[] = "- Div (alpha T^beta Grad T) = 0.\n\
Different matrices for the Jacobian and the preconditioner.\n\
Demonstrates the use of matrix-free Newton-Krylov methods in conjucntion\n\
with a user-provided precondintioner. Input arguments are:\n\
	-snes_mf : use matrix-free Newton methods (default)\n\
	-snes_view : print the SNES data structure\n\
	-draw_pause <sec> : set time (in seconds) that the program pauses\n\
	-nox : disable all x-windows output\n";
/*T
 *    Concepts: SNES^matrix-free/finite-difference Jacobian methods
 *    Concepts: SNES^user-provided preconditioner;
 *    Processors: n (parallel)
 *    
T*/

/* 
 *    Include "petscsnes.h" so that we can use SNES solvers.  Note that this
 *    file automatically includes:
 *    petsc.h       - base PETSc routines   petscvec.h - vectors
 *    petscsys.h    - system routines       petscmat.h - matrices
 *    petscis.h     - index sets            petscksp.h - Krylov subspace methods
 *    petscviewer.h - viewers               petscpc.h  - preconditioners
 *    petscsles.h   - linear solvers
 *    
*/

#include "petscsnes.h"
#include "petscda.h"
#include "fortran.h"

typedef struct {
  Scalar var[NVAR];
} Field;

#ifdef absoft
extern void FORTRAN_NAME(EVALUATENONLINEARRESIDUAL) 
                      (Field*, Field*, int*, int*, int*, int*, int*, int*,
		       int*, int*, int*, int*, int*, int*);
extern void FORTRAN_NAME(INITIALIZECALCULATION) (Field*, int*, int*, int*, int*, int*, int*);
extern void FORTRAN_NAME(PROCESSOLDSOLUTION) (Field*,int*,int*,int*,int*,int*,int*,int*,int*);
extern void FORTRAN_NAME(FORTRANDESTROY) ();
#else
extern void FORTRAN_NAME(evaluatenonlinearresidual) 
                      (Field*, Field*, int*, int*, int*, int*, int*, int*,
		       int*, int*, int*, int*, int*, int*);
extern void FORTRAN_NAME(initializecalculation) (Field*, int*, int*, int*, int*, int*, int*);
extern void FORTRAN_NAME(processoldsolution) (Field*,int*,int*,int*,int*,int*,int*,int*,int*);
extern void FORTRAN_NAME(fortrandestroy) ();
#endif

/* User-defined application context */
typedef struct {
  DA	    da;            /* distributed array data structure */
  int       nmax,nwits,gmits,maxitnwt,maxitgm;
  PetscReal atol,rtol,tolgm;
} AppCtx;

#define POWFLOP 5 /* assume a pow() takes five flops */
#define DIFF_NORM 1.0e-10

/* User-defined routines */
extern int FormFunction    (SNES, Vec, Vec, void*);
extern int FormInitialGuess(SNES, Vec, void*);
extern int Monitor         (SNES, int, double, void*);

int main(int argc, char **argv)
{
	SNES	snes;	        /* SNES context */
	SLES    sles;           /* SLES context */
	KSP     ksp;            /* KSP context  */
	PC      pc;             /* PC context   */
	Vec	x, r;	        /* Vec x: solution, r: residual */
	Mat	J;
	int	ierr, its, L = 10,  N = 10, M = 10, i, size, rank, steps;
	PetscTruth	matrix_free, coloring;
	
	AppCtx	user;

	ierr = PetscInitialize(&argc, &argv, (char *)0, help);CHKERRQ(ierr);
	
	/* MPI_Comm_size(PETSC_COMM_WORLD, &size);
	   MPI_Comm_rank(PETSC_COMM_WORLD, &rank); */

	ierr = PetscOptionsGetInt(PETSC_NULL,"-nx",&L,PETSC_NULL);CHKERRQ(ierr);
	ierr = PetscOptionsGetInt(PETSC_NULL,"-ny",&N,PETSC_NULL);CHKERRQ(ierr);
	ierr = PetscOptionsGetInt(PETSC_NULL,"-nz",&M,PETSC_NULL);CHKERRQ(ierr);
	
	/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	                    Set problem parameters
           - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  */

	/*ierr = PetscOptionsGetDouble(PETSC_NULL,"-tleft",&user.tleft,PETSC_NULL);CHKERRQ(ierr);
	  ierr = PetscOptionsGetDouble(PETSC_NULL,"-tright",&user.tright,PETSC_NULL);CHKERRQ(ierr);
	  ierr = PetscOptionsGetDouble(PETSC_NULL,"-beta",&user.beta,PETSC_NULL);CHKERRQ(ierr);
	  ierr = PetscOptionsGetDouble(PETSC_NULL,"-theta",&user.theta,PETSC_NULL);CHKERRQ(ierr);
	  ierr = PetscOptionsGetDouble(PETSC_NULL,"-dt",&user.dt,PETSC_NULL);CHKERRQ(ierr); */

	ierr = PetscOptionsGetInt(PETSC_NULL,"-nmax",&user.nmax,PETSC_NULL);CHKERRQ(ierr);

	user.atol  = PETSC_DEFAULT;
	user.rtol  = 1e-4;
	user.tolgm = 1e-2;

	user.maxitgm = 15;
	user.maxitnwt= 15;
	
	
	/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	                    Create SNES context 
           - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  */

	ierr = SNESCreate(PETSC_COMM_WORLD, SNES_NONLINEAR_EQUATIONS, &snes);CHKERRQ(ierr);
	
	ierr = DACreate3d(PETSC_COMM_WORLD,DA_NONPERIODIC,DA_STENCIL_BOX,L,N,M,\
			  PETSC_DECIDE,PETSC_DECIDE,PETSC_DECIDE,\
			  NVAR,1,PETSC_NULL,PETSC_NULL,PETSC_NULL,&user.da);CHKERRQ(ierr);

	ierr = DACreateGlobalVector(user.da,&x);CHKERRQ(ierr);
	ierr = VecDuplicate(x,&r);CHKERRQ(ierr);
	
	/*ierr = SNESSetApplicationContext(snes,(void*)&user);CHKERRQ(ierr);
	ierr = SNESSetFunction(snes,r,FormFunction,PETSC_NULL);CHKERRQ(ierr);*/

	ierr = SNESSetFunction(snes,r,FormFunction,(void*)&user);CHKERRQ(ierr);
	
	/*ierr = SNESSetMonitor(snes,Monitor,(void*)&user,PETSC_NULL);CHKERRQ(ierr);*/

        /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                Customize nonlinear solver; set runtime options
           - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

	/* 
	   Set linear solver defaults for this problem. By extracting the
	   SLES, KSP, and PC contexts from the SNES context, we can then
	   directly call any SLES, KSP, and PC routines to set various options.
	*/

	/* SNESSetTolerances(snes,atol,rtol,stol,maxit,maxf); */

	ierr = SNESSetTolerances(snes,user.atol,user.rtol,PETSC_DEFAULT,user.maxitnwt,PETSC_DEFAULT);
               CHKERRQ(ierr);

	/* KSPSetTolerances(ksp,rtol,atol,dtol,maxits); */

	ierr = SNESGetSLES(snes,&sles); CHKERRQ(ierr);
	ierr = SLESGetKSP(sles,&ksp); CHKERRQ(ierr);
	ierr = SLESGetPC(sles,&pc); CHKERRQ(ierr);
	ierr = PCSetType(pc,PCNONE); CHKERRQ(ierr);
	ierr = KSPSetTolerances(ksp,user.tolgm,PETSC_DEFAULT,PETSC_DEFAULT,user.maxitgm);
	       CHKERRQ(ierr);

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
	
	return(0);
}

/* -------------------- User-defined Monitor() routine ------------- */
/*
   Monitor - User-defined monitoring routine that views the
   current iterate with an x-window plot.

   Input Parameters:
   snes - the SNES context
   its - iteration number
   norm - 2-norm function value (may be estimated)
   ctx - optional user-defined context for private data for the
         monitor routine, as set by SNESSetMonitor()

  Note:
  See the manpage for PetscViewerDrawOpen() for useful runtime option
  such as -nox to deactivate all x-window output.
*/
int Monitor(SNES snes, int its, double fnorm, void *ctx)
{
	int	ierr, p, my_rank;
	Vec	x, local_x, xslice;
	Scalar  ***x_array, **x_array_slice;
	PetscDraw	draw;
	AppCtx  *user = (AppCtx*)ctx;

	int     i,j,k,xs,ys,zs,xm,ym,zm,mx,my,mz;
	DA      da;


	DAGetInfo(user->da, PETSC_NULL, &mx, &my, &mz, PETSC_NULL, PETSC_NULL,
			PETSC_NULL, PETSC_NULL, PETSC_NULL, PETSC_NULL, PETSC_NULL);

	DAGetCorners(user->da, &xs, &ys, &zs, &xm, &ym, &zm);

	/*PetscPrintf(PETSC_COMM_WORLD, "iter = %d, SNES Function norm %g\n", its, fnorm);*/
	MPI_Comm_size(PETSC_COMM_WORLD, &p);
	MPI_Comm_rank(PETSC_COMM_WORLD, &my_rank);
	
	if (1) {

	  ierr = SNESGetSolution(snes, &x);CHKERRQ(ierr);

	  ierr = DAVecGetArray(user->da, x, (void**)&x_array);CHKERRQ(ierr);
	  ierr = DACreate2d(PETSC_COMM_WORLD,DA_NONPERIODIC,DA_STENCIL_STAR,mx,my,\
			    PETSC_DECIDE,PETSC_DECIDE,1,1,PETSC_NULL,PETSC_NULL,&da);
	  ierr = DAGetGlobalVector(da, &xslice);CHKERRQ(ierr);
	  ierr = DAVecGetArray(da, xslice, (void**)&x_array_slice);CHKERRQ(ierr);

	  /*
	  ierr = DAGetLocalVector(user->da,&local_x);CHKERRQ(ierr);
	  ierr = DAGlobalToLocalBegin(user->da,x,INSERT_VALUES,local_x);CHKERRQ(ierr);
	  ierr = DAGlobalToLocalEnd(user->da,x,INSERT_VALUES,local_x);CHKERRQ(ierr)
	  ierr = DAVecGetArray(user->da,local_x,(void**)&x_array);CHKERRQ(ierr);

	  ierr = DACreate2d(PETSC_COMM_WORLD,DA_NONPERIODIC,DA_STENCIL_STAR,mx,my,\
			    PETSC_DECIDE,PETSC_DECIDE,1,1,PETSC_NULL,PETSC_NULL,&da);
	  ierr = DAGetGlobalVector(da, &xslice);CHKERRQ(ierr);
	  ierr = DAGetLocalVector(user->da,&local_xslice);CHKERRQ(ierr);
	  ierr = DAGlobalToLocalBegin(user->da,xslice,INSERT_VALUES,local_xslice);CHKERRQ(ierr);
	  ierr = DAGlobalToLocalEnd(user->da,xslice,INSERT_VALUES,local_xslice);CHKERRQ(ierr)

	  ierr = DAVecGetArray(user->da,local_xslice,(void**)&x_array_slice);CHKERRQ(ierr);
	  */

	  for (k = zs; k < zs+1; k++) {
	      for (j = 0; j < my; j++) {
		   for (i = 0; i < mx; i++) {
		     x_array_slice[j][i] = x_array[k][j][i];
		   }
	      }
	  }
	  /*
		     PetscPrintf(PETSC_COMM_WORLD, "%lf ", x_array[k][j][i]);
		   }
		   PetscPrintf(PETSC_COMM_WORLD, "\n");
	      }
	      PetscPrintf(PETSC_COMM_WORLD, "\n");
	  }
	  */

	  ierr = DAVecRestoreArray(da,xslice,(void**)&x_array_slice);CHKERRQ(ierr);

	  ierr = PetscViewerDrawGetDraw(PETSC_VIEWER_DRAW_(PETSC_COMM_WORLD), 0, &draw);CHKERRQ(ierr);
	  ierr = PetscDrawSetDoubleBuffer(draw);CHKERRQ(ierr);
	  ierr = VecView(xslice, PETSC_VIEWER_DRAW_(PETSC_COMM_WORLD));CHKERRQ(ierr);

	  /*
	  ierr = DAVecRestoreArray(user->da,local_x,(void**)&x_array);CHKERRQ(ierr);
	  ierr = DARestoreLocalVector(user->da,&local_x);CHKERRQ(ierr);
	  */
	  ierr = DAVecRestoreArray(user->da, x,(void**)&x_array);CHKERRQ(ierr);
	  }

	

	return (0);
}

/* --------------------  Form initial approximation ----------------- */
#undef __FUNCT__
#define __FUNCT__ "FormInitialGuess"
int FormInitialGuess(SNES snes,Vec X,void *ptr)
{
	AppCtx	*user = (AppCtx*)ptr;
	Field	***xvec;
	Scalar	hx, hy, hz, xp, yp;
	int	ierr, i, j, k, xs, ys, zs, xm, ym, zm,  mx, my, mz;
	int     ze,ye,xe;

	int     xs_g, ys_g, zs_g, ze_g,ye_g,xe_g, xm_g, ym_g, zm_g;


	/*
	DAGetInfo(user->da, PETSC_NULL, &mx, &my, &mz, PETSC_NULL, PETSC_NULL,
			PETSC_NULL, PETSC_NULL, PETSC_NULL, PETSC_NULL, PETSC_NULL);

	hx = 1.0/(Scalar)(mx - 1);
	hy = 1.0/(Scalar)(my - 1);
	hz = 1.0/(Scalar)(mz - 1);
	*/

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

	return(0);
}

/* --------------------  Process old time solution ----------------- */
#undef __FUNCT__
#define __FUNCT__ "ProcessOldSolution"
int ProcessOldSolution(SNES snes,Vec X,void *ptr)
{
	AppCtx	*user = (AppCtx*)ptr;
	Field	***xvec;
	Scalar	hx, hy, hz, xp, yp;
	int	ierr, i, j, k, xs, ys, zs, xm, ym, zm,  mx, my, mz;
	int     ze,ye,xe;
	
	int     xs_g, ys_g, zs_g, ze_g,ye_g,xe_g, xm_g, ym_g, zm_g;

	/*
	DAGetInfo(user->da, PETSC_NULL, &mx, &my, &mz, PETSC_NULL, PETSC_NULL,
			PETSC_NULL, PETSC_NULL, PETSC_NULL, PETSC_NULL, PETSC_NULL);
	*/

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
	return(0);
}

/* --------------------  Evaluate Function G(x,t) --------------------- */
#undef __FUNCT__
#define __FUNCT__ "FormFunction"
int FormFunction(SNES snes,Vec X,Vec F,void* ptr)
{
	AppCtx  *user = (AppCtx*)ptr;
	/*AppCtx  *user;*/
	
	int     ierr,i,j,k,xs,ys,zs,xm,ym,zm,mx,my,mz;
	Field   ***x,***f;
	Vec     localX, localX_old;
	double  dt, theta;

	int     ye,xe,ze;
	int     ivar,my_rank, xs_g,ys_g,zs_g,ze_g,ye_g,xe_g,xm_g,ym_g,zm_g;

	/*PetscFunctionBegin;*/

	/*dt = user->dt;		theta = user->theta;*/

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
	/*PetscFunctionReturn(0);*/

	return (0);

}
