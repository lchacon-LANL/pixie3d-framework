/* Fortran <--> C/C++ interfacing stuff */
/* $Id$ */

/* this header file is idempotent */
#ifndef FORTRAN_H_SEEN
#define FORTRAN_H_SEEN

/*****************************************************************************/

/*
 * C/C++ compatability:
 */

/*
 * Use this in prototypes like this:  extern FORTRAN_FUNCTION void foo(...)
 *
 * At present, this is set up to tell a C++ compiler that  foo()  uses
 * a C-compatible calling convention.
 */
#ifdef __cplusplus
  #define FORTRAN_FUNCTION	"C"
#else
  #define FORTRAN_FUNCTION	/* empty */
#endif

/*****************************************************************************/

/* array subscripting offset, i.e. C "arr[k]" is Fortran "arr(k+?)" */
#define FORTRAN_INDEX_ORIGIN	1

/*****************************************************************************/

/* C type of Fortran integer/logical variables */
/* also actual integers used for Fortran logical .true. and .false. */

/*
 * FIXME: these are what I (JT) used on a 32-bit SGI system with the
 *	  SGI C and Fortran compilers in 1992-1993; they should be
 *	  checked to see if they're still valid for current compilers
 *	  and/or for 64-bit systems.
 */
#if   defined(sgi) || defined(SGI) || defined(__sgi__) || defined(__SGI__)
  #define FORTRAN_INTEGER_IS_INT	TRUE
  typedef int integer;
  typedef unsigned int logical;
  #define FORTRAN_LOGICAL_TRUE	1
  #define FORTRAN_LOGICAL_FALSE	0

/* see FIXME above for validity of these */
#elif defined(alpha) || defined(ALPHA) || defined(__alpha__) || defined(__ALPHA__)
  #define FORTRAN_INTEGER_IS_INT	TRUE
  typedef int integer;
  typedef unsigned int logical;
  #define FORTRAN_LOGICAL_TRUE	1
  #define FORTRAN_LOGICAL_FALSE	0

#elif defined(sun) || defined(SUN) || defined(__sun__) || defined(__SUN__)
  #define FORTRAN_INTEGER_IS_INT	TRUE
  typedef int integer;
  typedef unsigned int logical;
  #define FORTRAN_LOGICAL_TRUE	1
  #define FORTRAN_LOGICAL_FALSE	0

#elif defined(__linux__) && defined(__i386__) && defined(__GNUC__)
  #define FORTRAN_INTEGER_IS_INT	TRUE
  typedef int integer;
  typedef unsigned int logical;
  #define FORTRAN_LOGICAL_TRUE	1
  #define FORTRAN_LOGICAL_FALSE	0

#else
  #error "don't know Fortran integer/logical datatypes for this system!"
#endif

/* old (backwards compatible) names for Fortran integers/logicals */
typedef integer fortran_integer_t;
typedef logical fortran_logical_t;

/*****************************************************************************/

/*
 * Names of Fortran routines are often altered by the compiler/loader.  The
 * following macro should be used to call a Fortran routine from C code, i.e.
 *	call sgefa(...)			-- Fortran code
 *	FORTRAN_NAME(sgefa)(...);	-- C code to do the same thing
 *
 * Unfortunately, the "alterations" are generally at the token level, and this
 * can't be done portably in pre-ANSI C.  In ANSI C, the preprocessor "token
 * splicing" facility is designed to handle just this sort of thing, but in
 * pre-ANSI C we have to use rather ugly system-dependent hacks of the sort
 * exemplified below.
 */

/* see FIXME above for validity of these */
 #if   defined(sgi) || defined(SGI) || defined(__sgi__) || defined(__SGI__)
  /* C code should reference Fortran names in lower case */
  #ifdef __STDC__
    #define FORTRAN_NAME(n_)	n_ ## _
  #else
    #define FORTRAN_NAME(n_)	n_/**/_
  #endif

/* see FIXME above for validity of these */
#elif defined(alpha) || defined(ALPHA) || defined(__alpha__) || defined(__ALPHA__)
  /* C code should reference Fortran names in lower case */
  #ifdef __STDC__
    #define FORTRAN_NAME(n_)	n_ ## _
  #else
    #define FORTRAN_NAME(n_)	n_/**/_
  #endif

#elif defined(sun) || defined(SUN) || defined(__sun__) || defined(__SUN__)
  /* C code should reference Fortran names in lower case */
  #ifdef __STDC__
    #define FORTRAN_NAME(n_)	n_ ## _
  #else
    #define FORTRAN_NAME(n_)	n_/**/_
  #endif

#elif defined(__linux__) && defined(__i386__) && defined(__GNUC__) && !defined(absoft) 
  /* C code should reference Fortran names in lower case */
  #ifdef __STDC__
    #define FORTRAN_NAME(n_)	n_ ## _
  #else
    #define FORTRAN_NAME(n_)	n_/**/_
  #endif
#elif defined(__linux__) && defined(__i386__) && defined(__GNUC__) && defined(absoft) 
  /* C code should reference Fortran names in lower case */
  #define FORTRAN_NAME(n_)	n_
#else
  #error "don't know Fortran integer/logical datatypes for this system!"
#endif 

/*****************************************************************************/

double **
d2d(int size1, int size2)
{
	double *ptr;
	double **ret;
	int n = size1*size2;
	int i,j,k;

/*
** Allocate a block of memory of length size1*size2.
** Each memory location is of type double.
** In FORTRAN terms this is a one dimensional array 0:(size1-1)*(size2-1).
** This will form the columns of the 2-d array.
** Check for allocation failure.
** Note that calloc zeros the memory area.
*/
	ptr = (double *)calloc(n,sizeof(double));

	/*	if (ptr == (double *)NULL)
		error("First stage memory allocation failure in c2d", EMEM); */

/*
** Allocate a block of memory of length size1.
** Each memory location is of type double *.
** In FORTRAN terms this means that each location
** can hold an entire one dimensional array (of arbitrary length)
** of type double.
** malloc does not zero the memory area, but this is ok because
** these elements are assigned to below.
*/
	ret = (double **)malloc(size2*sizeof(double *));

	/* if (ret == (dcomplex **)NULL)
	   error("Second stage memory allocation failure in d2d", EMEM); */

/*
** The body of the loop assigns the address of the
** 0,size,2*size,...,size*size element of ptr to the 0,1,2,...,size
** element of ret. Thus forming a valid 2-d array each of whose elements
** is of type double.
*/
	for (i=0,j=0; j<=size2; i+=size1,j++)
		ret[j] = &ptr[i];

	return(ret);
}


void
free_d2d(ptr, size1, size2)
double **ptr;
int size1,size2;
{
	int i,j;

	for (i=size1,j=0; j<size2; i-=size1,j++)
		if (ptr[j] != (double *)NULL)
				free((void *)ptr[j]);

	if (ptr != (double **)NULL)
		free((void *)ptr);

	return;
}

#endif	/* FORTRAN_H_SEEN */
