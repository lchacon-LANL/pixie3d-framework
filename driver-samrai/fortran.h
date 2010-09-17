/* Fortran <--> C/C++ interfacing stuff */
/* $Id: fortran.h,v 1.2 2007/06/21 16:21:47 lchacon Exp $ */

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

#if defined(absoft) 
  /* C code should reference Fortran names verbatim      */
  #define FORTRAN_NAME(n_)	n_
#else
  /* C code should reference Fortran names in lower case */
  #ifdef __STDC__
    #define FORTRAN_NAME(n_)	n_ ## _
  #else
    #define FORTRAN_NAME(n_)	n_/**/_
  #endif
#endif 

/*****************************************************************************/

 #endif	/* FORTRAN_H_SEEN */
