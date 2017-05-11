#define _Dim 3
#define _Dimx 3

c Norms (L1 fast)

#define L1_NORM

c Gfortran

#ifdef gfortran
#define __GNU
#endif

c Precision (_Npg defined in Makefile)

#if _Npg == 8
#define _prec 4
#elif _Npg == 4 || _Npg ==2 || _Npg == 1
#define _prec 8
#endif

c Vectorization

#if _Npg == 1
#define VECTOR 
#else
#define VECTOR $OMP
#endif
