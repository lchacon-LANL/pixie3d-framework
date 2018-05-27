#define _Dim 3
#define _Dimx 2

c Norms (L1 fast)

#define L1_NORM

c Gfortran

#ifdef gfortran
#define __GNU
#endif

c Vectorization

#if _Npg == 1
#define VECTOR 
#else
#define VECTOR $OMP
#endif
