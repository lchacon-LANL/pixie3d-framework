#ifdef _DP
#define c_float_kind  c_double
#define float_p double
#define float2_p double2
#define _fconst(t)  t
#define _abs(t) fabs(t)
#define _max(t,s) fmax(t,s)
#define _min(t,s) fmin(t,s)
#define _rsqrt(t) rsqrt(t)
#define _sqrt(t) __dsqrt_rn(t)
#define _fma(a,b,c) __fma_rn(a,b,c)
#define _floor(t) floor(t)
#define _rcp(a) __drcp_rn(a)
#define _div(a,b) __ddiv_rn(a,b)
#define _acos(a) acos(a)
#define _cos(a) cos(a)
#else
#define c_float_kind  c_float
#define float2_p float2
#define float_p float
#define _fconst(t)  t##f
#define _abs(t) fabsf(t)
#define _max(t,s) fmaxf(t,s)
#define _min(t,s) fminf(t,s)
#define _rsqrt(t) rsqrtf(t)
#define _sqrt(t) __fsqrt_rn(t)
#define _fma(a,b,c) __fmaf_rn(a,b,c)
#define _floor(t) floorf(t)
#define _rcp(a) __frcp_rn(a)
#define _div(a,b) __fdiv_rn(a,b)
#define _acos(a) acosf(a)
#define _cos(a) cosf(a)
#endif
