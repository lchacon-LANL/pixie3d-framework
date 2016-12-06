#include "cprecision.h"

#include "util/cuPrintf.cu"
//typedef unsigned int   int;
typedef unsigned short ushort;
typedef unsigned char  uchar;

#define binSize_v (13) //(4)
#define UMUL(a, b) ( (a) * (b) ) 
#define UMAD(a, b, c) ( UMUL((a), (b)) + (c) )  
#define WARP_SIZE 32
#define SHARED_MEMORY_BANKS 16
#define RHO_BIN_COUNT 32
#define GRID_SIZE 64 //16 //4 //one block first(no partial sum of rho)
#define G_SHARED_BANKS 4
#define G_WARP 8 //4 //(G_SHARED_BANKS>>1) //8 //4 //2
#define LOG_WARP 5
#define LOG_Nx 5
#define BLOCK_SIZE (SHARED_MEMORY_BANKS*G_SHARED_BANKS) //(WARP_SIZE*3) //power of 2 (>32) careful
#define SMEMSIZE (RHO_BIN_COUNT)*BLOCK_SIZE //12k float_p shared array s_rho
#define MS (32*2*GRID_SIZE) //nx=32 for d_j; 2 species

#define maxThreads (64)  // number of threads per block
#define maxBlocks  (16)

#ifdef _USE_TEX
texture<float_p> d_xxTex;
texture<float_p> d_ETex;
#define FETCH(t, i) tex1Dfetch(t##Tex, i)
#else
#define FETCH(t, i) t[i]
#endif

#define MAXNSP (4)
__constant__ float dc_q_m_sp[MAXNSP];
__constant__ float_p dc_q_m[MAXNSP];
__constant__ float_p dc_m[MAXNSP];
__constant__ float_p dc_q[MAXNSP];
__constant__ int  dc_acNpcls[MAXNSP+1];
__constant__ int  dc_Npcls[MAXNSP];
__constant__ float_p dc_vth[MAXNSP];

__constant__ float etolr_sp;
__constant__ float etola_sp;
__constant__ float hx_sp;
__constant__ float ihx_sp;

__constant__ float_p etolr;
__constant__ float_p etola;
__constant__ float_p eps;
__constant__ float_p PI;
__constant__ float_p LL;
__constant__ float_p hx;
__constant__ float_p ihx;
__constant__ float_p xf0;
__constant__ float_p xx0;
__constant__ int  Nx;
__constant__ int  totNp;
__constant__ int Nsp;

//new 
struct Particles{
  float_p *d_x_n;  //xn
  float_p *d_x_np; //xnp
  float_p *d_v_n;  //xn
  float_p *d_v_np; //xnp
  int *d_ci_n; //cin
  int *d_ci_np;//cinp
} pcles;

#ifndef MIN
#define MIN(x,y) ((x < y) ? x : y)
#endif

template <class T>
__device__ T sign(T a)
{ 
  T b = (T) 1;
  return a<0 ? -b : b;
}


unsigned int nextPow2( unsigned int x ) {
    --x;
    x |= x >> 1;
    x |= x >> 2;
    x |= x >> 4;
    x |= x >> 8;
    x |= x >> 16;
    return ++x;
}

//Cuda  C programming guide v4.0 B.11 p119 
#ifdef _DP

__device__ void atomicAdd2(float_p * address, float_p val) { 
   unsigned long long int* address_as_ull = (unsigned long long int*)address; 
   unsigned long long int old, assumed, val2;
   do{ 
     assumed = * address_as_ull; //old; 
     val2 = __double_as_longlong(val +__longlong_as_double(assumed));
     old = atomicCAS(address_as_ull, assumed, val2);
   } while (assumed != old); 
   //return __longlong_as_double(old); 
}
/*
__device__ void atomicAdd2(float_p * address, float_p val) { 
   unsigned long long int* address_as_ull = (unsigned long long int*)address; 
   unsigned long long int old=*address_as_ull, assumed; 
   do{ 
     assumed = old; 
     old = atomicCAS(address_as_ull, assumed, __double_as_longlong(val +__longlong_as_double(assumed)));
   } while (assumed != old); 
   //return __longlong_as_double(old); 
}
*/
#endif

/* for testing */
__global__ void device_greetings(float_p *d_xx,Particles pcles,float_p dt)
{
  for(int pi=UMAD(blockIdx.x, blockDim.x, threadIdx.x); pi<totNp; pi+=UMUL(blockDim.x,gridDim.x) ){
  float_p x_n   = pcles.d_x_np [pi];
  float_p v_n   = pcles.d_v_np [pi];
  float_p x_np;
  x_np = x_n+dt;
  pcles.d_x_np[pi] = x_np;


  float_p v_np;
  v_np = v_n + dt;
  pcles.d_v_np[pi] = v_np;

   int ci_n  = pcles.d_ci_n [pi];
   int ci_np;
   ci_np = ci_n;
   pcles.d_ci_np[pi]= ci_np;
  }

}
/* end device_greetings() */

__device__ inline uchar select_species_dc(int pi)
{
  uchar isp = 0;
  #pragma unroll
  for(int i=Nsp-1; i>=0; i--){
    isp += min(pi/dc_acNpcls[i+1],1);
  }
  return isp;  
}
/* end select_species_dc() */

__device__ inline uchar select_species(int pi,int *acNpcls)
{
  uchar isp = 0;
  for(int i=Nsp-1; i>=0; i--){
    isp += min(pi/acNpcls[i+1],1);
  }
  return isp;  
}
/* end select_species() */

//particle B.C.
__device__ inline float_p par_bc_pcles(float_p xp)
{
#if   _SPORDER == 1  
  return fmod((xp-xx0+LL),LL)+xx0;  
#elif _SPORDER == 2
  xp = fmod((xp-xf0+LL),LL)+xf0;  //[xf0,xf0+LL)
  int is_xf0 = (xp==xf0);
  return xp + (float_p)(is_xf0)*LL; //(xf0,xf0+LL] 
#endif

}
/* end par_bc_pcles() */

/* locate the particle cell index */
__device__ inline int locate(float_p xp)
{
#if   _SPORDER == 1  
  return floor((xp-xx0)*ihx);
#elif _SPORDER == 2
  return ceil((xp-xf0)*ihx);
#endif
}
/* end locate() */

/* find charge density on shared memeory*/
__device__ inline void s_gather(float_p *garr, float_p q, int i, float_p x)
{
#if   _SPORDER == 1  
  garr[UMUL(i  ,BLOCK_SIZE)] += (_fconst(1.0)-x)*q;
  garr[UMUL((i+1)&(Nx-1),BLOCK_SIZE)] += x *q;
#elif _SPORDER == 2
  garr[UMUL(i-1         ,BLOCK_SIZE)]  += _fconst(0.5)*(_fconst(0.5)-x)*(_fconst(0.5)-x)*q;
  garr[UMUL(i&(Nx-1)    ,BLOCK_SIZE)]  += (_fconst(0.75)-x*x)*q;
  garr[UMUL((i+1)&(Nx-1),BLOCK_SIZE)]  += _fconst(0.5)*(_fconst(0.5)+x)*(_fconst(0.5)+x)*q;
#endif
}
/* end s_gather() */

/* find current density on shared memeory*/
__device__ inline void s_gatherj(float_p *garr, float_p q, int i, float_p x)
{
#if   _SPORDER == 2
  garr[UMUL(i-1,BLOCK_SIZE)] += (_fconst(1.0)-x)*q;
  //  garr[UMUL(i  ,BLOCK_SIZE)] +=     x *q;
  garr[UMUL(i&(Nx-1),BLOCK_SIZE)] +=     x *q;
#endif
}
/* end s_gatherj() */

/* find charge density */
__device__ void gather(float_p *garr, float_p q, int i, float_p x)
{
  garr[i]  += (_fconst(1.0)-x)*q;
  garr[i+1]+= x*q;
}
/* end gather() */

//find current density
__global__ void find_current2(Particles pcles, int isp, float_p *d_j)
{
  const int threadPos = threadIdx.x;
  __shared__ float_p s_j[SMEMSIZE];
  //initialize shared memory
  for(int xi=0; xi<RHO_BIN_COUNT; xi++ ){ //careful
    s_j[threadIdx.x + xi*BLOCK_SIZE ] = 0.;
  }
  __syncthreads();
  
  for(int pi=dc_acNpcls[isp]+UMAD(blockIdx.x, blockDim.x, threadIdx.x); pi<dc_acNpcls[isp+1]; pi+=UMUL(blockDim.x,gridDim.x) ){
    float_p x_n = pcles.d_x_n [pi];
    float_p x_np= pcles.d_x_np[pi];
    float_p vh  = (pcles.d_v_np[pi]+pcles.d_v_n[pi])*0.5;
    int  i_n  = pcles.d_ci_n [pi];
    int  i_np = pcles.d_ci_np[pi];
#ifdef _OFFSET
    s_gatherj(s_j+threadPos,vh,i_n ,(x_n *ihx));  
    s_gatherj(s_j+threadPos,vh,i_np,(x_np*ihx));      
#else
    s_gatherj(s_j+threadPos,vh,i_n ,(x_n *ihx-i_n +1.5));  
    s_gatherj(s_j+threadPos,vh,i_np,(x_np*ihx-i_np+1.5));  
#endif
  }
  __syncthreads();
  //reduce
  if(threadIdx.x < RHO_BIN_COUNT){ //careful
    float_p *s_rhoBase = s_j + UMUL(threadIdx.x, BLOCK_SIZE);
    float_p sum = 0;
    int pos;
    for(int j=0; j<G_SHARED_BANKS; j++){
      pos = threadIdx.x & (SHARED_MEMORY_BANKS-1);//to avoid bank conflict
      for(int i=0; i<SHARED_MEMORY_BANKS; i++){
	sum += s_rhoBase[j*SHARED_MEMORY_BANKS+pos];
	pos = (pos + 1) & (SHARED_MEMORY_BANKS-1);
      }
    }
    
    d_j[(isp*gridDim.x+blockIdx.x)*(RHO_BIN_COUNT) + threadIdx.x] =dc_q[isp]*sum*0.5*ihx;
  }  
}
/* end find_current() */

//find charge density
__global__ void find_rho(Particles pcles, int isp, float_p *d_rho)
{  
  const int threadPos = threadIdx.x;
  __shared__ float_p s_rho[SMEMSIZE]; //fixed for the time being
  //initialize shared memory
  for(int xi=0; xi<RHO_BIN_COUNT; xi++ ){ //careful
    s_rho[threadIdx.x + xi*BLOCK_SIZE ] = 0.;
  }
  __syncthreads();

  for(int pi=dc_acNpcls[isp]+UMAD(blockIdx.x, blockDim.x, threadIdx.x); pi<dc_acNpcls[isp+1]; pi+=UMUL(blockDim.x,gridDim.x) ){
    float_p x_np  = pcles.d_x_np [pi];
    int   i_np  = pcles.d_ci_np[pi]; 
    //    float_p q     = dc_q  [isp];
#ifdef _OFFSET
# if _SPORDER%2==0         
    s_gather(s_rho+threadPos,dc_q[isp],i_np,x_np*ihx-0.5);  
# else
    s_gather(s_rho+threadPos,dc_q[isp],i_np,x_np*ihx);  
# endif
#else    
    s_gather(s_rho+threadPos,dc_q[isp],i_np,(x_np*ihx-i_np+1));
#endif
  }
  __syncthreads();

  //reduce
  if(threadIdx.x < RHO_BIN_COUNT){ //careful
    float_p *s_rhoBase = s_rho + UMUL(threadIdx.x, BLOCK_SIZE);
    float_p sum = 0;
    int pos;
    for(int j=0; j<G_SHARED_BANKS; j++){
      pos = threadIdx.x & (SHARED_MEMORY_BANKS-1);//to avoid bank conflict
      for(int i=0; i<SHARED_MEMORY_BANKS; i++){
	sum += s_rhoBase[j*SHARED_MEMORY_BANKS+pos];
	pos = (pos + 1) & (SHARED_MEMORY_BANKS-1);
      }
    }

    d_rho[(isp*gridDim.x+blockIdx.x)*(RHO_BIN_COUNT) + threadIdx.x] = sum;
  }


    // if(index==0){
    //   cuPrintf("\n\n");
    //   for(int i=0;i<Nx+2;i++){
    // 	cuPrintf("%f %f\n",xx[i],rho[i]);
    //   }
    // }

  //}
  
}
/* end find_rho() */



/* field on particle */
__device__ float_p fld_on_pcle(float_p *d_E, int i, float_p x)
{
  return FETCH(d_E,i)*(_fconst(1.0)-x) + FETCH(d_E,i+1)*x;
  //  return d_E[i]*(1.-x) + d_E[i+1]*x;
}
/* end fld_on_pcle() */



/* push velocity backward initially */
__global__ void init_vp_exp(Particles pcles, int *acNpcls,float_p *d_xx,float_p *d_q_m,float_p *d_E,float_p dt)
{
  for(int pi=UMAD(blockIdx.x, blockDim.x, threadIdx.x); pi<totNp; pi+=UMUL(blockDim.x,gridDim.x) ){
    //select species
    int isp = select_species_dc(pi);//select_species(pi,acNpcls);

  //save old positions and velocities
    float_p x_n  = pcles.d_x_np [pi];
    float_p v_n  = pcles.d_v_np [pi];
    int i_n  = pcles.d_ci_np[pi];
   
    float_p q_m=d_q_m[isp];
#ifdef _OFFSET
    float_p Ep = fld_on_pcle(d_E,i_n,x_n);
#else
    float_p Ep = fld_on_pcle(d_E,i_n,(x_n-FETCH(d_xx,i_n))*ihx);
#endif    
    pcles.d_v_np[pi] = v_n - 0.5*dt*q_m*Ep;
    // if(pi==62)
    //   cuPrintf("#in init_vp_exp: %d %e %e\n%d %e %e %e %e",pi,pcles.d_v_n[pi],pcles.d_v_np[pi],i_n,d_E[i_n],d_E[i_n+1],x_n,Ep);
  }
}
/* end init_vp_exp() */


//perturb particles
__global__ void perturb_pcles2(short nh1,float_p eps_pic, Particles pcles, int isp, float_p *d_xx,float_p dt,float_p *d_rho)
{
  if(totNp==2) return;
  const int threadPos = threadIdx.x;
  __shared__ float_p s_rho[SMEMSIZE]; //fixed for the time being

  //initialize shared memory
  for(int xi=0; xi<RHO_BIN_COUNT; xi++ ){
    s_rho[threadIdx.x + xi*BLOCK_SIZE ] = _fconst(0.0);
  }
  __syncthreads();

  float_p k0 = _fconst(2.0)*_acos(_fconst(-1.0))*float_p(nh1)/LL;
  float_p aa = eps_pic/k0;
  float_p dxp;
  for(int pi=dc_acNpcls[isp]+UMAD(blockIdx.x, blockDim.x, threadIdx.x); pi<dc_acNpcls[isp+1]; pi+=UMUL(blockDim.x,gridDim.x) ){
#ifdef _OFFSET

    //save old positions and velocities
    int i_n  = pcles.d_ci_n[pi];
# if _SPORDER%2==0         
    float_p xp = pcles.d_x_n[pi] - _fconst(0.5)*hx + d_xx[i_n] ;
# else
    float_p xp = pcles.d_x_n[pi] + d_xx[i_n];
# endif
    dxp = k0*xp;
    float_p x_np  = pcles.d_x_n [pi] + aa*_cos(dxp);
   
    //Locate particle in mesh (between 0 and nx)
    // int di = _floor(x_np*ihx);
    // int i_np = i_n + di;
    // x_np -= di*hx;
    int i_np = i_n;
    if(x_np<0){
      do{
	x_np += hx;
	i_np--;
      }while(x_np<0);
    }else if(x_np>hx){
      do{
	x_np -= hx;
	i_np++;
      }while(x_np>hx);
    }

    if(x_np==_fconst(0.0)){
      float_p v_np=pcles.d_v_np[pi];
      if(v_np<_fconst(0.0)){
	x_np = hx;
	i_np--;
      }
    }else if(x_np==hx){
      float_p v_np=pcles.d_v_np[pi];
      if(v_np>_fconst(0.0)){
        x_np = _fconst(0.0);
        i_np++;
      }
    }
    //B.C.
    //peroidic, assuming nxx power of 2  
# if _SPORDER%2==0         
    i_np = ((i_np-1) & (Nx-1)) + 1;
# else
    i_np = i_np & (Nx-1);
# endif    
    //write back to global memory
    pcles.d_x_n[pi] = pcles.d_x_np[pi] = x_np;
    pcles.d_ci_n[pi]= pcles.d_ci_np[pi]= i_np;
 
    //collect charge density
# if _SPORDER%2==0         
    s_gather(s_rho+threadPos,dc_q[isp],i_np,x_np*ihx-_fconst(0.5));  
# else
    s_gather(s_rho+threadPos,dc_q[isp],i_np,x_np*ihx);  
# endif
#else

    //save old positions and velocities
    // float_p x_np  = pcles.d_x_np [pi];    
    // x_np = x_np + eps/k0*cos(k0*x_np);
    float_p xp  = pcles.d_x_np [pi];    
    float_p x_np = xp + eps/k0*cos(k0*xp);
    //B.C.
    x_np = par_bc_pcles(x_np);
    int i_np = locate(x_np);

    //write back to global memory
    pcles.d_x_n[pi] = pcles.d_x_np[pi] = x_np;
    pcles.d_ci_n[pi]= pcles.d_ci_np[pi]= i_np;

    //collect charge density
    s_gather(s_rho+threadPos,dc_q[isp],i_np,(x_np*ihx-i_np+1));  
#endif

    //perturb velocity
    // float_p v_np  = pcles.d_v_np [pi];
    // v_np = v_np + eps*sin(k0*xp);
    // pcles.d_v_n[pi] = pcles.d_v_np[pi] = v_np;
      
  }

  //Accummulate per-thread rho into per-block and write to global memory
  __syncthreads();


  if(threadIdx.x < RHO_BIN_COUNT){ //careful
    float_p *s_rhoBase = s_rho + UMUL(threadIdx.x, BLOCK_SIZE);
    float_p sum = 0;
    int pos;
    for(int j=0; j<G_SHARED_BANKS; j++){
      pos = threadIdx.x & (SHARED_MEMORY_BANKS-1);//to avoid bank conflict
      for(int i=0; i<SHARED_MEMORY_BANKS; i++){
	sum += s_rhoBase[j*SHARED_MEMORY_BANKS+pos];
	pos = (pos + 1) & (SHARED_MEMORY_BANKS-1);
      }
    }

    d_rho[(isp*gridDim.x+blockIdx.x)*(RHO_BIN_COUNT) + threadIdx.x] = sum;
  }
    
}
/* end perturb_pcles2() */


/* explicit leap-frog mover */
__global__ void move_pcles2(Particles pcles, int isp, float_p *d_xx,float_p *d_q_m,float_p *d_E,float_p dt,float_p *d_q,float_p *d_rho)
{
  const int threadPos = threadIdx.x;
  __shared__ float_p s_rho[SMEMSIZE]; //fixed for the time being
  //initialize shared memory
  for(int xi=0; xi<RHO_BIN_COUNT; xi++ ){ //careful
    s_rho[threadIdx.x + xi*BLOCK_SIZE ] = 0.;
  }
  __syncthreads();

  for(int pi=dc_acNpcls[isp]+UMAD(blockIdx.x, blockDim.x, threadIdx.x); pi<dc_acNpcls[isp+1]; pi+=UMUL(blockDim.x,gridDim.x) ){

  //save old positions and velocities
   float_p x_n  = pcles.d_x_np [pi];
   float_p v_n  = pcles.d_v_np [pi];
   int   i_n  = pcles.d_ci_np[pi]; //locate(x_n); //

   pcles.d_x_n [pi]  = x_n; //careful
   pcles.d_v_n [pi]  = v_n;
   pcles.d_ci_n [pi] = i_n;

   float_p x_np, v_np;
   int   i_np;
   float_p q   = dc_q  [isp]; //d_q  [isp]; //-0.001259; //
   float_p q_m = dc_q_m[isp]; //d_q_m[isp]; //-1.; //
   float_p Ep   = fld_on_pcle(d_E,i_n,(x_n*ihx-(i_n-1)));

   //Find new velocities first
   v_np = v_n + dt*q_m*Ep;
   x_np = x_n + dt*v_np;
   //B.C.
   x_np = par_bc_pcles(x_np);
   i_np = locate(x_np);
  
   //write back to global memory  
   pcles.d_x_np[pi] = x_np;
   pcles.d_v_np[pi] = v_np;
   pcles.d_ci_np[pi]= i_np;

   s_gather(s_rho+threadPos,q,i_np,(x_np*ihx-i_np+1));  
  }
  __syncthreads();

  if(threadIdx.x < RHO_BIN_COUNT){ //careful
    float_p *s_rhoBase = s_rho + UMUL(threadIdx.x, BLOCK_SIZE);
    float_p sum = 0;
    int pos;
    for(int j=0; j<G_SHARED_BANKS; j++){
      pos = threadIdx.x & (SHARED_MEMORY_BANKS-1);//to avoid bank conflict
      for(int i=0; i<SHARED_MEMORY_BANKS; i++){
	sum += s_rhoBase[j*SHARED_MEMORY_BANKS+pos];
	pos = (pos + 1) & (SHARED_MEMORY_BANKS-1);
      }
    }

    d_rho[(isp*gridDim.x+blockIdx.x)*(RHO_BIN_COUNT) + threadIdx.x] = sum;
  }



}
/* end move_pcles2() */

__global__ void merge_d_rho(float_p * d_rho)
{
  float_p mySum=_fconst(0.0);
  float_p error=_fconst(0.0);
  float_p a,b;
  for(int i=threadIdx.x; i<(Nx)*Nsp*GRID_SIZE; i+=UMUL(blockDim.x,gridDim.x) ){
    //mySum +=  d_rho[i]; //*num[i];
    a = mySum;
    b = error+d_rho[i];
    mySum = a+b;
    error = b-(mySum-a);
  }
  __syncthreads(); //careful
  d_rho[threadIdx.x] = mySum;
  
}
/* end merge_d_rho() */


/*--- diagnostics ---*/
__global__ void output_mv(float_p *d_v_n, int isp, float_p *d_moments,int blockSize)
{
  extern __shared__ float_p shared[];

  int tid      = threadIdx.x;
  int i        = blockIdx.x*blockSize*2 + threadIdx.x + dc_acNpcls[isp];//careful
  int gridSize = blockSize*2*gridDim.x;
  
  float_p mySum   = _fconst(0.0);
  float_p error   = _fconst(0.0);
  float_p a,b;
  while (i < dc_acNpcls[isp+1]) {     //careful
    //mySum +=  dc_m[isp]*d_v_np[i]; 
    //mySum +=  d_v_n[i]; 
    a = mySum;
    b = error+d_v_n[i]; 
    mySum = a+b;
    error = b-(mySum-a);
    if (i + blockSize < dc_acNpcls[isp+1]){ //careful
      //mySum +=  dc_m[isp]*d_v_np[i+blockSize];
      //mySum +=  d_v_n[i+blockSize];
      a = mySum;
      b = error+d_v_n[i+blockSize]; 
      mySum = a+b;
      error = b-(mySum-a);
    }
    i += gridSize;               
  }

  shared[tid] = mySum;
  __syncthreads();

  //do reduction in shared memory
  for(int s=blockDim.x/2; s>32; s>>=1) {
      if(tid<s){
	shared[tid] = mySum = mySum+shared[tid+s];
      }
      __syncthreads();
    }

  if (tid < 32) {
         float_p *smem = shared;
         if (blockSize >=  64) { smem[tid] = mySum = mySum + smem[tid + 32]; __syncthreads(); }
         if (blockSize >=  32) { smem[tid] = mySum = mySum + smem[tid + 16]; __syncthreads(); }
         if (blockSize >=  16) { smem[tid] = mySum = mySum + smem[tid +  8]; __syncthreads(); }
         if (blockSize >=   8) { smem[tid] = mySum = mySum + smem[tid +  4]; __syncthreads(); }
         if (blockSize >=   4) { smem[tid] = mySum = mySum + smem[tid +  2]; __syncthreads(); }
         if (blockSize >=   2) { smem[tid] = mySum = mySum + smem[tid +  1]; __syncthreads(); }
    
  }
      
  if(tid == 0) {
        d_moments[isp*gridDim.x + blockIdx.x] = dc_m[isp]*shared[0];  	
  }
  
}
/* output_mv() */

__global__ void output_v2(float_p *d_v_n,float_p *d_v_np, int isp, float_p *d_moments,int blockSize)
{
  extern __shared__ float_p shared[];

  int tid      = threadIdx.x;
  int i        = blockIdx.x*blockSize*2 + threadIdx.x + dc_acNpcls[isp];//careful
  int gridSize = blockSize*2*gridDim.x;
  
  float_p mySum   = 0.;
  while (i < dc_acNpcls[isp+1]) {     //careful
    mySum +=  d_v_n[i]*d_v_np[i]; 
    if (i + blockSize < dc_acNpcls[isp+1])  //careful
      mySum +=  d_v_n[i+blockSize]*d_v_np[i+blockSize];
    i += gridSize;               
  }

  shared[tid] = mySum;
  __syncthreads();

  //do reduction in shared memory
  for(int s=blockDim.x/2; s>32; s>>=1) {
      if(tid<s){
	shared[tid] = mySum = mySum+shared[tid+s];
      }
      __syncthreads();
    }

  if (tid < 32) {
         float_p *smem = shared;
         if (blockSize >=  64) { smem[tid] = mySum = mySum + smem[tid + 32]; __syncthreads(); }
         if (blockSize >=  32) { smem[tid] = mySum = mySum + smem[tid + 16]; __syncthreads(); }
         if (blockSize >=  16) { smem[tid] = mySum = mySum + smem[tid +  8]; __syncthreads(); }
         if (blockSize >=   8) { smem[tid] = mySum = mySum + smem[tid +  4]; __syncthreads(); }
         if (blockSize >=   4) { smem[tid] = mySum = mySum + smem[tid +  2]; __syncthreads(); }
         if (blockSize >=   2) { smem[tid] = mySum = mySum + smem[tid +  1]; __syncthreads(); }
    
  }
      
  if(tid == 0) {
        d_moments[isp*gridDim.x + blockIdx.x] = shared[0];  	
	//cuPrintf("%f\n",shared[0]);	
  }
  
}
/* output_v2() */


/*-----offset functions-----*/
// //perturb particles_offset
// __global__ void perturb_pcles_offset(short nh1,float_p eps, Particles pcles, int isp, float_p *d_xx,float_p dt,float_p *d_rho)
// {
//   if(totNp==2) return;
//   const int threadPos = threadIdx.x;
//   __shared__ float_p s_rho[SMEMSIZE]; //fixed for the time being

//   //initialize shared memory
//   for(int xi=0; xi<RHO_BIN_COUNT; xi++ ){
//     s_rho[threadIdx.x + xi*BLOCK_SIZE ] = 0.;
//   }
//   __syncthreads();

//   float_p k0 = 2*PI*nh1/LL;
//   for(int pi=dc_acNpcls[isp]+UMAD(blockIdx.x, blockDim.x, threadIdx.x); pi<dc_acNpcls[isp+1]; pi+=UMUL(blockDim.x,gridDim.x) ){

//     //save old positions and velocities
//     int i_n  = pcles.d_ci_n[pi];
//     float_p xp = pcles.d_x_n[pi] + d_xx[i_n];

//     float_p x_np  = pcles.d_x_np [pi] + eps/k0*cos(k0*xp);
   
//     //Locate particle in mesh (between 0 and nx)
//     int di = floor(x_np*ihx);
//     int i_np = i_n + di;
//     x_np -= di*hx;
//     //B.C.
//     //peroidic, assuming nxx power of 2  
// #if _SPORDER%2==0         
//     i_np = ((i_np-1) & (Nx-1)) + 1;
// #else
//     i_np = i_np & (Nx-1);
// #endif    
//     //write back to global memory
//     pcles.d_x_n[pi] = pcles.d_x_np[pi] = x_np;
//     pcles.d_ci_n[pi]= pcles.d_ci_np[pi]= i_np;
 
//     //collect charge density
//     s_gather(s_rho+threadPos,dc_q[isp],i_np,x_np*ihx);  


//   }

//   //Accummulate per-thread rho into per-block and write to global memory
//   __syncthreads();


//   if(threadIdx.x < RHO_BIN_COUNT){ //careful
//     float_p *s_rhoBase = s_rho + UMUL(threadIdx.x, BLOCK_SIZE);
//     float_p sum = 0;
//     int pos;
//     for(int j=0; j<G_SHARED_BANKS; j++){
//       pos = threadIdx.x & (SHARED_MEMORY_BANKS-1);//to avoid bank conflict
//       for(int i=0; i<SHARED_MEMORY_BANKS; i++){
// 	sum += s_rhoBase[j*SHARED_MEMORY_BANKS+pos];
// 	pos = (pos + 1) & (SHARED_MEMORY_BANKS-1);
//       }
//     }
//     d_rho[(isp*gridDim.x+blockIdx.x)*(RHO_BIN_COUNT) + threadIdx.x] = sum;
//   }
    
// }
// /* end perturb_pcles_offset() */



/* explicit leap-frog mover using offset */
__global__ void move_pcles_offset(Particles pcles, int isp, float_p *d_xx,float_p *d_E,float_p dt,float_p *d_rho)
{
  const int threadPos = threadIdx.x;
  __shared__ float_p s_rho[SMEMSIZE]; //fixed for the time being
  //initialize shared memory
  for(int xi=0; xi<RHO_BIN_COUNT; xi++ ){ //careful
    s_rho[threadIdx.x + xi*BLOCK_SIZE ] = 0.;
  }
  __syncthreads();

  for(int pi=dc_acNpcls[isp]+UMAD(blockIdx.x, blockDim.x, threadIdx.x); pi<dc_acNpcls[isp+1]; pi+=UMUL(blockDim.x,gridDim.x) ){
  //save old positions and velocities    
   float_p x_n  = pcles.d_x_np [pi];
   float_p v_n  = pcles.d_v_np [pi];
   int i_n    = pcles.d_ci_np[pi]; 

   pcles.d_x_n [pi]  = x_n; //careful
   pcles.d_v_n [pi]  = v_n;
   pcles.d_ci_n [pi] = i_n;


   float_p x_np, v_np;

   float_p q   = dc_q  [isp]; //d_q  [isp]; //-0.001259; //
   float_p q_m = dc_q_m[isp]; //d_q_m[isp]; //-1.; //
   float_p Ep   = fld_on_pcle(d_E,i_n,(x_n*ihx));

   //Find new velocities first
   v_np = v_n + dt*q_m*Ep;
   x_np = x_n + dt*v_np;

    int di = floor(x_np*ihx);
    int i_np = i_n + di;
    x_np -= di*hx;
    //B.C.
    //peroidic, assuming nxx power of 2           
    i_np = i_np & (Nx-1);
      
   //write back to global memory  
   pcles.d_x_np[pi] = x_np;
   pcles.d_v_np[pi] = v_np;
   pcles.d_ci_np[pi]= i_np;
   
   s_gather(s_rho+threadPos,dc_q[isp],i_np,x_np*ihx);  

  }
  __syncthreads();


  if(threadIdx.x < RHO_BIN_COUNT){ //careful
    float_p *s_rhoBase = s_rho + UMUL(threadIdx.x, BLOCK_SIZE);
    float_p sum = 0;
    int pos;
    for(int j=0; j<G_SHARED_BANKS; j++){
      pos = threadIdx.x & (SHARED_MEMORY_BANKS-1);//to avoid bank conflict
      for(int i=0; i<SHARED_MEMORY_BANKS; i++){
	sum += s_rhoBase[j*SHARED_MEMORY_BANKS+pos];
	// if(threadIdx.x==Nx){
	//   cuPrintf("%d %d %f\n",blockIdx.x,pos,sum);
	// }
	pos = (pos + 1) & (SHARED_MEMORY_BANKS-1);
      }
    }

    //d_rho[(blockIdx.x*Nsp+isp)*(RHO_BIN_COUNT+1) + threadIdx.x] = sum;
    d_rho[(isp*gridDim.x+blockIdx.x)*(RHO_BIN_COUNT) + threadIdx.x] = sum;
  }

  // __syncthreads();


}
/* end move_pcles_offset() */


/*------ implicit functions -----*/

/* field on particle (at half grid points)*/
__device__ float_p fld_on_pcleh(float_p *d_E, int i, float_p x)
{
  //#if   _SPORDER == 2  
  return FETCH(d_E,i-1)*(_fconst(1.0)-x) + FETCH(d_E,i)*x;
  //#endif
}
/* end fld_on_pcleh() */

//particle move for a sub-time step
__device__ void acc_cn_incell(int isp, float_p x_n,float_p v_n,int i_n, float_p dt, float_p *d_E, float_p &x_np, float_p &v_np, int &i_np, float_p &x_nh, int &i_nh,float_p &dt_sub)
{
  dt_sub = dt;
  //Source terms from last time step
  float_p a_n  = dc_q_m[isp]*fld_on_pcleh(d_E,i_n,(x_n*ihx-i_n+1.5));
  if(v_n==0 && a_n==0.) return;
  //estimate sub-time-step 
  float_p dadx = dc_q_m[isp]*(FETCH(d_E,i_n)-FETCH(d_E,i_n-1))*ihx;
  float_p aa = sqrt(dadx*dadx*v_n*v_n+a_n*a_n);
  float_p bb = etolr * sqrt(v_n*v_n+a_n*a_n);
  float_p cc = etola;
  
  //dt_sub = ( bb + sqrt(bb*bb+2.0*aa*cc) )/ (aa+eps);
  
  //dt_sub = fmin(dt_sub,dt);

  
  //start with dt_sub
  //converge
  v_np = v_n;
  // Initialize nonlinear iteration
  float_p dx_pn = fabs(dt*a_n);
  float_p dx_pk = dx_pn;
  ushort   it = 0;
  float_p x_k, dx;
  float_p am;

  while(dx_pk>eps*dx_pn && it<200){
  //while(it<5){
    x_k  = v_np;
    dx   = (v_np+v_n)*dt_sub*0.5; //careful
    x_np = x_n + dx*0.5;
    //B.C.
    x_np = par_bc_pcles(x_np);
    i_np = locate(x_np);
    //am = dc_q_m[isp]*fld_on_pcleh(d_E,i_np,(x_np*ihx-i_np+1.5));
    am   = a_n + dadx*dx*0.5;
    v_np = v_n + am*dt_sub; //careful
    
    //Check convergence
    dx_pk = fabs(v_np - x_k);
    it++;
  }
  i_nh = i_np;
  x_nh = x_np;
  x_np = x_n + dx;
  //B.C.
  x_np = par_bc_pcles(x_np);
  i_np = locate(x_np);
 

  int cross = (fabs(dx)>hx) || (i_np != i_n);
  if(cross){
   //  i_npc = i_n + sign_int(dx);
     int di = (dx<0) ? -1 : 1 ;
     i_np = i_n + di ;


#if   _SPORDER == 2
     //x_np = hx*(i_n-1+i_np-1)*0.5; //careful   
     //x_np = hx*(i_np-0.5) - ((1+di)/2)*hx;   
     x_np = hx*(i_np- (1+di)/2) -0.5*hx;   
     // if(x_np != x_np2){
     //   cuPrintf("%e %e %e\n",x_np,x_np2,x_np-x_np2);
     // }
#endif  
     dx = x_np - x_n;
     x_nh = x_n + dx*0.5;
   //B.C
   //x_nh = par_bc_pcles(x_nh);
     i_nh = i_n;

   //am = dc_q_m[isp]*fld_on_pcleh(d_E,i_nhc,(x_nhc*ihx-i_nhc+1.5));
     am   = a_n + dadx*dx*0.5;
     v_np = sign(dx)*sqrt(fabs(2.*am*dx+v_n*v_n));
     
     dt_sub = (dx+v_np-v_n)/(0.5*(v_np+v_n)+am);
#if   _SPORDER == 2
     i_np = ((i_np-1) & (Nx-1)) +1;
     //x_np = hx*(i_np-0.5) - ((1+di)/2)*hx;   
     x_np = hx*(i_np- (1+di)/2) -0.5*hx;   
      // if(i_np == 33){
      //   i_np = 1;
      //   x_np = -0.5*hx;
      // }else if(i_np == 0){
      //   i_np = 32;
      //   x_np = LL - 0.5*hx;
      // }
#endif
   }
  
}
/* end acc_cn_incell() */

//particle move for a sub-time step
__device__ void acc_cn_incell_offset(int isp, float_p x_n,float_p v_n,int i_n, float_p dt, float_p *d_E, float_p &x_np, float_p &v_np, int &i_np, float_p &x_nh, int &i_nh,float_p &dt_sub,int pi)
{
  dt_sub = dt;
  //Source terms from last time step
  //float_p a_n  = dc_q_m[isp]*fld_on_pcleh(d_E,i_n,(x_n*ihx));
  float_p E0 = FETCH(d_E,i_n-1);
  float_p E1 = FETCH(d_E,i_n  );
  //  float_p x = x_n*ihx;                           //1

  float_p N = x_n;
  float_p D = hx;
  float_p R = ihx;
  float_p x = N*R;                           //1
  float_p r = N - x*D;
  x += r*R;

  float_p a_n;
  //a_n = _fma(-x,E0,E0);                          //2(1)
  //a_n = _fma(x,E1,a_n);                          //2(1)
  a_n = E0 - E0*x + E1*x;
  a_n*= dc_q_m[isp];
  float_p a_nd = a_n*_fconst(2.0);
  float_p v_nd = v_n*_fconst(2.0);
  //if(v_n==0 && a_n==0.) return;
  //estimate sub-time-step 
  //float_p dadx = dc_q_m[isp]*(FETCH(d_E,i_n)-FETCH(d_E,i_n-1))*ihx;
  //float_p dadx = dc_q_m[isp]*ihx*(E1-E0);
  
  N = (E1-E0);
  R = ihx;
  D = hx;
  float_p dadx = N*R;
  r = N - dadx*D;
  dadx += r*R;
  dadx *= dc_q_m[isp];

  float_p  aa = _abs(a_n)+_abs(dadx*v_n);
  float_p  bb = etolr*(_abs(v_n) + _abs(a_n));
  float_p  cc = etola;

  //float_p dt_sube = _max(cc/_sqrt(aa), bb/aa*_fconst(2.0));

  
#ifdef _DP   
  bb = bb/aa;
#else
  //bb = bb/aa;
  //  bb = __fdividef(bb,aa);

  //aa = _rsqrt(aa+eps);
  aa = _rsqrt(aa);
  dt_sub = cc*(aa);
  bb = bb*aa*aa;

#endif
  float_p dt_sube = _max(cc, bb*_fconst(2.0));


  /*
  float_p aa = __fsqrt_ru(dadx*dadx*v_n*v_n+a_n*a_n);  //4(1)
  float_p bb = etolr * __fsqrt_ru(v_n*v_n+a_n*a_n);    //3(1)
  float_p cc = etola;
  
  float_p dt_sube = ( bb + _sqrt(bb*bb+_fconst(2.0)*aa*cc) )/ (aa+eps);  //6(1)
  */
  dt_sub = _min(dt_sube,dt);
  //dt_sub = _min(dt_sube,_fconst(0.5)*hx/v_n);

  //if(pi==167736){
  //if(dt_sub==0){
  //   printf("%d %e %e %e %e %e %e %e %e %e\n",pi,dt, dt_sub,dt_sube,a_n,v_n,dadx,aa,bb,cc);
  // }
  
  //single step cn move                                                     
  float_p dth_sub = _fconst(0.5)*dt_sub;            //1
  float_p dth_sub2 = dth_sub*dth_sub;
  D = _fconst(1.0)-dadx*dth_sub2;
  //D = _fconst(1.0)-dadx*dth_sub*dth_sub;
  //N = v_n + dadx*dth_sub*dth_sub*v_n + a_n*dt_sub;
  N = dadx*dth_sub*v_n + a_n;
  //float_p N = _fma(a_nd,dth_sub,v_nd)*dth_sub;

#ifdef _DP 
  float_p dvbydt = _div(N,D);
  v_np = v_n + dvbydt*dt_sub;

  //v_np = _div(N,D);
  /*
  float Nf = float(N);
  float Df = float(D);
  float Rf;
  asm("rcp.approx.f32 %0, %1;\n\t"   //a = rcp(c) approximate       
      : "=f"(Rf) :  "f"(Df));
  float dxf = Nf*Rf;
  float r = Nf - dxf*Df;
  float_p dx = dxf + r*Rf;
  */

#else
  float_p dvbydt = _div(N,D);
  v_np = v_n + dvbydt*dt_sub;

  //v_np = N/D;
  
  /*  
  asm("rcp.approx.f32 %0, %1;\n\t"   //a = rcp(c) approximate       
  : "=f"(R) :  "f"(D));
  //R= _rcp(D);

  v_np = N*R;
  r = N - v_np*D; //Newton-Raphson     
  v_np = v_np + r*R;      //necessary
  */
  /*
  float_p dx = N*R;
  float_p r = N - dx*D; //Newton-Raphson
  dx = dx + r*R;
  */
#endif


  //float_p am;                             //5(1)
  //float_p dx = dth_sub*(v_n + v_np);                 //2(1)
  //float_p fd = _fma(dadx,dx,a_nd);
  //v_np = _fma(fd,dth_sub,v_n);


  //x_np = x_n + dx; 
  x_np = x_n + dth_sub*v_n + dth_sub*v_np; 
  float_p dx = x_np - x_n;
  i_np = i_n;

  //float_p amd = a_n + E0*dc_q_m[isp] + dadx*x_np;
  float_p amd = a_n + E0*dc_q_m[isp] + dadx*x_np;

  //int di = _floor(x_np*ihx);
  // N = x_np;
  // R = ihx;
  // D = hx;
  // x = N*R;
  // r = N - x*D;
  // x+= r*R;
  // int di = _floor(x);

  //if(pi==167736){
    //printf("%d %e %e %e %e %e %e %e %e %e\n",pi,dt, dt_sub,dt_sube,a_n,v_n,dadx,aa,bb,cc);
    //printf("%d %e  %e %d %e %e %e %e %e %e %e %e",pi,N,D,di,x_n,x_np,dx,v_n,v_np,a_n,dadx,dt_sub);
    //printf("%24.12e %24.12e %24.12e %24.12e %24.12e %24.12e",v_n,v_np,x_n,x_np,dt_sub,dx);
  //}

  /*
  //test without crossing  i_np += di; //(di<0)?(-1):(1); //careful
  i_np = ((i_np-1) & (Nx-1)) +1;
  x_np = (di<0 ? hx : _fconst(0.0));
  */
  float_p vmd,sign_dx;

  /*
   if( (x_np==hx)&&(v_np==0) ){
   }else if( ((x_np==0)&&(v_np>0)) || ((x_np==hx)&&(v_np<0)) ){
   }else if( (di!=0) ){
  */
  //if( (x_np<0) || (x_np>hx) ){
  //bool cross = (x_np<0) || (x_np>hx) || ( (x_np==0)&&(v_np>0) )|| ( (x_np==hx)&&(v_np<0) ) ;
  bool cross = (x_np<0) || (x_np>hx)  ;
  dt_sube = dt_sub;
    if( cross){

     //i_np += sign(di);                        //,2

     if(x_np<=0){
        dx = -x_n;
        x_np = hx;
	i_np--;
        sign_dx = _fconst(-1.0);
	amd = a_n+E0*dc_q_m[isp];
     }else{
        dx = hx-x_n;
        x_np = _fconst(0.0);
	i_np++;
        sign_dx = _fconst(1.0);
	amd = a_n+E1*dc_q_m[isp];
      }
     i_np = ((i_np-1) & (Nx-1)) +1;           //,4
     //x_np = (di>0 ? _fconst(0.0) : hx);      //,1    
     //float_p sign_dx = sign(dx);
     //dx = (di<=0 ?  _fconst(0.0) : hx) - x_n;  //1,1

     //am   = a_n + dadx*_fconst(0.5)*dx;      //3(1)
     //amd   = _fma(dadx,dx,a_nd);      //3(1)
     
     float_p vn2 = v_n*v_n;
     //float_p vnp2 = _fma(amd,dx,vn2) + _fma(v_n,v_n,-vn2);
     float_p vnp2 = vn2 + amd*dx;
     float_p vnp2_2 = vn2 - v_n*v_n;
     vnp2 -= vnp2_2;

     //float_p vnp2 = amd*dx+v_n*v_n;
     if(vnp2<=_fconst(0.0)){ //careful
       //vnp2 = _fconst(0.0);
       //}
       //reset
       v_np = _fconst(0.0); 
       i_np = i_n;
       if(x_np==_fconst(0.0)){
	 x_np = hx;
       }else{
	 x_np = _fconst(0.0);
       }
       // vmd = (v_n+v_np)*_fconst(0.5);
       // if(vmd==_fconst(0.0)){                                   //1
       // 	 dt_sub = _div((v_np-v_n),am);
       // }else{
       // 	 dt_sub = _div(dx,vmd);       
       // } 
     }else{

 #ifdef _DP 
     v_np = sign_dx*_sqrt(vnp2);
     vmd = (v_n+v_np);
     //if(vmd==_fconst(0.0) || dx==0){                                   //1
     if(vmd==_fconst(0.0)){                                   //1
     //if(am!=_fconst(0.0)){
       dt_sub = _div((v_np-v_n),amd)*_fconst(2.0);
     }else{
       dt_sub = _div(dx,vmd)*_fconst(2.0);     
     } 

 #else
     v_np = sign_dx*_sqrt(vnp2);
     vmd = (v_n+v_np);
     //if(vmd==_fconst(0.0) || dx==0){                                   //1
     if(vmd==_fconst(0.0)){                                   //1
     //if(am!=_fconst(0.0)){
       dt_sub = _div((v_np-v_n),amd)*_fconst(2.0);
     }else{
       dt_sub = _div(dx,vmd)*_fconst(2.0);     
     } 

     /*
     float_p ivnp = _rsqrt(vnp2);
     dth_sub = __fdividef(-v_n*ivnp + sign(dx), amd*ivnp);
     v_np = _fma(amd,dth_sub,v_n);
     float_p F =  _fma(dth_sub,v_np, -dx);    //4(1)
     F = _fma(dth_sub,v_n,F);
     dth_sub -= __fdividef(F, v_np)*_fconst(0.5);
     v_np = _fma(amd,dth_sub,v_n);

     dt_sub = _fconst(2.0)*dth_sub;
     */
 #endif

     }

   }else if( (x_np==0)&&(v_np<0) ){
     i_np--; 
     i_np = ((i_np-1) & (Nx-1)) +1; 
     x_np = hx;
   }else if( (x_np==hx)&&(v_np>0) ){
     i_np++;
     i_np = ((i_np-1) & (Nx-1)) +1;
     x_np = 0;
  }

  i_nh = i_n;
  x_nh = x_n + dx*_fconst(0.5);
  //am = a_n + dadx*dx*_fconst(0.5);
  dt_sub = _min(dt_sub, dt_sube);
  dth_sub = dt_sub*_fconst(0.5);
  v_np = v_n + amd*dth_sub;

    //     printf(" %e %e %e\n",x_n,x_np,dx);
  if(v_np!=v_np)
    printf("%d %24.12e %24.12e %24.12e %24.12e %24.12e %d %d\n",pi,v_n,v_np,x_n, x_np, dt_sub, i_n,i_np);
  //}

}
/* end acc_cn_incell_offset() */

//private domain per thread
__global__ void move_pcles_acc_cn2(Particles pcles, int isp, float_p *d_E,float_p dt,float_p *d_j)
{

  const int threadPos = threadIdx.x;
  __shared__ float_p s_j[SMEMSIZE]; //fixed for the time being
  __shared__ float_p s_e[SMEMSIZE];
  //initialize shared memory
  for(int xi=0; xi<RHO_BIN_COUNT; xi++ ){ //careful
    s_j[threadIdx.x + xi*BLOCK_SIZE ] = _fconst(0.0);
    s_e[threadIdx.x + xi*BLOCK_SIZE ] = _fconst(0.0);

  }
  __syncthreads();

  float_p *t_j = s_j + threadPos;
  float_p *t_e = s_e + threadPos;

  for(int pi=dc_acNpcls[isp]+UMAD(blockIdx.x, blockDim.x, threadIdx.x); pi<dc_acNpcls[isp+1]; pi+=UMUL(blockDim.x,gridDim.x) ){
     float_p dtp = dt;
     float_p dt_sub;
     //save old positions and velocities
     float_p x_n  = pcles.d_x_n [pi];
     float_p v_n  = pcles.d_v_n [pi];
     int   i_n  = pcles.d_ci_n[pi]; 

     float_p x_np,v_np,x_nh;
     int i_np,i_nh;
     float_p vh;
     float_p j0=_fconst(0.0),j1=_fconst(0.0);
     float_p e0=_fconst(0.0),e1=_fconst(0.0);
     float_p e00=_fconst(0.0),e11=_fconst(0.0);
     float_p a,b,c;
     int it = 0;
     //while(1){
     while(dtp>0){
       acc_cn_incell_offset(isp,x_n,v_n,i_n,dtp,d_E,x_np,v_np,i_np,x_nh,i_nh,dt_sub,pi);
       vh = (v_n+v_np)*_fconst(0.5);
       
       float_p q = vh*dt_sub;                     //1                  
       //float_p x = x_nh*ihx; 
       float_p N = x_nh;
       float_p R = ihx;
       float_p D = hx;
       float_p x = N*R; //x_nh*ihx;
       float_p r = N - x*D;
       x += r*R;



        int ie0 = (i_nh-1)*BLOCK_SIZE;
        int ie1 = ((i_nh)&(Nx-1))*BLOCK_SIZE;
        int i0 = (i_nh-1)*BLOCK_SIZE;
        int i1 = ((i_nh)&(Nx-1))*BLOCK_SIZE;

        e0 = t_e[ie0];
        e1 = t_e[ie1];
        j0 = t_j[i0];
        j1 = t_j[i1];

	c = q-q*x;
	if(_abs(c)>_abs(j0)){
	  a = c;
	  c = j0;
	}else{
	  a = j0;
	}
        b = e0 + c;
        j0 = a+b;
        e0 = b-(j0-a);

	c = q*x;
	if(_abs(c)>_abs(j1)){
	  a = c;
	  c = j1;
	}else{
	  a = j1;
	}
        b = e1+c;
        j1 = a+b;
        e1 = b-(j1-a);

        // j0 +=  q-q*x;
        // j1 += q*x;

        t_e[ie0] = e0;
        t_e[ie1] = e1;
        t_j[i0] = j0;
        t_j[i1] = j1;

       
       // j0 = q-q*x;   
       // j1 = q*x;

       
       // t_j[       (i_nh-1)*BLOCK_SIZE] += j0;
       // t_j[((i_nh)&(Nx-1))*BLOCK_SIZE] += j1;

       x_n = x_np;
       v_n = v_np;
       i_n = i_np;
       dtp -= dt_sub;
     }
     //write back to global memory  
     pcles.d_x_np[pi] = x_np;
     pcles.d_v_np[pi] = v_np;
     pcles.d_ci_np[pi]= i_np;
     
  }

  __syncthreads();
  //reduce
   if(threadIdx.x < RHO_BIN_COUNT){ //careful
     float_p *s_rhoBase = s_j + UMUL(threadIdx.x, BLOCK_SIZE);
     float_p *s_eBase = s_e + UMUL(threadIdx.x, BLOCK_SIZE);
     float_p sum = _fconst(0.0);
     float_p error=_fconst(0.0);
     float_p a,b,c;
     int pos;
     for(int j=0; j<G_SHARED_BANKS; j++){
       pos = threadIdx.x & (SHARED_MEMORY_BANKS-1);//to avoid bank conflict
       for(int i=0; i<SHARED_MEMORY_BANKS; i++){
	 //sum += s_rhoBase[j*SHARED_MEMORY_BANKS+pos];
	 //error += s_eBase[j*SHARED_MEMORY_BANKS+pos];
	 c = s_rhoBase[j*SHARED_MEMORY_BANKS+pos];
	 if(_abs(c)>_abs(sum)){
	   a = c;
	   c = sum;
	 }else{
	   a = sum;
	 }
	 b = error+s_eBase[j*SHARED_MEMORY_BANKS+pos]+c;
	 sum = a+b;
	 error = b-(sum-a);
	 
	 pos = (pos + 1) & (SHARED_MEMORY_BANKS-1);
	// if(threadIdx.x==31)
	//    cuPrintf("%d %d %d %d %f\n",threadIdx.x,j,i,pos,sum);
       }
     }

     d_j[(isp*gridDim.x+blockIdx.x)*(RHO_BIN_COUNT) + threadIdx.x] = sum;

   }  

}
/* end  move_pcles_acc_cn2() */

/* Crank-Nicolson mover (adaptive-charge-conserving)*/
__global__ void move_pcles_acc_cn(Particles pcles, int isp, float_p *d_E,float_p dt,float_p *d_j)
 {
   const float_p idt = _rcp(dt);

   //initialize shared memory
   extern __shared__ float_p s_j[];
   __shared__ float_p s_m[G_WARP*32];

   if(threadIdx.x<Nx*G_WARP){
     s_j[threadIdx.x] = _fconst(0.0);
     s_m[threadIdx.x] = _fconst(0.0);
   }
   __syncthreads();
   float_p *w_j = s_j + (threadIdx.x>>LOG_Nx)*Nx;
   float_p *w_m = s_m + (threadIdx.x>>LOG_Nx)*Nx;
   //float *w_j = s_j + (threadIdx.x>>LOG_WARP)*WARP_SIZE;

   for(int pi=dc_acNpcls[isp]+UMAD(blockIdx.x, blockDim.x, threadIdx.x); pi<dc_acNpcls[isp+1]; pi+=UMUL(blockDim.x,gridDim.x) ){
     float_p dtp = dt;
     float_p dt_sub;
     //save old positions and velocities
     float_p x_n  = pcles.d_x_n [pi];
     float_p v_n  = pcles.d_v_n [pi];
     int   i_n  = pcles.d_ci_n[pi]; 

     float_p x_np,v_np,x_nh;
     int i_np,i_nh;
     float_p vh;
     float_p j0=_fconst(0.0),j1=_fconst(0.0);
     float_p e0=_fconst(0.0),e1=_fconst(0.0);
     int it = 0;
     //while(1){
     while(dtp>0){
       acc_cn_incell_offset(isp,x_n,v_n,i_n,dtp,d_E,x_np,v_np,i_np,x_nh,i_nh,dt_sub,pi);
       vh = (v_n+v_np)*_fconst(0.5);
       
       float_p q = vh*dt_sub;                     //1                  
       //float_p x = x_nh*ihx; 
       float_p N = x_nh;
       float_p R = ihx;
       float_p D = hx;
       float_p x = N*R; //x_nh*ihx;
       float_p r = N - x*D;
       x += r*R;
       /*
       float_p qx,a,b,c;

	c = q-q*x;
	if(_abs(c)>_abs(j0)){
	  a = c;
	  c = j0;
	}else{
	  a = j0;
	}
        b = e0 + c;
        j0 = a+b;
        e0 = b-(j0-a);

	c = q*x;
	if(_abs(c)>_abs(j1)){
	  a = c;
	  c = j1;
	}else{
	  a = j1;
	}
        b = e1+c;
        j1 = a+b;
        e1 = b-(j1-a);
       */
       j0 += (q-q*x);   
       j1 += q*x;

       if(i_np!=i_n){      
#ifdef _DP
	 if(j0>=0){
	   atomicAdd2(w_j+i_nh-1       , j0);
	 }else{
	   atomicAdd2(w_m+i_nh-1       , j0);
	 }
	 if(j1>=0){	 
	   atomicAdd2(w_j+(i_nh&(Nx-1)), j1);
	 }else{
	   atomicAdd2(w_m+(i_nh&(Nx-1)), j1);
	 }
       //atomicAdd(w_j+i_nh-1       , float(j0));
       //atomicAdd(w_j+(i_nh&(Nx-1)), float(j1));
#else
       atomicAdd(w_j+i_nh-1       , j0);   
       atomicAdd(w_j+(i_nh&(Nx-1)), j1);  
#endif  
       j0 = _fconst(0.0); 
       j1 = _fconst(0.0);
       e0 = _fconst(0.0);
       e1 = _fconst(0.0);
       }
       //cuPrintf("[%d %d] %f %f %d, %f %f %d %f\n",isp,pi,x_n,v_n,i_n,x_np,v_np,pcles.d_ci_np[pi],dtp);
       // if(pi==79736){
       // 	 printf("%d %e  %e %d %d %e %e\n",pi,v_n, v_np, i_n,i_np, d_E[i_n],d_E[i_np]);
       // }
       it++;
       //if(it>390){
       //if(pi==167736){
       //printf("%d %d %e %e  %e %e %d %e %e %d\n",it, pi,dtp,dt_sub,x_n,v_n,i_n,x_np,v_np,i_np);
       if(it>400) break;
       //}

       x_n = x_np;
       v_n = v_np;
       i_n = i_np;
       dtp -= dt_sub;

       //if(pi==9416 && it>=400){
       // if(dtp>4.0f){
       //  	 printf("%d %e %e\n",pi,dtp,dt_sub);
       // }
       //  if(__all(dtp<=_fconst(0.0) || it>400))break;
       //if(__all(dtp<=_fconst(0.0)))break;

     }
  
    //write back to global memory  
     pcles.d_x_np[pi] = x_np;
     pcles.d_v_np[pi] = v_np;
     pcles.d_ci_np[pi]= i_np;

 #ifdef _DP
     if(j0>=0){
        atomicAdd2(w_j+i_nh-1       , j0);
     }else{
        atomicAdd2(w_m+i_nh-1       , j0);
     }
     if(j1>=0){
        atomicAdd2(w_j+(i_nh&(Nx-1)), j1);	 
     }else{
        atomicAdd2(w_m+(i_nh&(Nx-1)), j1);	 
     }
 #else
        atomicAdd(w_j+i_nh-1       , j0);
        atomicAdd(w_j+(i_nh&(Nx-1)), j1);
 #endif

   }
   __syncthreads();
   //reduce
   if(threadIdx.x < RHO_BIN_COUNT){ //careful

     float_p sum = _fconst(0.0);
     float_p mum = _fconst(0.0);
     float_p error=_fconst(0.0);
     float_p mrror=_fconst(0.0);
     float_p a,b,c;
     for(int i=0; i<G_WARP; i++){
     //for(int i=G_WARP-1; i>=0; i--){
       sum +=  s_j[threadIdx.x+i*Nx];
       mum +=  s_m[threadIdx.x+i*Nx];
       
       /*
       c = s_j[threadIdx.x+i*Nx]; 
       if(_abs(c)>_abs(sum)){
	 a = c;
	 c = sum;
       }else{
	 a = sum;
       }
       b = error+c;
       sum = a+b;
       error = b-(sum-a);

       c = s_m[threadIdx.x+i*Nx]; 
       if(_abs(c)>_abs(mum)){
	 a = c;
	 c = mum;
       }else{
	 a = mum;
       }
       b = mrror+c;
       mum = a+b;
       mrror = b-(mum-a);
       */
     }     

     __syncthreads();

     //d_j[(isp*gridDim.x+blockIdx.x)*(RHO_BIN_COUNT) + threadIdx.x] =dc_q[isp]*sum*ihx*idt;
     //d_j[(isp*gridDim.x+blockIdx.x)*(RHO_BIN_COUNT) + threadIdx.x] =dc_q[isp]*sum;
     d_j[(isp*gridDim.x+blockIdx.x)*(RHO_BIN_COUNT) + threadIdx.x] =sum;
     d_j[(isp*gridDim.x+blockIdx.x)*(RHO_BIN_COUNT)+MS + threadIdx.x] =mum;
     //if(threadIdx.x==0)
     //  printf("%e, %e\n",sum,dc_q[isp]*sum*ihx*idt);
     
   }  

 }
/* end  move_pcles_acc_cn() */


//particle move for a sub-time step (single-precision)
__device__ void acc_cn_incell_offset_sp(int isp, float x_n,float v_n,int i_n, float dt, float E0, float E1, float &x_np, float &v_np, int &i_np, float &x_nh, int &i_nh,float &dt_sub,int pi)
{
  dt_sub = dt;
  //Source terms from last time step

  float N = x_n;
  float D = hx_sp;
  float R = ihx_sp;
  float x = N*R;        
  float r = N - x*D;
  x += r*R;

  float a_n;
  a_n = E0 - E0*x + E1*x;
  a_n*= dc_q_m_sp[isp]; //!
  float a_nd = a_n*2.0f;
  float v_nd = v_n*2.0f;

  //estimate sub-time-step 
  N = (E1-E0);
  R = ihx_sp;
  D = hx_sp;
  float dadx = N*R;
  r = N - dadx*D;
  dadx += r*R;
  dadx *= dc_q_m_sp[isp]; //!

  float  aa = fabsf(a_n)+fabsf(dadx*v_n);
  float  bb = etolr_sp*(fabsf(v_n) + fabsf(a_n));
  float  cc = etola_sp;

      
  aa = rsqrtf(aa);
  dt_sub = cc*(aa);
  bb = bb*aa*aa*2.0f;  
  float dt_sube = fmaxf(dt_sub, bb);

  /*
  bb = bb/aa;
  //bb = __fdividef(bb,aa);
  float dt_sube = fmaxf(cc, bb*2.0f);
  */
  dt_sub = fminf(dt_sube,dt);
  
  //single step cn move                                                     
  float dth_sub = 0.5f*dt_sub;            //1
  float dth_sub2 = dth_sub*dth_sub;
  D = 1.0f-dadx*dth_sub2;
  N = dadx*dth_sub*v_n + a_n;

  float dvbydt = __fdiv_rn(N,D);
  v_np = v_n + dvbydt*dt_sub;

  
  /*  
  asm("rcp.approx.f32 %0, %1;\n\t"   //a = rcp(c) approximate       
  : "=f"(R) :  "f"(D));
  //R= _rcp(D);

  v_np = N*R;
  r = N - v_np*D; //Newton-Raphson     
  v_np = v_np + r*R;      //necessary
  */
  /*
  float dx = N*R;
  float r = N - dx*D; //Newton-Raphson
  dx = dx + r*R;
  */



  x_np = x_n + dth_sub*v_n + dth_sub*v_np; 
  float dx = x_np - x_n;
  i_np = i_n;

  float amd = a_n + E0*dc_q_m_sp[isp] + dadx*x_np; //!

  float vmd,sign_dx;

  bool cross = (x_np<0) || (x_np>hx_sp)  ;
  dt_sube = dt_sub;
  if( cross){

     if(x_np<=0){
        dx = -x_n;
        x_np = hx_sp;
	i_np--;
        sign_dx = -1.0f;
	amd = a_n+E0*dc_q_m_sp[isp]; //!
     }else{
        dx = hx_sp-x_n;
        x_np = 0.0f;
	i_np++;
        sign_dx = 1.0f;
	amd = a_n+E1*dc_q_m_sp[isp]; //!
      }
     i_np = ((i_np-1) & (Nx-1)) +1;           //,4
     
     float vn2 = v_n*v_n;
     float vnp2 = vn2 + amd*dx;
     float vnp2_2 = vn2 - v_n*v_n;
     vnp2 -= vnp2_2;

     if(vnp2<=0.0f){ //careful
       v_np = 0.0f; 
       i_np = i_n;
       if(x_np==0.0f){
	 x_np = hx_sp;
       }else{
	 x_np = 0.0f;
       }
     }else{
       v_np = sign_dx*__fsqrt_rn(vnp2);
       vmd = (v_n+v_np);
       if(vmd==0.0f){                                   //1
	 dt_sub = __fdiv_rn((v_np-v_n),amd)*2.0f;
       }else{
	 dt_sub = __fdiv_rn(dx,vmd)*2.0f;     
       } 

     /*
     float ivnp = _rsqrt(vnp2);
     dth_sub = __fdividef(-v_n*ivnp + sign(dx), amd*ivnp);
     v_np = _fma(amd,dth_sub,v_n);
     float F =  _fma(dth_sub,v_np, -dx);    //4(1)
     F = _fma(dth_sub,v_n,F);
     dth_sub -= __fdividef(F, v_np)*_fconst(0.5);
     v_np = _fma(amd,dth_sub,v_n);

     dt_sub = _fconst(2.0)*dth_sub;
     */


     }

   }else if( (x_np==0)&&(v_np<0) ){
     i_np--; 
     i_np = ((i_np-1) & (Nx-1)) +1; 
     x_np = hx_sp;
   }else if( (x_np==hx_sp)&&(v_np>0) ){
     i_np++;
     i_np = ((i_np-1) & (Nx-1)) +1;
     x_np = 0;
  }

  i_nh = i_n;
  x_nh = x_n + dx*0.5f;

  dt_sub = fmin(dt_sub, dt_sube);
  dth_sub = dt_sub*0.5f;
  v_np = v_n + amd*dth_sub;
  

  // if(v_np!=v_np)
  //   printf("%d %24.12e %24.12e %24.12e %24.12e %24.12e %d %d\n",pi,v_n,v_np,x_n, x_np, dt_sub, i_n,i_np);

}
/* end acc_cn_incell_offset_sp() */


/* Crank-Nicolson mover (adaptive-charge-conserving)*/
__global__ void move_pcles_acc_cn_sp(Particles pcles, int isp, float_p *d_E,float_p dt,float_p *d_j)
 {
   const float idt = __frcp_rn(__double2float_rn(dt));

   //initialize shared memory
   extern __shared__ float ss_j[];


   if(threadIdx.x<Nx*G_WARP){
     ss_j[threadIdx.x] = _fconst(0.0);
   }
   __syncthreads();
   float *w_j = ss_j + (threadIdx.x>>LOG_Nx)*Nx;

   for(int pi=dc_acNpcls[isp]+UMAD(blockIdx.x, blockDim.x, threadIdx.x); pi<dc_acNpcls[isp+1]; pi+=UMUL(blockDim.x,gridDim.x) ){
     float dtp = dt;
     float dt_sub;
     //save old positions and velocities
     float x_n  = __double2float_rn(pcles.d_x_n [pi]);
     float v_n  = __double2float_rn(pcles.d_v_n [pi]);
     int   i_n  = pcles.d_ci_n[pi]; 

     float x_np,v_np,x_nh;
     int i_np,i_nh;
     float vh;

     float j0=0.0f,j1=0.0f;

     float E0 =  __double2float_rn(d_E[i_n-1]);
     float E1 =  __double2float_rn(d_E[i_n  ]);

     int it = 0;
     //while(1){
     while(dtp>0){
       acc_cn_incell_offset_sp(isp,x_n,v_n,i_n,dtp,E0,E1,x_np,v_np,i_np,x_nh,i_nh,dt_sub,pi);
       vh = (v_n+v_np)*0.5f;
       
       float q = vh*dt_sub;                     //1                  
       //float_p x = x_nh*ihx_sp; 
       float N = x_nh;
       float R = ihx_sp;
       float D = hx_sp;
       float x = N*R; //x_nh*ihx;
       float r = N - x*D;
       x += r*R;

       j0 += (q-q*x);   
       j1 += (q*x);

       if(i_np!=i_n){      
	 atomicAdd(w_j+i_nh-1       , j0);   
	 atomicAdd(w_j+(i_nh&(Nx-1)), j1);  

	 j0 = 0.0f; 
	 j1 = 0.0f;

	 E0 =  __double2float_rn(d_E[i_np-1]);
	 E1 =  __double2float_rn(d_E[i_np  ]);

       }

       //it++;
       //if(it>400) break;


       x_n = x_np;
       v_n = v_np;
       i_n = i_np;
       dtp -= dt_sub;


     }
  
    //write back to global memory  
     pcles.d_x_np[pi] = double(x_np);
     pcles.d_v_np[pi] = double(v_np);
     pcles.d_ci_np[pi]= i_np;

     atomicAdd(w_j+i_nh-1       , j0);
     atomicAdd(w_j+(i_nh&(Nx-1)), j1);


   }
   __syncthreads();
   //reduce
   if(threadIdx.x < RHO_BIN_COUNT){ //careful

     float sum = 0.0f;
     for(int i=0; i<G_WARP; i++){
       sum +=  ss_j[threadIdx.x+i*Nx];
     }     

     __syncthreads();

     //d_j[(isp*gridDim.x+blockIdx.x)*(RHO_BIN_COUNT) + threadIdx.x] =dc_q[isp]*sum*ihx*idt;
     //d_j[(isp*gridDim.x+blockIdx.x)*(RHO_BIN_COUNT) + threadIdx.x] =dc_q[isp]*sum;
     d_j[(isp*gridDim.x+blockIdx.x)*(RHO_BIN_COUNT) + threadIdx.x] =double(sum);
     d_j[(isp*gridDim.x+blockIdx.x)*(RHO_BIN_COUNT)+MS + threadIdx.x] =0.0;
     //if(threadIdx.x==0)
     //  printf("%e, %e\n",sum,dc_q[isp]*sum*ihx*idt);
     
   }  

 }
/* end  move_pcles_acc_cn_sp() */


/* Crank-Nicolson mover */
 __global__ void move_pcles_cn(Particles pcles, int isp, float_p *d_E,float_p dt,float_p *d_j)
 {
   const int threadPos = threadIdx.x;
   __shared__ float_p s_j[SMEMSIZE];
   //initialize shared memory
   // for(int xi=0; xi<RHO_BIN_COUNT; xi++ ){ //careful
   //   s_j[threadIdx.x + xi*BLOCK_SIZE ] = 0.;
   // }
   for(int xi=threadIdx.x; xi<SMEMSIZE; xi++){
     s_j[xi] = 0.;     
   }
   __syncthreads();

   for(int pi=dc_acNpcls[isp]+UMAD(blockIdx.x, blockDim.x, threadIdx.x); pi<dc_acNpcls[isp+1]; pi+=UMUL(blockDim.x,gridDim.x) ){
     //save old positions and velocities
     float_p x_n  = pcles.d_x_n [pi];
     float_p v_n  = pcles.d_v_n [pi];
     int   i_n  = pcles.d_ci_n[pi]; 

     float_p a_n  = dc_q_m[isp]*fld_on_pcleh(d_E,i_n,(x_n*ihx-i_n+1.5));

     float_p x_np,v_np,a_np;
     int i_np;

     // Initialize nonlinear iteration
     float_p dx_pn = fabs(dt*a_n);
     float_p dx_pk = dx_pn;
     ushort   it = 0;
     float_p x_k, dx;
     v_np = v_n;
     //     while(dx_pk>eps*dx_pn && it<200){
     while(it<5){
       x_k  = v_np;
       dx   = (v_np+v_n)*dt*0.5;
       x_np = x_n + dx;
       //B.C.
       x_np = par_bc_pcles(x_np);
       i_np = locate(x_np);
       a_np = dc_q_m[isp]*fld_on_pcleh(d_E,i_np,(x_np*ihx-i_np+1.5));
       v_np = v_n + (a_n+a_np)*0.5*dt;

       //Check convergence
       dx_pk = fabs(v_np - x_k);
       it++;
     }
  
    //write back to global memory  
     pcles.d_x_np[pi] = x_np;
     pcles.d_v_np[pi] = v_np;
     pcles.d_ci_np[pi]= i_np;
     // if(i_np!=i_n){

     // }

     float_p vh  = (v_n+v_np)*0.5;

     s_gatherj(s_j+threadPos,vh,i_n ,(x_n *ihx-i_n +1.5));  
     s_gatherj(s_j+threadPos,vh,i_np,(x_np*ihx-i_np+1.5));      
     //cuPrintf("[%d %d] %f %f %d, %f %f %d (%d)\n",isp,pi,x_n,v_n,i_n,x_np,v_np,pcles.d_ci_np[pi],it);
     //cuPrintf("%d %d %d %f %f %d %f %f\n",isp,pi,i_n,*(s_j+threadPos+(i_n-1)*BLOCK_SIZE),*(s_j+threadPos+(i_n&(Nx-1))*BLOCK_SIZE),i_np,*(s_j+threadPos+(i_np-1)*BLOCK_SIZE),*(s_j+threadPos+(i_np&(Nx-1))*BLOCK_SIZE) );
   }
   __syncthreads(); //careful
   
   //reduce
   if(threadIdx.x < RHO_BIN_COUNT){ //careful
     float_p *s_rhoBase = s_j + UMUL(threadIdx.x, BLOCK_SIZE);
     float_p sum = 0;
     int pos;
     for(int j=0; j<G_SHARED_BANKS; j++){
       pos = threadIdx.x & (SHARED_MEMORY_BANKS-1);//to avoid bank conflict
       for(int i=0; i<SHARED_MEMORY_BANKS; i++){
 	sum += s_rhoBase[j*SHARED_MEMORY_BANKS+pos];
 	pos = (pos + 1) & (SHARED_MEMORY_BANKS-1);
	// if(threadIdx.x==31)
	//    cuPrintf("%d %d %d %d %f\n",threadIdx.x,j,i,pos,sum);
       }
     }
  
     d_j[(isp*gridDim.x+blockIdx.x)*(RHO_BIN_COUNT) + threadIdx.x] =dc_q[isp]*sum*0.5*ihx;
   }  

 }
/* end  move_pcles_cn() */

//update particle quantities for the next time step
__global__ void update_tnp_pcles(Particles pcles)
{
  //process all particles
  for(int pi=UMAD(blockIdx.x, blockDim.x, threadIdx.x); pi<totNp; pi+=UMUL(blockDim.x,gridDim.x) ){
    pcles.d_x_n [pi] = pcles.d_x_np [pi];
    pcles.d_v_n [pi] = pcles.d_v_np [pi];
    pcles.d_ci_n[pi] = pcles.d_ci_np[pi];
  }
}
/* end update_tnp_pcles() */


//sort v only
 __global__ void sortv_particles(Particles pcles,int isp) 
 {//after tmp/pic_kernel.cuh-tmp15

   __shared__ int part_cell_index[binSize_v*(maxThreads+1)];
   float_p idv = _fconst(4.0); //(1.0); //
   /*
  0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 
----:---:---:---:---:---:---:-----
   -3  -2  -1   0   1   2   3
    */
   int tot_binSize_v = binSize_v*(maxThreads+1);
   for(int bi=threadIdx.x; bi< tot_binSize_v; bi+=blockDim.x){
     part_cell_index[bi] = 0;
   }
   __syncthreads();

   //scan first
   float_p v_np;
   int ci;
   int *sbin = part_cell_index + threadIdx.x*binSize_v; //each thread has its own domain (row oriented)
   
  for(int pi=dc_acNpcls[isp]+UMAD(blockIdx.x, blockDim.x, threadIdx.x); pi<dc_acNpcls[isp+1]; pi+=UMUL(blockDim.x,gridDim.x) ){

    v_np = _abs(pcles.d_v_np[pi])/dc_vth[isp];

    if(v_np>_fconst(3.0)){
      ci = binSize_v-1; 
    }else{
      ci = _floor(v_np*idv);
    }
    sbin[ci]++;
    //printf("%e\t %d %d %d\n",v_np,pi,ci,sbin[ci]);
  }
  __syncthreads();

  

  //prefix sum
  int *accum_bin;
  if(threadIdx.x<binSize_v){
    accum_bin = part_cell_index + threadIdx.x; 
    int ti_0 = accum_bin[0];
    int ti_1;
    accum_bin[0] = 0;
    //#pragma unroll 
    for(int i=1; i<maxThreads+1; i++){ //column oriented
      ci = i*binSize_v;
      ti_1 = accum_bin[ci];
      accum_bin[ci] = ti_0 + accum_bin[ci-binSize_v];
      ti_0 = ti_1;
      //if(threadIdx.x==3)printf("%d\t %d %d\n",ci, ti_1,  accum_bin[ci]);
    }
  }
  __syncthreads();

  
  if(threadIdx.x==0){
    accum_bin = part_cell_index + (maxThreads)*binSize_v;
    int ti_0 = accum_bin[0];
    int ti_1;
    ci = 0;
    accum_bin[0] = 0;

    //#pragma unroll
    for(int i=1; i<binSize_v; i++){ //row oriented
      ci = i;
      ti_1 = accum_bin[ci];
      accum_bin[ci] = ti_0 + accum_bin[ci-1];
      ti_0 = ti_1;

    }    

  }
  __syncthreads();


  //copy next
  int pi_old;
  for(int pi=dc_acNpcls[isp]+UMAD(blockIdx.x, blockDim.x, threadIdx.x); pi<dc_acNpcls[isp+1]; pi+=UMUL(blockDim.x,gridDim.x) ){
    v_np = _abs(pcles.d_v_np[pi])/dc_vth[isp];
    if(v_np>_fconst(3.0)){
      ci = binSize_v-1; 
    }else{
      ci = _floor(v_np*idv);
    }

    pi_old = sbin[ci] + part_cell_index[(maxThreads)*(binSize_v)+ci] + dc_acNpcls[isp];


    pcles.d_x_n [pi_old] = pcles.d_x_np [pi];
    pcles.d_v_n [pi_old] = pcles.d_v_np [pi];
    pcles.d_ci_n[pi_old] = pcles.d_ci_np[pi];

    
    sbin[ci]++;

  }  


 }
/* end sortv_particles() */
