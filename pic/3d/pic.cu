#include <unistd.h>
#include <assert.h>
#include <stdio.h>
//#include "cudpp/cudpp.h"
#include "pic_kernel.cuh"
#include "pic.h"
#include "rand.h"

float_p *h_charge;
float_p *h_mass;
int   h_Nsp;
int tot_Np;
float_p h_LL ;
float_p h_hx;
float_p h_ihx;
float_p h_xx0;
float_p h_xf0;
const float_p zero = 0.;


//const float_p h_etola = _fconst(0.015625); //22
//const float_p h_etola = _fconst(0.03125); //22
//const float_p h_etola = _fconst(0.0625); 
//const float_p h_etola = _fconst(0.0078125); 
const float_p h_etola = _fconst(0.0001); 
//const float_p h_etola = _fconst(0.0009765625); 
//const float_p h_etolr = _fconst(0.03); 
//const float_p h_etolr = _fconst(0.03125); 
const float_p h_etolr = _fconst(0.0625); 
//const float_p h_etolr = _fconst(0.015625); 
//const float_p h_etolr = _fconst(0.0078125); 

//const float_p h_etolr = _fconst(0.01); 

const float_p h_eps   =  1e-300; //1e-20; //1e-8; //_fconst(0.0000000001);
const float_p h_PI = _fconst(3.14159265358979323846);

//pointers in gpu global memory
float_p *d_mv,*d_v2; 
float_p *d_E;     //electric field
float_p *d_rho;   //charge density
float_p *d_j;    //current density
float_p *d_xx;    //spatial grid 

int   *d_acNpcls; //number of particles of species (accumulated)
float_p *d_q_m;   //q/m of species
float_p *d_q;     //q of species
float_p *d_m;     //m of species
cudaEvent_t     start, stop;

extern "C" {

// void h_move_pcles_acc_cn(Particles pcles, int isp, float_p *E,float_p dt,float_p *j)
// {

//}


  void checkCUDAError(const char *msg)
  {
    cudaError_t err = cudaGetLastError();
    if( cudaSuccess != err) 
      {
        fprintf(stderr, "Cuda error: %s: %s.\n", msg, 
		cudaGetErrorString( err) );
        exit(EXIT_FAILURE);
      }                         
  }
  /* end checkCUDAError() */

  void device_start_timing()
  {
    cudaEventCreate( &start );
    cudaEventCreate( &stop  );
    cudaEventRecord( start, 0 );    
  }
  /* end device_start_timing() */

  void device_stop_timing()
  {
    cudaEventRecord( stop, 0 );
    cudaEventSynchronize( stop );
    float   elapsedTime;
    cudaEventElapsedTime( &elapsedTime, start, stop );
    printf( "PIC GPU Time:  %3.1f ms\n", elapsedTime );

    cudaEventDestroy( start ) ;
    cudaEventDestroy( stop  ) ;
  }
  /* end device_stop_timing() */

void getNBT_re(int n, int &blocks, int &threads)
{        
  if(n==2){
    blocks  = 1;
    threads = 2;
  }else{
    threads = (n < maxThreads*2) ? nextPow2((n + 1)/ 2) : maxThreads;
    blocks  = (n + (threads * 2 - 1)) / (threads * 2);
    blocks  = MIN(maxBlocks, blocks);
  }
}
/* end getNBT_re() */

  int h_locate(float_p xp)
  {
#if   _SPORDER == 1
    return      floor((xp-h_xx0)*h_ihx);
#elif _SPORDER == 2
    return       ceil((xp-h_xf0)*h_ihx);
#else 
    printf("Spline order (%d) unavailable in h_locate(),\nEXIT.\n",_SPORDER);
    exit(1);
#endif
  }

//particle B.C.
float_p h_par_bc_pcles(float_p xp)
{
#if   _SPORDER == 1  
  return fmod((xp-h_xx0+h_LL),h_LL)+h_xx0;  
#elif _SPORDER == 2
  xp = fmod((xp-h_xf0+h_LL),h_LL)+h_xf0;  //[xf0,xf0+LL)
  int is_xf0 = (xp==h_xf0);
  return xp + (float_p)(is_xf0)*h_LL; //(xf0,xf0+LL] 
#endif

}
/* end h_par_bc_pcles() */

  void alloc_particles(Particles &pcles,int memSize)
  {
    cudaMalloc((void **)&pcles.d_x_n  , memSize);
    cudaMalloc((void **)&pcles.d_x_np , memSize);
    cudaMalloc((void **)&pcles.d_v_n  , memSize);
    cudaMalloc((void **)&pcles.d_v_np , memSize);
    cudaMalloc((void **)&pcles.d_ci_n , memSize);
    cudaMalloc((void **)&pcles.d_ci_np, memSize);    
  }

  void free_particles(Particles &pcles)
  {
    cudaFree(pcles.d_x_n  );
    cudaFree(pcles.d_x_np );
    cudaFree(pcles.d_v_n  );
    cudaFree(pcles.d_v_np );
    cudaFree(pcles.d_ci_n );
    cudaFree(pcles.d_ci_np);
  }

  void device_finalize()
  {    
    free(h_charge);
    free(h_mass);
#ifdef _USE_TEX
    cudaUnbindTexture(d_xxTex);
    cudaUnbindTexture(d_ETex);
#endif
    cudaFree(d_j    );
    cudaFree(d_E    );
    cudaFree(d_rho  );
    cudaFree(d_xx   );

    //new 
    free_particles(pcles);

    cudaFree(d_acNpcls);
    cudaFree(d_q_m  );
    cudaFree(d_q    );
    cudaFree(d_m    );

    cudaFree(d_mv);
    checkCUDAError("Free device pointers");
    
  }
  /* end device_finalize() */

  void  device_print_msg()
  {
    cudaPrintfInit();
    int tBlockSize,tGridSize;
    if(tot_Np==2){
      tGridSize  = 1;
      tBlockSize = 2;
    }else{
      tGridSize  = GRID_SIZE;
      tBlockSize = BLOCK_SIZE;
    }

    device_greetings<<<tGridSize,tBlockSize>>>(d_xx,pcles,0.1);
    checkCUDAError("kernel invocation");
    cudaPrintfDisplay();
    cudaPrintfEnd();   
  }
  /*end device_print_msg() */



  //setup grid and some parameters
  void device_setup_pic(int nx,float_p xx[],bool &explicit_pic,bool &quiet_start, float_p q_m[],float_p v0[],float_p vth[], int npcles[],int &n_sp)
  //  void device_setup_pic(int nxin,float_p xin[],bool explicit,bool quiet_start, float_p q_m[], int npcls[],int nsp)
  {
    h_Nsp = n_sp;

    //printf("#explicit=%s,quiet_start=%s\n",explicit_pic?"t":"f",quiet_start?"t":"f");
    // for(int i=0;i<n_sp;i++){
    //   printf("#spes:%d, q/m=%f, npcles=%d, v0=%f, vth=%f\n\n\n",i,q_m[i],npcles[i],v0[i],vth[i]);
    // }

    
    //copy xx to device
    assert(nx%2==0); //ensure power of 2 (periodic B.C.)
    int memSize = sizeof(float_p)*(nx+2);
    cudaMalloc((void **)&d_xx, memSize);        
    cudaMemcpy(d_xx,xx,memSize,cudaMemcpyHostToDevice);      
    checkCUDAError("copy grid (xx) memory");
#ifdef _USE_TEX
    //texture for grid 
    cudaBindTexture( NULL, d_xxTex, d_xx, memSize );
    checkCUDAError("bind grid (xx) texture memory");
#endif

    //set grid parameters
    h_LL = xx[nx] - xx[0];
    h_hx  = xx[1]-xx[0];
    h_ihx = _fconst(1.0)/h_hx;
    h_xf0= _fconst(0.5)*(xx[1]+xx[0]);
    h_xx0= xx[0];
    //printf("#h=%26.16e, ih=%26.16e, xf0=%26.16e, LL=%26.16e\n",h_hx,h_ihx,h_xf0,h_LL);
    //exit(1);
    //single-precision const.
    float h_etolr_sp = float (h_etolr);
    float h_etola_sp = float (h_etola);
    float h_hx_sp = float (h_hx);
    float h_ihx_sp = float (h_ihx);
    cudaMemcpyToSymbol("etolr_sp",&h_etolr_sp, sizeof(float),
		       0,cudaMemcpyHostToDevice);
    cudaMemcpyToSymbol("etola_sp",&h_etola_sp, sizeof(float),
		       0,cudaMemcpyHostToDevice);
    cudaMemcpyToSymbol("hx_sp",&h_hx_sp, sizeof(float),
		       0,cudaMemcpyHostToDevice);
    cudaMemcpyToSymbol("ihx_sp",&h_ihx_sp, sizeof(float),
		       0,cudaMemcpyHostToDevice);

    cudaMemcpyToSymbol("etolr",&h_etolr, sizeof(float_p),
		       0,cudaMemcpyHostToDevice);
    cudaMemcpyToSymbol("etola",&h_etola, sizeof(float_p),
		       0,cudaMemcpyHostToDevice);
    cudaMemcpyToSymbol("eps",&h_eps, sizeof(float_p),
		       0,cudaMemcpyHostToDevice);
    cudaMemcpyToSymbol("PI",&h_PI, sizeof(float_p),
		       0,cudaMemcpyHostToDevice);
    cudaMemcpyToSymbol("LL",&h_LL, sizeof(float_p),
		       0,cudaMemcpyHostToDevice);
    cudaMemcpyToSymbol("hx",&h_hx, sizeof(float_p),
		       0,cudaMemcpyHostToDevice);
    cudaMemcpyToSymbol("ihx",&h_ihx, sizeof(float_p),
		       0,cudaMemcpyHostToDevice);
    cudaMemcpyToSymbol("xf0",&h_xf0, sizeof(float_p),
		       0,cudaMemcpyHostToDevice);
    cudaMemcpyToSymbol("xx0",&h_xx0, sizeof(float_p),
		       0,cudaMemcpyHostToDevice);
    cudaMemcpyToSymbol("Nx",&nx, sizeof(int),
		       0,cudaMemcpyHostToDevice);
    checkCUDAError("set constant memory");

    //Find PIC physical quantities
    float_p omega, n0;
    h_charge=(float_p *)malloc(n_sp*sizeof(float_p));
    h_mass  =(float_p *)malloc(n_sp*sizeof(float_p));
    
    for(int i=0; i<n_sp; i++){
      //omega = sqrt(abs(q_m[i]));
      //n0        = float_p(npcles[i])/h_LL;
      h_charge[i] = (q_m[i]>0?_fconst(1.0):-_fconst(1.0))*h_LL/npcles[i];
      h_mass[i]   = h_charge[i]/q_m[i];
      printf("sp%d: charge=%26.16e, mass=%26.16e \n",i, h_charge[i],h_mass[i]);
    }
    //exit(1);
    //single-precision constants
    float *q_m_sp = (float *) malloc(n_sp*sizeof(float));
    for(int i=0; i<n_sp; i++){
      q_m_sp[i] =  float(q_m[i]);
    }
    printf("%f %f\n",q_m_sp[0],q_m_sp[1]);
    
    cudaMemcpyToSymbol("dc_q_m_sp",q_m_sp, sizeof(float)*n_sp,
		       0,cudaMemcpyHostToDevice);
    free(q_m_sp);

    //q/m of particle species
    cudaMemcpyToSymbol("dc_q_m",q_m, sizeof(float_p)*n_sp,
		       0,cudaMemcpyHostToDevice);
    cudaMemcpyToSymbol("dc_q",h_charge, sizeof(float_p)*n_sp,
		       0,cudaMemcpyHostToDevice);
    cudaMemcpyToSymbol("dc_m",h_mass, sizeof(float_p)*n_sp,
		       0,cudaMemcpyHostToDevice);
    cudaMemcpyToSymbol("dc_vth",vth, sizeof(float_p)*n_sp,
    		       0,cudaMemcpyHostToDevice);


     memSize = sizeof(float_p)*n_sp;
     cudaMalloc((void**)&d_q_m, memSize);
     cudaMemcpy(d_q_m,q_m   ,memSize,cudaMemcpyHostToDevice);      
     //q
     cudaMalloc((void**)&d_q, memSize);
     cudaMemcpy(d_q,h_charge,memSize,cudaMemcpyHostToDevice);      
     //m
     cudaMalloc((void**)&d_m, memSize);
     cudaMemcpy(d_m,h_mass  ,memSize,cudaMemcpyHostToDevice);      


    //load particles
     device_load_pcles(v0,vth,npcles,n_sp,xx,nx);


    
  }
  /* end device_setup_pic() */


  void device_load_pcles(float_p v0[],float_p vth[],int numParticles[], int n_sp, float_p xx[], int nx) 
   {
     //set particle numbers
     tot_Np = 0;
     int *acNpcles=(int *)malloc((n_sp+1)*sizeof(int));
     acNpcles[0] = 0;
     for(int i=0; i<n_sp; i++){
       tot_Np += numParticles[i];
       acNpcles[i+1] = tot_Np; //accumulated number
     }

     //number of particles of species
     int memSize = sizeof(int)*(n_sp+1);
     cudaMalloc((void**)&d_acNpcls, memSize);
     cudaMemcpy(d_acNpcls,acNpcles,memSize,cudaMemcpyHostToDevice); 
     cudaMemcpyToSymbol("dc_acNpcls",acNpcles, memSize,
			0,cudaMemcpyHostToDevice);

     
     cudaMemcpyToSymbol("dc_Npcls",numParticles, sizeof(int)*n_sp,
			0,cudaMemcpyHostToDevice);
     
     
     cudaMemcpyToSymbol("totNp",&tot_Np, sizeof(int),
			0,cudaMemcpyHostToDevice);
     cudaMemcpyToSymbol("Nsp",&n_sp, sizeof(int),
			0,cudaMemcpyHostToDevice);

     checkCUDAError("set constant memory");
     
     //allocate GPU memory             
     //new 
     memSize = sizeof(float_p)*tot_Np*DIM_PIC; 
     int memSize_i = sizeof(int)*tot_Np*DIM_PIC; //cin, cinp
     alloc_particles(pcles,memSize);
     checkCUDAError("Allocate Particles on device");

     //allocate host memory
     float_p *h_x=(float_p *)malloc(memSize);
     float_p *h_v=(float_p *)malloc(memSize);
     int *h_i=(int *)malloc(memSize_i);


     //FILE *fp = fopen("xv.dat","w");
     int ip ,ip_x,ip_v, ci;     
     float_p v_avg = 0;
     float_p dx,xp_off,xp_off0;
     int seed, seedn;

     if(tot_Np==2){
       //2 particles
       h_x[0] = 0.6; //x
       h_x[0] = 0.0; //v
       h_v[1] = 0.8; //x
       h_v[1] = 0.0; //v       
       h_i[0] = h_locate(h_x[0]); //ci_n, ci_np
       h_i[1] = h_locate(h_x[1]);
     }else{
       double rx;
       double rx0=0.0; //1.0;
       ip = 0;
       for(int is=0; is<n_sp; is++){
       seed  = 1;
       seedn = 1;

#if _SPORDER%2==0
	  ci = 1;
#else
	  ci = 0;
#endif
	  int Ni=numParticles[is]/nx;
	  dx = h_LL/numParticles[is];
	  xp_off0 = xp_off = 0.;

	  for(int i=0; i<numParticles[is]; i++){

#ifndef _OFFSET
	   printf("use _OFFSET. exit.\n");
	   exit(1);
#else
	   	   
	   //if(((i+1) & 1) == 0){
	     xp_off+= dx;
	     //}

	   if(xp_off>h_hx){ //careful
	     xp_off -= h_hx;
	     ci++;
	   //B.C. (peroidic, assuming nxx power of 2)
# if _SPORDER%2==0
	     ci = ((ci-1)&(nx-1)) + 1; //from 1 to nx
# else
	     ci = (ci)&(nx-1); //from 0 to nx-1
# endif
	   }
	   /*
# if _SPORDER%2==0
 	     ci = ((ci-1)&(nx-1)) + 1; //from 1 to nx
# else
  	     ci = (ci)&(nx-1); //from 0 to nx-1
# endif
	     //LCG_rand(rx,seed);
	     //Randoms(rx,seed);	     
	     //lcg_rand(rx0,seed);
	     rx0 = float_p(rand())/float_p(RAND_MAX+1);
	     xp_off = h_hx*float_p(rx0);
	   */

 	   h_i[ip] = ci; 
 	   h_x[ip] = xp_off;	   	   
#endif
	   //set velocity
	   //if((i & 1) == 0){
	     randn(rx,seedn);
	     h_v[ip] = v0[is] + vth[is]*float_p(rx); 
	   // }else{ //quiet
	   //   h_v[ip] = h_v[ip-1];
	   // }

	    v_avg += h_v[ip]*h_v[ip];

	    //ci = rand()%nx + 1; //1 to nx
	    //if((i & 1) == 1){
	    // if(((ip+1)%Ni) == 0){
	    //    ci++;
	    // }
	   ip++;
	  }
	  //printf("sum_cpu = %e\n",v_avg*0.5*h_mass[is]);   
	  //v_avg = 0;
       } 



       /*
       for(int is=0; is<n_sp; is++){
	 for(int i=0; i<numParticles[is]; i++){
	   h_x[ip] = h_xx0 + 1e-3 + i*h_LL/numParticles[is]; //x (careful)
	   randn(rx);
	   h_v[ip] = v0[is] + vth[is]*float_p(rx); //v
	   //h_v[ip] = float_p(rx); //v
	   //B.C.
	   h_x[ip] = h_par_bc_pcles(h_x[ip]);
	   h_i[ip] = h_locate(h_x[ip]);//ci_n&np
#ifdef _OFFSET
# if _SPORDER%2==0
	   h_i[ip] = ((h_i[ip]-1) & (nx-1)) + 1; //from 1 to nx
	   h_x[ip] = h_x[ip] - xx[h_i[ip]] + 0.5*h_hx;
# else
	   //peroidic, assuming nxx power of 2           
	   h_i[ip] = (h_i[ip] & (nx-1)); //from 0 to nx-1
	   //offset to the nearest grid point
	   h_x[ip]-= xx[h_i[ip]];
# endif
#endif
	   //fprintf(fp,"%d %d %f %f %d\n",i, ip,h_x[ip],h_v[ip],h_i[ip]);
	   // v_avg += h_v[ip]*h_v[ip];
	   ip++;	   
	 }
	 //fprintf(fp,"\n\n");
	 // printf("sum_cpu = %e\n",v_avg*0.5*h_mass[is]);
	 // v_avg = 0.;
       }
       */
       //printf("#Initial x&v(tot_Np>2) completed.\n");
       
     }

     //copy the particle x&v to device (only this once)
     //new
     cudaMemcpy(pcles.d_x_n  ,h_x,memSize  ,cudaMemcpyHostToDevice);
     cudaMemcpy(pcles.d_x_np ,h_x,memSize  ,cudaMemcpyHostToDevice);
     cudaMemcpy(pcles.d_v_n  ,h_v,memSize  ,cudaMemcpyHostToDevice);
     cudaMemcpy(pcles.d_v_np ,h_v,memSize  ,cudaMemcpyHostToDevice);
     cudaMemcpy(pcles.d_ci_n ,h_i,memSize_i,cudaMemcpyHostToDevice);
     cudaMemcpy(pcles.d_ci_np,h_i,memSize_i,cudaMemcpyHostToDevice);

     checkCUDAError("copy particle (ci) memory");
     
     //new 
     free(h_x );
     free(h_v );
     free(h_i);

     free(acNpcles);
     //fclose(fp);

     //sleep(3);

   }
   /* end device_load_pcles() */

  void device_allocate_field_memory(int nx)
  {

    size_t memSize = sizeof(float_p)*(nx+2);
    cudaMalloc((void **)&d_E  , memSize);        
#ifdef _USE_TEX
    //texture for field
    cudaBindTexture( NULL, d_ETex, d_E, memSize );
    checkCUDAError("bind field (E) texture memory");
#endif
    checkCUDAError("allocate field (E) device memory");

    //each species has a set of rho for all blocks (careful)
    //    memSize = sizeof(float_p)*(nx+1)*h_Nsp*GRID_SIZE;
    memSize = sizeof(float_p)*(nx)*h_Nsp*GRID_SIZE;

    cudaMalloc((void **)&d_rho, memSize);            
    cudaMemset( d_rho, zero, memSize );
    checkCUDAError("allocate charge density (rho) device memory");

    memSize *= 2; //p-m
    cudaMalloc((void **)&d_j, memSize);        
    cudaMemset( d_j, zero, memSize );
    checkCUDAError("allocate current density (j) device memory");

  }
  /* end device_allocate_field_memory() */

  void device_find_current(float_p j[], int nx)
  {
    int tBlockSize,tGridSize;
    if(tot_Np==2){
      tGridSize  = 1;
      tBlockSize = 2;
    }else{
      tGridSize  = GRID_SIZE;
      tBlockSize = BLOCK_SIZE;
    }

    for(int isp=0; isp<h_Nsp; isp++){    
      find_current2<<<tGridSize,tBlockSize>>>(pcles, isp, d_j);
    }

    //merge d_j
    merge_d_rho<<<1,nx>>>(d_j);
    cudaMemcpy(j, d_j, sizeof(float_p)*(nx), cudaMemcpyDeviceToHost);
    checkCUDAError("find_current");

    //B.C.
    j[nx]   = j[0];
    j[nx+1] = j[1];


    // for(int i=0; i<nx+2; i++){
    //   printf("%d %f\n",i,j[i]);
    // }
    // exit(0);

  }
  /* end device_find_current() */

  void device_find_rho(float_p rho[], int nx)
  {
    //    printf("#in device_find_rho\n");
    cudaPrintfInit();
    int tBlockSize,tGridSize;
    if(tot_Np==2){
      tGridSize  = 1;
      tBlockSize = 2;
    }else{
      tGridSize  = GRID_SIZE;
      tBlockSize = BLOCK_SIZE;
    }

    for(int isp=0; isp<h_Nsp; isp++){    
      find_rho<<<tGridSize,tBlockSize>>>(pcles, isp, d_rho);
    }

    //merge d_rho
    merge_d_rho<<<1,nx>>>(d_rho);
    
    cudaMemcpy(rho, d_rho, sizeof(float_p)*(nx), cudaMemcpyDeviceToHost);
    
    //B.C.
    rho[nx]  = rho[0];
    rho[nx+1]= rho[1];
    

    cudaPrintfDisplay();
    cudaPrintfEnd();  

     // for(int i=0;i<nx+2;i++){
     //   printf("%d %f\n",i,rho[i]);
     // }
     // printf("\n\n");
     // exit(1);
  }
  /* end device_find_rho() */

  //perturb the particle positions
  void device_perturb_pcles(short nh1,float_p eps,float_p dt,float_p rho[],int nx)
  {
    //cudaPrintfInit();
    int tBlockSize,tGridSize;
    if(tot_Np==2){
      tGridSize  = 1;
      tBlockSize = 2;
    }else{
      tGridSize  = 1; //GRID_SIZE;
      tBlockSize = 32; //BLOCK_SIZE;
    }
    //eps = 0.;
    //perturb_pcles<<<tGridSize,tBlockSize>>>( nh1, eps, pcles, d_acNpcls, d_xx, dt, d_q, d_rho);
    size_t memSize = sizeof(float_p)*(nx)*h_Nsp*GRID_SIZE;

    for(int isp=0; isp<h_Nsp; isp++){
      perturb_pcles2<<<tGridSize,tBlockSize>>>( nh1, eps, pcles, isp, d_xx, dt, d_rho);
    }

    //merge d_rho
    merge_d_rho<<<1,nx>>>(d_rho);

    // cudaPrintfDisplay();
    // cudaPrintfEnd();  

    cudaMemcpy(rho, d_rho, sizeof(float_p)*(nx), cudaMemcpyDeviceToHost);
    cudaDeviceSynchronize();
    checkCUDAError("perturb particles");
    //B.C.
    rho[nx  ] = rho[0];
    rho[nx+1] = rho[1];



    //       for(int i=0; i<nx+2; i++){
    //   printf("%d %e\n",i,rho[i]);
    // }
    //exit(0);

  }
   /* end device__perturb_pcles() */

void device_set_field(float_p E[], int nx, int itime) 
  {
    //printf("#set device memory for E ...\n");
    int memSize = sizeof(float_p)*(nx+2);
    //copy E to device
    cudaMemcpy(d_E  ,E  ,memSize,cudaMemcpyHostToDevice);     
    checkCUDAError("copy E to d_E");
    
    //cudaMemcpy(d_rho,rho,memSize,cudaMemcpyHostToDevice);      
    //cudaMemcpy(rho, d_rho, memSize, cudaMemcpyDeviceToHost);
    //reset_mem<<<1,nx+2>>>( d_rho, nxp2, 0. );
				                         
    //    if(itime == 713){
    //for(int i=0;i<nx+2;i++){
    // printf("%d %24.16e\n",i,E[i]);
    //}
    // printf("#setupNonlinearFunction after mover.\n");     
    //printf("\n\n");
    //}
  }
  /* end device_set_field() */

  void device_move_pcles(float_p dt, float_p rho[], int nx)
  {
    //cudaPrintfInit();
    int tBlockSize,tGridSize;
    if(tot_Np==2){
      tGridSize  = 1;
      tBlockSize = 2;
    }else{
      tGridSize  = GRID_SIZE;
      tBlockSize = BLOCK_SIZE;
    }


    size_t memSize = sizeof(float_p)*(nx)*h_Nsp*GRID_SIZE;

    for(int isp=0; isp<h_Nsp; isp++){
#ifdef _OFFSET
      move_pcles_offset<<<tGridSize,tBlockSize>>>(pcles, isp, d_xx, d_E,dt, d_rho);
#else
      move_pcles2<<<tGridSize,tBlockSize>>>(pcles, isp, d_xx,d_q_m, d_E,dt, d_q, d_rho);
#endif
    }
      //merge d_rho
    merge_d_rho<<<1,nx>>>(d_rho);


    cudaMemcpy(rho, d_rho, sizeof(float_p)*(nx), cudaMemcpyDeviceToHost);


    //free(part_rho);    
    //B.C.
    rho[nx]  = rho[0];
    rho[nx+1]= rho[1];
    
    // printf("\n\n");

    // for(int i=0;i<nx+2;i++){
    //    printf("%d %f\n",i,rho[i]);
    // }
    //  printf("\n\n");

    // exit(1);
  }
  /* end device_move_pcles() */

  //push velocity backward half timestep 
  void device_init_vp_exp(float_p dt)
  {
    cudaPrintfInit();
    int tBlockSize,tGridSize;
    if(tot_Np==2){
      tGridSize  = 1;
      tBlockSize = 2;
    }else{
      tGridSize  = GRID_SIZE;
      tBlockSize = BLOCK_SIZE;
    }

    init_vp_exp<<<tGridSize,tBlockSize>>>(pcles, d_acNpcls, d_xx,d_q_m, d_E,dt);
    
    cudaPrintfDisplay();
    cudaPrintfEnd();    
    
  }
  /* end device_init_vp_exp() */

  //find moments of each species
  void device_moments(float_p *moments,float_p *ki_en)
  {
     int numBlocks=16, numThreads=64;
     //getNBT_re(tot_Np, numBlocks, numThreads);  
     

     //for toal momentum
     cudaMalloc((void **)&d_mv, sizeof(float_p)*numBlocks*h_Nsp); 

     for(int isp=0; isp<h_Nsp; isp++){
       output_mv<<<numBlocks, numThreads,
    	 numThreads * sizeof(float_p)>>>(pcles.d_v_n,isp,d_mv,numThreads); //careful
     }


     float_p h_part_sum[numBlocks*h_Nsp];
     cudaMemcpy(&h_part_sum, d_mv, sizeof(float_p) * numBlocks*h_Nsp,
     		cudaMemcpyDeviceToHost);
     float_p final_sum = 0;
     for(int i = 0; i < numBlocks*h_Nsp; i++) {
       final_sum += h_part_sum[i];
     }
     moments[1] = final_sum;


    //sum up kineitc energy
     cudaMalloc((void **)&d_v2, sizeof(float_p)*numBlocks*h_Nsp); 
     moments[2] = 0;
     for(int isp=0; isp<h_Nsp; isp++){
       output_v2<<<numBlocks, numThreads,
	 numThreads * sizeof(float_p)>>>(pcles.d_v_n,pcles.d_v_n,isp,d_v2,numThreads); //careful
       cudaMemcpy(&h_part_sum[isp*numBlocks], d_v2+isp*numBlocks, sizeof(float_p) * numBlocks, cudaMemcpyDeviceToHost);
       final_sum = 0;
       for(int i = isp*numBlocks; i < (isp+1)*numBlocks; i++) {
	 final_sum += h_part_sum[i];
       }
       ki_en[isp] = 0.5*h_mass[isp]*final_sum;
       moments[2]+= ki_en[isp];
     }

  }




/*----- implicit functions ------*/

void device_move_pcles_acc_cn(float_p dt, float_p j[], int nx, int itime,float_p dE[],bool use_sp)
  {
    //cudaPrintfInit();

    int tBlockSize,tGridSize;
    if(tot_Np==2){
      tGridSize  = 1;
      tBlockSize = 2;
    }else{
      tGridSize  = GRID_SIZE;
      tBlockSize = G_WARP*WARP_SIZE; //BLOCK_SIZE;
    }

    if(use_sp){
      for(int isp=0; isp<h_Nsp; isp++){
	move_pcles_acc_cn_sp<<<tGridSize,tBlockSize,G_WARP*nx*sizeof(float)>>>(pcles, isp, d_E,dt, d_j);
      }
    }else{
      for(int isp=0; isp<h_Nsp; isp++){
      //move_pcles_cn<<<tGridSize,tBlockSize>>>(pcles, isp, d_E,dt, d_j);
	move_pcles_acc_cn<<<tGridSize,tBlockSize,G_WARP*nx*sizeof(float_p)>>>(pcles, isp, d_E,dt, d_j);
      //move_pcles_acc_cn2<<<tGridSize,tBlockSize,G_WARP*nx*sizeof(float_p)>>>(pcles, isp, d_E,dt, d_j);
      //cudaPrintfDisplay();
      }
    }

    //cudaPrintfEnd();    
    cudaDeviceSynchronize();
    //merge d_j
    /*
    merge_d_rho<<<1,nx>>>(d_j);

    cudaMemcpy(j, d_j, sizeof(float_p)*(nx), cudaMemcpyDeviceToHost);

    for(int i=0; i<nx; i++){
      j[i] = j[i]/(h_hx*dt);
    }
    */
    size_t memSize = sizeof(float_p)*(nx)*h_Nsp*GRID_SIZE*2;

    float_p *h_j = (float_p *)malloc(memSize);
    cudaMemcpy(h_j, d_j, memSize, cudaMemcpyDeviceToHost);
    for(int ix=0; ix<nx; ix++){
      j[ix] = _fconst(0.0);
    }
    //merge
    float_p a,b,c;
    for(int ix=0; ix<nx; ix++){
      float_p e0=_fconst(0.0),e1=_fconst(0.0);
      for(int ib=1; ib<GRID_SIZE; ib++){
	h_j[ix] += h_j[ix+ib*nx];
	h_j[ix+GRID_SIZE*nx] +=  h_j[ix+GRID_SIZE*nx+ib*nx];
	h_j[ix+MS] += h_j[ix+MS+ib*nx];
	h_j[ix+MS+GRID_SIZE*nx] +=  h_j[ix+MS+GRID_SIZE*nx+ib*nx];
	//	printf("%d %d %24.16e  %24.16e \n",ix, ib,h_j[ix+ib*nx],h_j[ix+ib*nx+MS]);

	/*
	c = h_j[ix+ib*nx]; 
	if(_abs(c)>_abs(h_j[ix])){
	  a = c;
	  c = h_j[ix]; 
	}else{
	  a = h_j[ix];
	}
	b = e0 + c;
	h_j[ix] = a+b;
	e0 = b-(h_j[ix]-a);

	c = h_j[ix+GRID_SIZE*nx+ib*nx]; 
	if(_abs(c)>_abs(h_j[ix+GRID_SIZE*nx])){
	  a = c;
	  c = h_j[ix+GRID_SIZE*nx]; 
	}else{	  
	  a = h_j[ix+GRID_SIZE*nx];
	}
	b = e1 + c;
	h_j[ix+GRID_SIZE*nx] = a+b;
	e1 = b-(h_j[ix+GRID_SIZE*nx]-a);
	*/
      }
    }
    for(int ix=0; ix<nx; ix++){
      for(int is=0; is<h_Nsp; is++){
	j[ix] += h_j[ix+is*GRID_SIZE*nx]*h_charge[is];
	//printf("%d %d %24.16e  %24.16e \n",is, ix,h_j[ix+is*GRID_SIZE],j[ix]);
      }
      for(int is=0; is<h_Nsp; is++){
	j[ix] += h_j[MS+ix+is*GRID_SIZE*nx]*h_charge[is];
	//printf("%d %d %24.16e  %24.16e \n",is, ix,h_j[MS+ix+is*GRID_SIZE],j[ix]);
      }

      j[ix] = j[ix]/(h_hx*dt);
    }
    //B.C.
    j[nx]  = j[0];
    j[nx+1]= j[1];
    
    checkCUDAError("move(acc_cn) particles");
    free(h_j);
    //    if(j[20]!=j[20]){
      //if(itime > 713){
    //for(int i=0;i<nx+2;i++){
    // printf("%d %24.16e %24.16e \n",i,dE[i],j[i]);
    //}
    //printf("#%d setupNonlinearFunction after mover.\n",itime);	
    //printf("\n\n");
    // exit(1);
    //    }

  }
  /* end device_move_pcles_acc_cn() */


  void device_update_tnp_pcles()
  {
    int tBlockSize,tGridSize;
    if(tot_Np==2){
      tGridSize  = 1;
      tBlockSize = 2;
    }else{
      tGridSize  = 1; //maxBlocks; //GRID_SIZE;
      tBlockSize = 64; //maxThreads; //BLOCK_SIZE;
    }
    

#ifdef _SORT
    for(int isp=0; isp<h_Nsp; isp++){
      sortv_particles<<<tGridSize,tBlockSize>>>(pcles,isp);
    }
#else
    update_tnp_pcles<<<tGridSize,tBlockSize>>>(pcles);
#endif
  }
  /* end device_update_tnp_pcles() */
}
