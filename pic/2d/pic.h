extern "C" {
#define DIM_PIC 1
#define MAXNSP_cpu (4)

  class Particles_cpu{
  public:
    Particles_cpu(const int     Nsp, 
		  const int     totNp,
		  const int     Nx,
		  const float_p xx0,
		  const float_p xf0,
		  const float_p ihx,
		  const float_p hx,
		  const float_p LL,
		  const float_p eps,
		  const float_p etola,
		  const float_p etolr,

		  const float_p *c_vth,
		  const float_p *c_Npcles,
		  const float_p *c_acNpcles,
		  const float_p *q,
		  const float_p *m,
		  const float_p *q_m
		  );
    ~Particles_cpu();
    
    float_p *x_n;  //xn
    float_p *x_np; //xnp
    float_p *v_n;  //xn
    float_p *v_np; //xnp
    int     *ci_n; //cin
    int     *ci_np;//cinp


    const float_p  c_q_m[MAXNSP_cpu];
    const float_p    c_m[MAXNSP_cpu];
    const float_p    c_q[MAXNSP_cpu];
    const int  c_acNpcls[MAXNSP_cpu+1];
    const int    c_Npcls[MAXNSP_cpu];
    const float_p  c_vth[MAXNSP_cpu];
    
    const float_p etolr;
    const float_p etola;
    const float_p eps;
    const float_p LL;
    const float_p hx;
    const float_p ihx;
    const float_p xf0;
    const float_p xx0;
    const int  Nx;
    const int  totNp;
    const int Nsp;    

    void cpu_acc_cn_incell_offset(int isp, float_p x_n,float_p v_n,int i_n, float_p dt, float_p *h_E, float_p &x_np, float_p &v_np, int &i_np, float_p &x_nh, int &i_nh,float_p &dt_sub,int pi);    

  };


  void device_load_pcles(float_p v0[],float_p vth[],int numParticles[], int n_sp, float_p xx[], int nx);



}
