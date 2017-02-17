void  cpu_acc_cn_incell_offset(int isp, float_p x_n,float_p v_n,int i_n, float_p dt, float_p *h_E, float_p &x_np, float_p &v_np, int &i_np, float_p &x_nh, int &i_nh,float_p &dt_sub,int pi)
{
  dt_sub = dt;
  float_p E0 = h_E[i_n-1];
  float_p E1 = h_E[i_n  ];

  float_p N = x_n;
  float_p D = hx;
  float_p R = ihx;
  float_p x = N*R;                           
  float_p r = N - x*D;
  x += r*R;

  float_p a_n;

}
