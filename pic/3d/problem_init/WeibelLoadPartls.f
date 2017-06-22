
       seed = 1
       v_tot = 0d0

       do is=1,size(spcs)
         pcles => spcs(is)%pcles
         vx2 = 0d0
         vt2 = 0d0 

         do ip=1,size(pcles)       !for every mpi-proc

           !note: SEED in parallel!
           call unif_rng(rx)
           call unif_rng(ry)
           call gauss2(rx1)
           call gauss2(rx2)
           call gauss2(rx3)

           !Uniform particle allocation
           cii = mod(ip-1,nxg*nyg) 
           i  = mod(cii, nxg) + 1
           j = (cii-1)/(nxg) + 1

           xp = hx*rx;
           pcles(ip)%x_np(1) = xp
           pcles(ip)%i_np    = i
                   
           yp = hy*ry
           pcles(ip)%x_np(2) = yp
           pcles(ip)%j_np    = j

           pcles(ip)%x_np(3) = 0d0

           pcles(ip)%w  = 1d0

           pcles(ip)%x_n = pcles(ip)%x_np
           pcles(ip)%i_n = pcles(ip)%i_np
           pcles(ip)%j_n = pcles(ip)%j_np

           !Initialize particle velocity (anisotropic Maxwellian)   
           vt(1) = v_thx(is)
           vt(2) = v_thy(is)
           vt(3) = v_thz(is)
           v0    = 0d0                    

           pcles(ip)%v_n (1) = vt(1)*rx1
           pcles(ip)%v_n (2) = vt(2)*rx2
           pcles(ip)%v_n (3) = vt(3)*rx3

           pcles(ip)%v_np = pcles(ip)%v_n

           v_tot(is,:) = v_tot(is,:) + pcles(ip)%v_n 

           vx2 = vx2 + pcles(ip)%v_n(1)*pcles(ip)%v_n(1)
           vt2 = vt2 + pcles(ip)%v_n(2)*pcles(ip)%v_n(2)
     .               + pcles(ip)%v_n(3)*pcles(ip)%v_n(3)

         end do                 !particles
       end do                   !species

