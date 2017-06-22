
       seed = 1
       v_tot = 0d0

       do is=1,size(spcs)

         vx2 = 0d0
         vt2 = 0d0 

         i = 1
         xp = 0d0
         tot_np = size(spcs(is)%pcles)*_Npg
         ipb = int(tot_np*r_b)

!$OMP PARALLEL DEFAULT(SHARED) private(ip,ip_ng,rx,ry,rx1,rx2,rx3,cii
!$OMP.                        ,xp,yp,vt,v0)
!$OMP DO
c$$$!$OMP.REDUCTION(+:v_tot,vx2,vt2)
        do ip=1,size(spcs(is)%pcles)     !for every mpi-proc

           !note: SEED in parallel!
           do ip_ng=1,_Npg
             call unif_rng(rx(ip_ng))
             call unif_rng(ry(ip_ng))
             call gauss2(rx1(ip_ng))
             call gauss2(rx2(ip_ng))
             call gauss2(rx3(ip_ng))
           enddo

           !Uniform particle allocation
           cii = mod(ip-1,nxg*nyg) + 1
           spcs(is)%pcles(ip)%ijk_np = cii

           xp = hx*rx;
           spcs(is)%pcles(ip)%x_np(:,1) = xp
                
           yp = hy*ry
           spcs(is)%pcles(ip)%x_np(:,2) = yp

           spcs(is)%pcles(ip)%x_np(:,3) = 0d0
        
           spcs(is)%pcles(ip)%w_np(:)   = 1d0
                
           spcs(is)%pcles(ip)%x_n   = spcs(is)%pcles(ip)%x_np
           spcs(is)%pcles(ip)%ijk_n = spcs(is)%pcles(ip)%ijk_np
           spcs(is)%pcles(ip)%w_n   = spcs(is)%pcles(ip)%w_np

           !Initialize particle velocity (drifting Maxwellian)   
           if(is==i_b) then
c            call unif_rng(rx)
c            if(rx>r_b) then
             if(ip*_Npg>ipb) then
               vt(1) = v_thx(is+1)
               vt(2) = v_thy(is+1)
               vt(3) = v_thz(is+1)
               v0(1) = v_0x(is+1)
               v0(2) = v_0y(is+1)
               v0(3) = v_0z(is+1)
             else
               vt(1) = v_thx(is)
               vt(2) = v_thy(is)
               vt(3) = v_thz(is)
               v0(1) = v_0x (is)
               v0(2) = v_0y (is)
               v0(3) = v_0z (is)
             end if
           else
             vt(1) = v_thx(is)
             vt(2) = v_thy(is)
             vt(3) = v_thz(is)
             v0    = 0d0                    
           end if

           spcs(is)%pcles(ip)%v_n(:,1) = vt(1)*rx1 + v0(1)
           spcs(is)%pcles(ip)%v_n(:,2) = vt(2)*rx2 + v0(2)
           spcs(is)%pcles(ip)%v_n(:,3) = vt(3)*rx3 + v0(3)

           spcs(is)%pcles(ip)%v_np = spcs(is)%pcles(ip)%v_n

c$$$           v_tot(is,:) = v_tot(is,:) + spcs(is)%pcles(ip)%v_n 
c$$$
c$$$           vx2 = vx2 + spcs(is)%pcles(ip)%v_n(1)*spcs(is)%pcles(ip)%v_n(1)
c$$$           vt2 = vt2 + spcs(is)%pcles(ip)%v_n(2)*spcs(is)%pcles(ip)%v_n(2)
c$$$     .               + spcs(is)%pcles(ip)%v_n(3)*spcs(is)%pcles(ip)%v_n(3)

         end do                 !particles
!$OMP END DO
!$OMP END PARALLEL
       end do                   !species

