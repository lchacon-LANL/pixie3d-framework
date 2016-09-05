
       seed = 1
c$$$       v_tot = 0d0

       do is=1,size(spcs)
cc         pcles => spcs(is)%pcles
         vx2 = 0d0
         vt2 = 0d0 

         numpcles = size(spcs(is)%pcles)*_Npg

!$OMP PARALLEL DEFAULT(SHARED) private(ip,ip_ng,ipl,rx,ry
!$OMP.                        ,rx1,rx2,rx3,cii,xp,yp,vt)
!$OMP DO
c$$$!$OMP.REDUCTION(+:v_tot,vx2,vt2)
         do ip=1,size(spcs(is)%pcles)       !for every mpi-proc

           do ip_ng=1,_Npg
             ipl = (ip-1)*_Npg + ip_ng

             !Uniform particle allocation
             cii = mod(ipl-1,nxg*nyg) + 1
             spcs(is)%pcles(ip)%ijk_np(ip_ng) = cii

c$$$             call unif_rng(rx(ip_ng))
c$$$             call unif_rng(ry(ip_ng))
             call gauss2(rx1(ip_ng))
             call gauss2(rx2(ip_ng))
             call gauss2(rx3(ip_ng))

             !Only works in 1D
             xp(ip_ng) = mod(Lx/numpcles*(ipl-0.5),hx)
cc             xp(ip_ng) = hx/npc_int(cii,is)*((ipl-1)/(nxg*nyg)+0.5)
cc             write (*,*) hx,'position',xp(ip_ng)
           enddo

c$$$           xp = hx*rx
           spcs(is)%pcles(ip)%x_np(:,1) = xp
                   
c$$$           yp = hy*ry
           yp = 0.5*hy
           spcs(is)%pcles(ip)%x_np(:,2) = yp

           spcs(is)%pcles(ip)%x_np(:,3) = 0d0

           spcs(is)%pcles(ip)%w_np(:)   = 1d0

           spcs(is)%pcles(ip)%x_n   = spcs(is)%pcles(ip)%x_np
           spcs(is)%pcles(ip)%ijk_n = spcs(is)%pcles(ip)%ijk_np
           spcs(is)%pcles(ip)%w_n   = spcs(is)%pcles(ip)%w_np

           !Initialize particle velocity (anisotropic Maxwellian)   
           vt(1) = v_thx(is)
           vt(2) = v_thy(is)
           vt(3) = v_thz(is)

           spcs(is)%pcles(ip)%v_n(:,1) = vt(1)*rx1 + v0_x(is)
           spcs(is)%pcles(ip)%v_n(:,2) = vt(2)*rx2 + v0_y(is)
           spcs(is)%pcles(ip)%v_n(:,3) = vt(3)*rx3 + v0_z(is)

           spcs(is)%pcles(ip)%v_np = spcs(is)%pcles(ip)%v_n

c$$$           v_tot(is,:) = v_tot(is,:) + pcles(ip)%v_n 
c$$$
c$$$           vx2 = vx2 + pcles(ip)%v_n(1)*pcles(ip)%v_n(1)
c$$$           vt2 = vt2 + pcles(ip)%v_n(2)*pcles(ip)%v_n(2)
c$$$     .               + pcles(ip)%v_n(3)*pcles(ip)%v_n(3)

         end do                 !particles
!$OMP END DO
!$OMP END PARALLEL
       end do                   !species

