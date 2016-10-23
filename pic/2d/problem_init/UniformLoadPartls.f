
c$$$       seed = 1
c$$$c$$$       v_tot = 0d0
c$$$
c$$$       do is=1,size(spcs)
c$$$cc         pcles => spcs(is)%pcles
c$$$         vx2 = 0d0
c$$$         vt2 = 0d0 
c$$$
c$$$         numpcles = size(spcs(is)%pcles)*_Npg
c$$$
c$$$!$OMP PARALLEL DEFAULT(SHARED) private(ip,ip_ng,ipl,rx,ry
c$$$!$OMP.                        ,rx1,rx2,rx3,cii,xp,yp,vt)
c$$$!$OMP DO
c$$$c$$$!$OMP.REDUCTION(+:v_tot,vx2,vt2)
c$$$         do ip=1,size(spcs(is)%pcles)       !for every mpi-proc
c$$$
c$$$           do ip_ng=1,_Npg
c$$$             ipl = (ip-1)*_Npg + ip_ng
c$$$
c$$$             !Uniform particle allocation
c$$$             cii = mod(ipl-1,nxg*nyg) + 1
c$$$             spcs(is)%pcles(ip)%ijk_np(ip_ng) = cii
c$$$
c$$$c$$$             call unif_rng(rx(ip_ng))
c$$$c$$$             call unif_rng(ry(ip_ng))
c$$$             call gauss2(rx1(ip_ng))
c$$$             call gauss2(rx2(ip_ng))
c$$$             call gauss2(rx3(ip_ng))
c$$$
c$$$             !Only works in 1D
c$$$             xp(ip_ng) = mod(Lx/numpcles*(ipl-0.5),hx)
c$$$cc             xp(ip_ng) = hx/npc_int(cii,is)*((ipl-1)/(nxg*nyg)+0.5)
c$$$cc             write (*,*) hx,'position',xp(ip_ng)
c$$$           enddo
c$$$
c$$$c$$$           xp = hx*rx
c$$$           spcs(is)%pcles(ip)%x_np(:,1) = xp
c$$$                   
c$$$c$$$           yp = hy*ry
c$$$           yp = 0.5*hy
c$$$           spcs(is)%pcles(ip)%x_np(:,2) = yp
c$$$
c$$$           spcs(is)%pcles(ip)%x_np(:,3) = 0d0
c$$$
c$$$           spcs(is)%pcles(ip)%w_np(:)   = 1d0
c$$$
c$$$           spcs(is)%pcles(ip)%x_n   = spcs(is)%pcles(ip)%x_np
c$$$           spcs(is)%pcles(ip)%ijk_n = spcs(is)%pcles(ip)%ijk_np
c$$$           spcs(is)%pcles(ip)%w_n   = spcs(is)%pcles(ip)%w_np
c$$$
c$$$           !Initialize particle velocity (anisotropic Maxwellian)   
c$$$           vt(1) = v_thx(is)
c$$$           vt(2) = v_thy(is)
c$$$           vt(3) = v_thz(is)
c$$$
c$$$           spcs(is)%pcles(ip)%v_n(:,1) = vt(1)*rx1 + v0_x(is)
c$$$           spcs(is)%pcles(ip)%v_n(:,2) = vt(2)*rx2 + v0_y(is)
c$$$           spcs(is)%pcles(ip)%v_n(:,3) = vt(3)*rx3 + v0_z(is)
c$$$
c$$$           spcs(is)%pcles(ip)%v_np = spcs(is)%pcles(ip)%v_n
c$$$
c$$$c$$$           v_tot(is,:) = v_tot(is,:) + pcles(ip)%v_n 
c$$$c$$$
c$$$c$$$           vx2 = vx2 + pcles(ip)%v_n(1)*pcles(ip)%v_n(1)
c$$$c$$$           vt2 = vt2 + pcles(ip)%v_n(2)*pcles(ip)%v_n(2)
c$$$c$$$     .               + pcles(ip)%v_n(3)*pcles(ip)%v_n(3)
c$$$
c$$$         end do                 !particles
c$$$!$OMP END DO
c$$$!$OMP END PARALLEL
c$$$       end do                   !species

cc       v_tot = 0d0

       do is=1,size(spcs)
         vx2 = 0d0
         vt2 = 0d0 
         v_tot = 0d0

         vt(1) = v_thx(is)
         vt(2) = v_thy(is)
         vt(3) = v_thz(is)

         ipc0 = 0

         do j = 1,nyg
           do i = 1,nxg
             ii = i + nxg*(j-1) 

!$OMP PARALLEL DEFAULT(SHARED) private(ip,ipc,ip_ng,xp,yp
!$OMP.                                ,rx1,rx2,rx3,vt)
!$OMP DO
cc!$OMP.REDUCTION(+:v_tot,vx2,vt2)
             do ipc = 1, npc_int(ii,is)

               xp = hx/npc_int(ii,is)*(ipc-0.5)
               yp = 0.5*hy

               !Fill particle properties
               ip = ipc0+ipc

               !note: SEED in parallel!
               do ip_ng=1,_Npg
                 call gauss2(rx1(ip_ng))
                 call gauss2(rx2(ip_ng))
                 call gauss2(rx3(ip_ng))
               enddo

               spcs(is)%pcles(ip)%ijk_np(:) = ii

               spcs(is)%pcles(ip)%x_np(:,1) = xp
                   
               spcs(is)%pcles(ip)%x_np(:,2) = yp

               spcs(is)%pcles(ip)%x_np(:,3) = 0d0
                   
               spcs(is)%pcles(ip)%w_np(:)   = 1d0
                  
               spcs(is)%pcles(ip)%x_n   = spcs(is)%pcles(ip)%x_np
               spcs(is)%pcles(ip)%ijk_n = spcs(is)%pcles(ip)%ijk_np
               spcs(is)%pcles(ip)%w_n   = spcs(is)%pcles(ip)%w_np

               !Initialize particle velocity (anisotropic Maxwellian)   
               spcs(is)%pcles(ip)%v_n(:,1) = vt(1)*rx1 + v0_x(is)
               spcs(is)%pcles(ip)%v_n(:,2) = vt(2)*rx2 + v0_y(is)
               spcs(is)%pcles(ip)%v_n(:,3) = vt(3)*rx3 + v0_z(is)

               spcs(is)%pcles(ip)%v_np = spcs(is)%pcles(ip)%v_n

               v_tot(is,:) = v_tot(is,:) + sum(spcs(is)%pcles(ip)%v_n,1)

               vx2 = vx2 + sum(spcs(is)%pcles(ip)%v_n(:,1)
     .                        *spcs(is)%pcles(ip)%v_n(:,1))
               vt2 = vt2 + sum(spcs(is)%pcles(ip)%v_n(:,2)
     .                        *spcs(is)%pcles(ip)%v_n(:,2))
     .                   + sum(spcs(is)%pcles(ip)%v_n(:,3)
     .                        *spcs(is)%pcles(ip)%v_n(:,3))
c$$$               v_tot(is,:) = v_tot(is,:) + spcs(is)%pcles(ip)%v_n 
c$$$
c$$$               vx2 = vx2 + spcs(is)%pcles(ip)%v_n(1)*spcs(is)%pcles(ip)%v_n(1)
c$$$               vt2 = vt2 + spcs(is)%pcles(ip)%v_n(2)*spcs(is)%pcles(ip)%v_n(2)
c$$$     .                   + spcs(is)%pcles(ip)%v_n(3)*spcs(is)%pcles(ip)%v_n(3)

             end do             !particles
!$OMP END DO
!$OMP END PARALLEL
             ipc0 = ipc0 + npc_int(ii,is)
           end do               !cell
         end do                 !cell y
       end do !species

