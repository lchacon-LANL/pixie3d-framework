
       do is=1,size(spcs)
         vx2 = 0d0
         vt2 = 0d0 
         v_tot = 0d0

         vt(1) = v_thx(is)
         vt(2) = v_thy(is)
         vt(3) = v_thz(is)

         ipc0 = 0

         do k = 1,nzg
           do j = 1,nyg
             do i = 1,nxg
               ii = i + nxg*(j-1) + nxg*nyg*(k-1)

!$OMP PARALLEL DEFAULT(SHARED) private(ip,ipc,ip_ng,xp,yp,zp
!$OMP.                                ,rx1,rx2,rx3)
!$OMP DO
cc!$OMP.REDUCTION(+:v_tot,vx2,vt2)
             do ipc = 1, npc_int(ii,is)

               xp = hx/npc_int(ii,is)*(ipc-0.5)
               yp = hy/npc_int(ii,is)*(ipc-0.5)
               zp = hz/npc_int(ii,is)*(ipc-0.5)

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

               spcs(is)%pcles(ip)%x_np(:,3) = zp
                   
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
           end do               !cell x
         end do                 !cell y
        end do                  !cell z
       end do !species

