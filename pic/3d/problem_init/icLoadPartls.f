
      seed = 1
 
      do is=1,size(spcs)

         vx2 = 0d0
         vt2 = 0d0 
         v_tot = 0d0

         ipc0 = 0

         do k = 1,nzg
            do j = 1,nyg
               do i = 1,nxg
                  ii = i + nxg*(j-1) + nxg*nyg*(k-1)

c$$$!$OMP PARALLEL DEFAULT(SHARED) private(ipc,ip,ip_ng,rx,ry,rz
c$$$!$OMP.                                ,rx1,rx2,rx3,xp,yp,zp)
c$$$!$OMP DO
c$$$!$OMP.REDUCTION(+:v_tot,vx2,vt2)
             do ipc = 1, npc_int(ii,is)

               !note: SEED in parallel!
               do ip_ng=1,_Npg
                 call unif_rng(rx(ip_ng))
                 call unif_rng(ry(ip_ng))
                 call unif_rng(rz(ip_ng))
                 call gauss2(rx1(ip_ng))
                 call gauss2(rx2(ip_ng))
                 call gauss2(rx3(ip_ng))
               enddo

               ip = ipc0+ipc

               spcs(is)%pcles(ip)%ijk_np(:) = ii

               xp = hx*rx
               spcs(is)%pcles(ip)%x_np(:,1) = xp
                   
               yp = hy*ry
               spcs(is)%pcles(ip)%x_np(:,2) = yp

               zp = hz*rz
               spcs(is)%pcles(ip)%x_np(:,3) = zp

               spcs(is)%pcles(ip)%w_np(:)   = 1d0
                  
               spcs(is)%pcles(ip)%x_n   = spcs(is)%pcles(ip)%x_np
               spcs(is)%pcles(ip)%ijk_n = spcs(is)%pcles(ip)%ijk_np
               spcs(is)%pcles(ip)%w_n   = spcs(is)%pcles(ip)%w_np

              !Initialize particle velocity (Maxwellian)            
              !This code requires similar mods in find_pcle_drift_velocity
              if(ipc<mnpc_int(ii,is)) then
                 spcs(is)%pcles(ip)%v_n(:,1) = v_thx(is)*rx1 
                 spcs(is)%pcles(ip)%v_n(:,2) = v_thy(is)*rx2 
                 spcs(is)%pcles(ip)%v_n(:,3) = v_thz(is)*rx3 
              else
                 spcs(is)%pcles(ip)%v_n(:,1) = v_thx(is)*rx1 +v0_x(is)
                 spcs(is)%pcles(ip)%v_n(:,2) = v_thy(is)*rx2 +v0_y(is)
                 spcs(is)%pcles(ip)%v_n(:,3) = v_thz(is)*rx3 +v0_z(is)
              end if

              spcs(is)%pcles(ip)%v_np = spcs(is)%pcles(ip)%v_n

               v_tot(is,:) = v_tot(is,:) + sum(spcs(is)%pcles(ip)%v_n,1)

               vx2 = vx2 + sum(spcs(is)%pcles(ip)%v_n(:,1)
     .                        *spcs(is)%pcles(ip)%v_n(:,1))
               vt2 = vt2 + sum(spcs(is)%pcles(ip)%v_n(:,2)
     .                        *spcs(is)%pcles(ip)%v_n(:,2))
     .                   + sum(spcs(is)%pcles(ip)%v_n(:,3)
     .                        *spcs(is)%pcles(ip)%v_n(:,3))

             end do             !particles
             ipc0 = ipc0 + npc_int(ii,is)
           end do               !cell x
         end do                 !cell y
       end do                   !cell z

c$$$       print *,v_tot,vx2,vt2

      end do                    !species

c$$$       stop


      
