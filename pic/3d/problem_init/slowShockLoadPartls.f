
       seed = 1
       v_tot = 0d0
       vx2 = 0d0
       vt2 = 0d0 

       do is=1,size(spcs)

         i = 1
         j = 1                !fix y-position
         xp = 0d0
         ipc0 = npc_int(i,is) + nip

         do ip=1,size(spcs(is)%pcles)    !for every mpi-proc

            if(ip > ipc0) then
cc           If(ip > nint(npc_acm(i)))then
cc                print *,i,gip, nint(npc_acm(i))*1
             i = i+1
             i = iand(i-1,nxg-1)+1
             ipc0 = ipc0 + npc_int(i,is)
           end if

           !note: SEED in parallel!
           do ip_ng=1,_Npg
             call unif_rng(rx(ip_ng))
             call unif_rng(ry(ip_ng))
             call gauss2(rx1(ip_ng))
             call gauss2(rx2(ip_ng))
             call gauss2(rx3(ip_ng))
           enddo

           if(ip.le.nip) then
cc                   xp = 0.00856d0*dtf*rx
             xp = 1d-16
             spcs(is)%pcles(ip)%x_np(:,1) = xp

             yp = hy*rx         !temporary
             spcs(is)%pcles(ip)%x_np(:,2) = yp

             i = 1
           else
             xp = hx*rx;
             spcs(is)%pcles(ip)%x_np(:,1) = xp

             yp = hy*rx         !temporary
             spcs(is)%pcles(ip)%x_np(:,2) = yp
           end if

           spcs(is)%pcles(ip)%w_np(:) = 1d0
                  
           spcs(is)%pcles(ip)%ijk_np = i + nxg*(j-1)

           spcs(is)%pcles(ip)%x_n   = spcs(is)%pcles(ip)%x_np
           spcs(is)%pcles(ip)%ijk_n = spcs(is)%pcles(ip)%ijk_np
           spcs(is)%pcles(ip)%w_n   = spcs(is)%pcles(ip)%w_np

           !Initialize particle velocity (drifting Maxwellian)            
           spcs(is)%pcles(ip)%v_n(:,1) = v_thx(is)*rx1 + v_0x(is)
           spcs(is)%pcles(ip)%v_n(:,2) = v_thy(is)*rx2 + v_0y(is)
           spcs(is)%pcles(ip)%v_n(:,3) = v_thz(is)*rx3 + v_0z(is)

           spcs(is)%pcles(ip)%v_np = spcs(is)%pcles(ip)%v_n

           v_tot(is,:) = v_tot(is,:) + sum(spcs(is)%pcles(ip)%v_n,1)

           vx2 = vx2 + sum(spcs(is)%pcles(ip)%v_n(:,1)
     .                    *spcs(is)%pcles(ip)%v_n(:,1))
           vt2 = vt2 + sum(spcs(is)%pcles(ip)%v_n(:,2)
     .                    *spcs(is)%pcles(ip)%v_n(:,2))
     .               + sum(spcs(is)%pcles(ip)%v_n(:,3)
     .                    *spcs(is)%pcles(ip)%v_n(:,3))

         end do                 !particles
       end do                   !species
