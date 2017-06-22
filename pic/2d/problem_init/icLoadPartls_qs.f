
       do is=1,size(spcs)
         allocate(r4(dimt,npc_max(is))) 

         vx2 = 0d0
         vt2 = 0d0 
         v_tot = 0d0

         ipc0 = 0

         !Initialize particle velocity (anisotropic Maxwellian)   
         vt(1) = v_thx(is)
         vt(2) = v_thy(is)
         vt(3) = v_thz(is)

         do j = 1,nyg
           do i = 1,nxg
             ii = i + nxg*(j-1) 

             call HamSeq(npc_int(ii,is),r4,dimt,npc_scan(ii,is)
     $                  ,rank=my_rank,nproc=np)

c$$$!$OMP PARALLEL DEFAULT(SHARED) private(ip,ipc,ipl,ip_ng,xp,yp
c$$$!$OMP.    ,rx,rx1,rx2,rx3,signx,signy,signz,ixyz,ipx,ipy,ipz)
c$$$!$OMP DO
c$$$!$OMP.REDUCTION(+:v_tot,vx2,vt2)
             do ipc = 1, npc_int(ii,is)
               rx = r4(1,ipc) 
               where (rx==0d0) rx = rx + 1d-16 !end with 1
               where (rx==1d0) rx = rx - 1d-16
               xp = hx*rx 

               rx = r4(2,ipc)
               where (rx==0d0) rx = rx + 1d-16 !end with 1
               where (rx==1d0) rx = rx - 1d-16
               yp = hy*rx 

               rx1= dinvnorm((r4(3,ipc)+1d0)*0.5d0)  
               rx2= dinvnorm((r4(4,ipc)+1d0)*0.5d0)
               rx3= dinvnorm((r4(5,ipc)+1d0)*0.5d0)

               !Fill particle properties
               ip = ipc0+ipc
                  
               do ip_ng=1,8/_Npg !Group particles

                 !Particle index
                 ipl = (ip-1)*8/_Npg + ip_ng

                 !Fill particle properties
                 spcs(is)%pcles(ipl)%ijk_np(:) = ii
                 spcs(is)%pcles(ipl)%x_np(:,1) = xp
                 spcs(is)%pcles(ipl)%x_np(:,2) = yp
                 spcs(is)%pcles(ipl)%x_np(:,3) = 0d0
                 spcs(is)%pcles(ipl)%w_np(:)   = 1d0

                 spcs(is)%pcles(ipl)%x_n   = spcs(is)%pcles(ipl)%x_np
                 spcs(is)%pcles(ipl)%ijk_n = spcs(is)%pcles(ipl)%ijk_np
                 spcs(is)%pcles(ipl)%w_n   = spcs(is)%pcles(ipl)%w_np

               enddo

               !Initialize particle velocity (anisotropic Maxwellian)   
               signx = -1 
               signy = -1
               signz = -1
               ixyz  = 0
               do ipx = 0,1
                 signx = -signx
                 do ipy = 0,1
                   signy = -signy
                   do ipz = 0,1
                     signz = -signz

                     ipl = (ip-1)*8/_Npg + 1 + ixyz/_Npg

                     ip_ng = mod(ixyz,_Npg) + 1

                     if (ipc.le.mnpc_int(ii,is)) then
                       spcs(is)%pcles(ipl)%v_n(ip_ng,1)
     .                    = signx*vt(1)*rx1(ip_ng)
                       spcs(is)%pcles(ipl)%v_n(ip_ng,2)
     .                    = signy*vt(2)*rx2(ip_ng)
                       spcs(is)%pcles(ipl)%v_n(ip_ng,3)
     .                    = signz*vt(3)*rx3(ip_ng)
                       spcs(is)%pcles(ipl)%w_n = w_bp(is)
                       spcs(is)%pcles(ipl)%w_np= w_bp(is)
                     else
                       spcs(is)%pcles(ipl)%v_n(ip_ng,1)
     .                    = signx*vt(1)*rx1(ip_ng) +v_0x(is)
                       spcs(is)%pcles(ipl)%v_n(ip_ng,2)
     .                    = signy*vt(2)*rx2(ip_ng) +v_0y(is)
                       spcs(is)%pcles(ipl)%v_n(ip_ng,3)
     .                    = signz*vt(3)*rx3(ip_ng) +v_0z(is)
                     endif

                     spcs(is)%pcles(ipl)%v_np(ip_ng,:)
     .                    = spcs(is)%pcles(ipl)%v_n(ip_ng,:)

                     ixyz = ixyz + 1
                   end do
                 end do
               end do

               do ip_ng=1,8/_Npg  !Group particles

                 !Particle index
                 ipl = (ip-1)*8/_Npg + ip_ng

                 v_tot(is,:) = v_tot(is,:)
     .                       + sum(spcs(is)%pcles(ipl)%v_n,1)

                 vx2 = vx2 + sum(spcs(is)%pcles(ipl)%v_n(:,1)
     .                          *spcs(is)%pcles(ipl)%v_n(:,1))
                 vt2 = vt2 + sum(spcs(is)%pcles(ipl)%v_n(:,2)
     .                          *spcs(is)%pcles(ipl)%v_n(:,2))
     .                     + sum(spcs(is)%pcles(ipl)%v_n(:,3)
     .                          *spcs(is)%pcles(ipl)%v_n(:,3))
               enddo

             end do             !particles
c$$$!$OMP END DO
c$$$!$OMP END PARALLEL
             ipc0 = ipc0 + npc_int(ii,is)
           end do               !cell
         end do                 !cell y

cc         write (*,*) v_thx,v_thy,v_thz,v0_x,v0_y,v0_z
cc         print *,is,v_tot,vx2,vt2

         deallocate(r4)
       end do !species

c$$$       stop
