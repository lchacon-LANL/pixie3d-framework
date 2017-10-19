
       do is=1,size(spcs)
         allocate(r4(dimt,npc_max(is))) 

         vx2 = 0d0
         vt2 = 0d0 
         v_tot = 0d0

         ipc0 = 0

         do k = 1,nzg
           do j = 1,nyg
             do i = 1,nxg
               ii = i + nxg*(j-1) + nxg*nyg*(k-1)

               call HamSeq(npc_int(ii,is),r4,dimt,npc_scan(ii,is))

c$$$!$OMP PARALLEL DEFAULT(SHARED) private(ip,ipc,ipl,ip_ng,xp,yp,zp
c$$$!$OMP.    ,rx,rx1,rx2,rx3,vpar,vperp,signx,signy,signz,ixyz,ipx
c$$$!$OMP.    ,ipy,ipz)
c$$$!$OMP DO
c$$$!$OMP.REDUCTION(+:v_tot,vx2,vt2)
             do ipc = 1, npc_int(ii,is)

!               print '(6F16.8)',r4(1,ipc),r4(2,ipc),r4(3,ipc),r4(4,ipc),
!     .               r4(5,ipc),r4(6,ipc)

               rx = r4(1,ipc) 
               where (rx==0d0) rx = rx + 1d-16 !end with 1
               where (rx==1d0) rx = rx - 1d-16
               xp = hx*rx 

               rx = r4(2,ipc)
               where (rx==0d0) rx = rx + 1d-16 !end with 1
               where (rx==1d0) rx = rx - 1d-16
               yp = hy*rx 

               rx = r4(3,ipc)
               where (rx==0d0) rx = rx + 1d-16 !end with 1
               where (rx==1d0) rx = rx - 1d-16
               zp = hz*rx

               rx1= dinvnorm((r4(4,ipc)+1d0)*0.5d0)  
               rx2= dinvnorm((r4(5,ipc)+1d0)*0.5d0)
               rx3= dinvnorm((r4(6,ipc)+1d0)*0.5d0)

               !Particle group index
               ip = ipc0+ipc

               do ip_ng=1,8/_Npg  !Group particles

                 !Particle index
                 ipl = (ip-1)*8/_Npg + ip_ng

                 !Fill particle properties
                 spcs(is)%pcles(ipl)%ijk_np(:) = ii
                 spcs(is)%pcles(ipl)%x_np(:,1) = xp
                 spcs(is)%pcles(ipl)%x_np(:,2) = yp
                 spcs(is)%pcles(ipl)%x_np(:,3) = zp
                 spcs(is)%pcles(ipl)%w_np(:)   = 1d0

                 spcs(is)%pcles(ipl)%x_n   = spcs(is)%pcles(ipl)%x_np
                 spcs(is)%pcles(ipl)%ijk_n = spcs(is)%pcles(ipl)%ijk_np
                 spcs(is)%pcles(ipl)%w_n   = spcs(is)%pcles(ipl)%w_np

               enddo

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

!                     vpar(ip_ng) = signx*v_thx(is)*rx1(ip_ng)! + v0_x(is) ! Assume x is //-dirn.
!                    vperp(ip_ng) = signy*v_thy(is)*rx2(ip_ng)! + v0_y(is) !        y is perp-dirn.
                     vpar(ip_ng) = v_thx(is)*rx1(ip_ng)! + v0_x(is) ! Assume x is //-dirn.
                    vperp(ip_ng) = v_thy(is)*rx2(ip_ng)! + v0_y(is) !        y is perp-dirn.


!                    print*,'theta=',theta
!                    print*,'vpar=',vpar
!                    print*,'vperp=',vperp

                     spcs(is)%pcles(ipl)%v_n(ip_ng,1)
     .             = signx*(vpar(ip_ng)*cos(theta*pi/180.)
     .                    -vperp(ip_ng)*sin(theta*pi/180.))

! (vpar(ip_ng)*B0(1)-vperp(ip_ng)*B0(2))
!     .                       /sqrt(B0(1)*B0(1)+B0(2)*B0(2))

                     spcs(is)%pcles(ipl)%v_n(ip_ng,2)
     .             = signy*(vpar(ip_ng)*sin(theta*pi/180.)
     .                    +vperp(ip_ng)*cos(theta*pi/180.))
!     .                    = (vpar(ip_ng)*B0(2)+vperp(ip_ng)*B0(1))
!     .                       /sqrt(B0(1)*B0(1)+B0(2)*B0(2))
!                     print*,'vx=',spcs(is)%pcles(ipl)%v_n(:,1)
!                     print*,'vy=',spcs(is)%pcles(ipl)%v_n(:,2)


                     spcs(is)%pcles(ipl)%v_n(ip_ng,3)
     .                    = signz*v_thz(is)*rx3(ip_ng) + v0_z(is)

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

             end do             !Particle groups
c$$$!$OMP END DO
c$$$!$OMP END PARALLEL
             ipc0 = ipc0 + npc_int(ii,is)             
           end do               !cell x
          end do                !cell y
         end do                 !cell z
c$$$         print *,v_thx(is),v_thy(is),v_thz(is)
c$$$     .          ,v0_x(is),v0_y(is),v0_z(is)
c$$$         print *,is,v_tot,vx2,vt2

         deallocate(r4)
       end do                   !species

c$$$       stop
