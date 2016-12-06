
       seed = 1
       v_tot = 0d0

       tot_np = size(spcs(1)%pcles)
cc       gtot_np = tot_np*nproc

       if(mod(tot_np,nnl*8).ne.0) then
         print *, "adjust Np (", tot_np*nproc
     $         ,") to be a multiple of nxg*nyg*nzg*8 (",nng*8,")"
c$$$       else
c$$$         npc = tot_np/(nnl*8) !8 samples for each Gaussian random number
       end if

       allocate(r4(dimt,npc_max)) 

       do is=1,size(spcs)
         pcles => spcs(is)%pcles
         vx2 = 0d0
         vt2 = 0d0 

         do j = 1,nyg
           do i = 1,nxg
             ii = i + nxg*(j-1) 

             ipc0 = nint(npc_acm(ii-1))
             call HamSeq(npc_int(ii),r4,dimt,ipc0*nproc,my_rank_pic
     $            ,nproc)

             do ipc = 1, npc_int(ii)
               !without worries about parallel run #
               rx = r4(1,ipc) 
               if(rx==0d0) rx = rx + 1d-16 !end with 1
               if(rx==1d0) rx = rx - 1d-16
               xp = hx*rx 

               rx = r4(2,ipc)
               if(rx==0d0) rx = rx + 1d-16 !end with 1
               if(rx==1d0) rx = rx - 1d-16
               yp = hy*rx 

               rx1= dinvnorm((r4(3,ipc)+1d0)*0.5d0)  
               rx2= dinvnorm((r4(4,ipc)+1d0)*0.5d0)
               rx3= dinvnorm((r4(5,ipc)+1d0)*0.5d0)

               if((rx1.ne.rx1) .or. (rx2.ne.rx2) .or. (rx3.ne.rx3))then
                 print *, "not a number",rx1,rx2,rx3
                 stop
               end if

               !Fill particle properties
               ip = ipc0+ipc
               ip = (ip-1)*8 + 1

               pcles(ip:ip+7)%x_np(1) = xp
               pcles(ip:ip+7)%i_np    = i

               pcles(ip:ip+7)%x_np(2) = yp
               pcles(ip:ip+7)%j_np    = j

               pcles(ip:ip+7)%x_np(3) = 0d0

               pcles(ip:ip+7)%w       = 1d0

               if(xp<=0d0 .or. xp>=hx .or. yp<=0d0 .or. yp>=hy) then
                 print *,"xp or yp is out of cell: ",xp,yp,hx,hy
                 print *,isp,ip 
                 stop
               end if

               !Initialize particle velocity (anisotropic Maxwellian)   
               vt(1) = v_thx(is)
               vt(2) = v_thy(is)
               vt(3) = v_thz(is)
               v0    = 0d0                    

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
                     pcles(ip+ixyz)%v_n (1) = signx*vt(1)*rx1
                     pcles(ip+ixyz)%v_n (2) = signy*vt(2)*rx2
                     pcles(ip+ixyz)%v_n (3) = signz*vt(3)*rx3

                     pcles(ip+ixyz)%x_n = pcles(ip+ixyz)%x_np
                     pcles(ip+ixyz)%i_n = pcles(ip+ixyz)%i_np
                     pcles(ip+ixyz)%j_n = pcles(ip+ixyz)%j_np

                     pcles(ip+ixyz)%v_np = pcles(ip+ixyz)%v_n

                     ixyz = ixyz + 1
                   end do
                 end do
               end do

               v_tot(is,:) = v_tot(is,:) + pcles(ip)%v_n 

               vx2 = vx2 + pcles(ip)%v_n(1)*pcles(ip)%v_n(1)
               vt2 = vt2 + pcles(ip)%v_n(2)*pcles(ip)%v_n(2)
     .                   + pcles(ip)%v_n(3)*pcles(ip)%v_n(3)

             end do             !particles
           end do               !cell
         end do                 !cell y

       end do !species

       deallocate(r4)
