!     only with quiet_start

      isp = 1

      kx = 2d0*acos(-1d0)*dble(nh1)/Lx

      do ip=1,size(spcs(isp)%pcles)

         call xform_pcle_idx(spcs(isp)%pcles(ip)%ijk_n,nxg,nyg,
     .                                          i_np,j_np,k_np)

         !perturb n
         do k=1,_Npg
           xn = spcs(isp)%pcles(ip)%x_n(k,1)

           xp = xn + pxx(i_np(k))  - 0.5d0*hx
 
           kxp = kx*xp

           x_np(k) = xn + eps_pic/kx*cos(kxp)

           do while(x_np(k)<0d0)
             x_np(k) = x_np(k) + hx
             i_np(k) = i_np(k) - 1
           enddo 
           
           do while(x_np(k)>hx)
             x_np(k) = x_np(k) - hx
             i_np(k) = i_np(k) + 1
           end do

           i_np(k) = iand(i_np(k)-1,nxg-1)+1 !careful

           spcs(isp)%pcles(ip)%x_np(k,1) = x_np(k)

         enddo

c$$$         if(x_np==0 .or. x_np==hx ) then
c$$$           print *,isp,ip,x_np,spcs(isp)%pcles(ip)%x_np(1),eps_pic/kx*cos(kxp)
c$$$     $             ,spcs(isp)%pcles(ip)%v_np(1)
c$$$           stop
c$$$         end if         

         spcs(isp)%pcles(ip)%ijk_np = i_np + nxg*(j_np-1)
         spcs(isp)%pcles(ip)%x_np(:,1) = x_np

         spcs(isp)%pcles(ip)%x_n(:,1) = spcs(isp)%pcles(ip)%x_np(:,1)
         spcs(isp)%pcles(ip)%ijk_n =  spcs(isp)%pcles(ip)%ijk_np

      end do
