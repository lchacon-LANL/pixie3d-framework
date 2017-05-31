      kx = 2d0*pi*nh1/Lx
      ky =     pi    /Ly

c$$$      if (bcPER(1)) kx = 2d0*kx
c$$$      if (bcPER(2)) ky = 2d0*ky

      if (eps_pic /= 0d0) then
        if (mass_matrix_solve) then

          !Perturb the momenta
          do isp=1,n_sp
            do i = 0, nxg+1
              do j = 0, nyg+1
                Jpz(i,j,isp) = Jpz(i,j,isp)
     .            + 0.5*eps_pic*sign(1d0,spcs(isp)%q)
     .                     *((2d0*dble(nh1)/Lx)**2
     $                     + (    dble(nh2)/Ly)**2)*pi*pi
     $             *sin(ky*(pyy(j)-ymin))
     .             *cos(kx*(pxx(i)-xmin))
              enddo
            end do
          end do

        else

          !Perturb particles
          do isp = 1,n_sp

            do ip=1,size(spcs(isp)%pcles)
             
              call xform_pcle_idx(spcs(isp)%pcles(ip)%ijk_n,nxg
     .                            ,i_np,j_np)
 
              !VECTOR SIMD
              do k=1,_Npg
                kxp = kx*((i_np(k)-1)*hx+spcs(isp)%pcles(ip)%x_n(k,1))
                kyp = ky*((j_np(k)-1)*hy+spcs(isp)%pcles(ip)%x_n(k,2))

                !perturb ve
                spcs(isp)%pcles(ip)%v_n(k,3)
     .                 = spcs(isp)%pcles(ip)%v_n(k,3) 
     .                  +eps_pic*v0_z(isp)*sin(kyp)*cos(kxp)
                spcs(isp)%pcles(ip)%v_np(k,3)
     .                 = spcs(isp)%pcles(ip)%v_n(k,3) 
              enddo
            end do
          end do

        endif
      endif
