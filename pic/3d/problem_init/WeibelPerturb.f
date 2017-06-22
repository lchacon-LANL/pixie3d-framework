!     only with quiet_start

      isp = 1 !electrons only

      pcles => spcs(isp)%pcles

      kx = 2d0*acos(-1d0)*dble(nh1)/Lx

      do ip=1,size(pcles)
        xn = pcles(ip)%x_n
        i_np = pcles(ip)%i_n 
        j_np = pcles(ip)%j_n 
        xp = xn + pxx(i_np)  - 0.5d0*hx
        kxp = kx*xp      


c$$$         if(pcles(ip)%v_n(2) ==0) then
c$$$            print *,ip,pcles(ip)%v_n(2),eps_pic*cos(kxp)
c$$$         end if
c$$$         pcles(ip)%v_n(2) = pcles(ip)%v_n(2) !jy only
c$$$     .        + eps_pic*cos(kxp)
c$$$         pcles(ip)%v_np(2)= pcles(ip)%v_n(2)            

!test y
c$$$         kxp = 2d0*acos(-1d0)*dble(nh1)/Ly*( pcles(ip)%x_n(2) +
c$$$     $        pyy(pcles(ip)%j_n)  - 0.5d0*hy )

c$$$            pcles(ip)%v_n(3) = pcles(ip)%v_n(3) !jz only
c$$$     .           + eps_pic *cos(kxp)
c$$$            pcles(ip)%v_np(3)= pcles(ip)%v_n(3)            

c         if(i_np==1.and.j_np==1) then
        if(i_np==nxg/2+1.and.j_np==nyg/2+1) then
          pcles(ip)%v_nz = pcles(ip)%v_nz !jz only
     .                      + eps_pic        !*cos(kxp)
          pcles(ip)%v_npz= pcles(ip)%v_nz           
        end if

      end do

      !Initialize density profiles for all species
c$$$      do isp=1,size(spcs)
c$$$        pcles => spcs(isp)%pcles
c$$$        ns(:,:,isp) = (size(pcles)-nip)/(nxg*nyg*nzg)
c$$$      enddo

cc      !Build 9-point stencil matrix
cc      MM2 = 0d0
cc      rr  = 0d0
cc      do isp=1,size(spcs)
cccc        aa = 0d0
cccc        bb = 0d0
cccc        cc = 0d0
cccc        MM = 0d0 
cc        pcles => spcs(isp)%pcles
cc
cc         do ip=1,size(pcles)
cc            if(ip>nip) then 
cc               xp = pcles(ip)%x_n(1) / hx
cc               yp = pcles(ip)%x_n(2) / hy
cc               Mx2(:,1) = shape2_offset(xp)
cc               My2(1,:) = shape2_offset(yp)
cc               Mxy2(:,:)= matmul(Mx2,My2)
cc
cc               ii = pcles(ip)%i_n  + nxg*(pcles(ip)%j_n-1)
cc
cc               MM2(ii,7:9,isp) = MM2(ii,7:9,isp)+Mxy2(1:3,3)
cc               MM2(ii,4:6,isp) = MM2(ii,4:6,isp)+Mxy2(1:3,2)
cc               MM2(ii,1:3,isp) = MM2(ii,1:3,isp)+Mxy2(1:3,1)
cc
ccc$$$               ci_ngp=pcles(ip)%i_n 
ccc$$$               ci_nm1 = ci_ngp-1
ccc$$$               ci_np1 = ci_ngp+1
ccc$$$              if(ci_ngp==1) then
ccc$$$                 cc(nxg   ,1)=cc(nxg   ,1)+Mxy2(1,1)+Mxy2(1,2)+Mxy2(1,3)
ccc$$$                 bb(ci_ngp,1)=bb(ci_ngp,1)+Mxy2(2,1)+Mxy2(2,2)+Mxy2(2,3)
ccc$$$                 aa(ci_np1,1)=aa(ci_np1,1)+Mxy2(3,1)+Mxy2(3,2)+Mxy2(3,3)
ccc$$$              else if(ci_ngp==nxg) then 
ccc$$$                 cc(ci_nm1,1)=cc(ci_nm1,1)+Mxy2(1,1)+Mxy2(1,2)+Mxy2(1,3)
ccc$$$                 bb(ci_ngp,1)=bb(ci_ngp,1)+Mxy2(2,1)+Mxy2(2,2)+Mxy2(2,3)
ccc$$$                 aa(1     ,1)=aa(1     ,1)+Mxy2(3,1)+Mxy2(3,2)+Mxy2(3,3)
ccc$$$              else
ccc$$$                 cc(ci_nm1,1)=cc(ci_nm1,1)+Mxy2(1,1)+Mxy2(1,2)+Mxy2(1,3)
ccc$$$                 bb(ci_ngp,1)=bb(ci_ngp,1)+Mxy2(2,1)+Mxy2(2,2)+Mxy2(2,3)
ccc$$$                 aa(ci_np1,1)=aa(ci_np1,1)+Mxy2(3,1)+Mxy2(3,2)+Mxy2(3,3)
ccc$$$              end if
cc
ccc$$$               ci_ngp=pcles(ip)%i_n 
ccc$$$               ci_nm1 = ci_ngp-1
ccc$$$               ci_np1 = ci_ngp+1
ccc$$$               if(ci_ngp==1) then
ccc$$$                  cc(nxg   ,1) = cc(nxg   ,1) + Mx2(1,1)
ccc$$$                  bb(ci_ngp,1) = bb(ci_ngp,1) + Mx2(2,1)
ccc$$$                  aa(ci_np1,1) = aa(ci_np1,1) + Mx2(3,1)
ccc$$$               else if(ci_ngp==nxg) then 
ccc$$$                  cc(ci_nm1,1) = cc(ci_nm1,1) + Mx2(1,1)
ccc$$$                  bb(ci_ngp,1) = bb(ci_ngp,1) + Mx2(2,1)
ccc$$$                  aa(1     ,1) = aa(1     ,1) + Mx2(3,1)
ccc$$$               else
ccc$$$                  cc(ci_nm1,1) = cc(ci_nm1,1) + Mx2(1,1)
ccc$$$                  bb(ci_ngp,1) = bb(ci_ngp,1) + Mx2(2,1)
ccc$$$                  aa(ci_np1,1) = aa(ci_np1,1) + Mx2(3,1)
ccc$$$               end if
cc
cc            end if
cc         end do
cc
cc         !Build rhs taking Jacobian into account
cc         do cj = 1,nyg
cc           do ci_ngp=1,nxg
cc             ii = ci_ngp + nxg*(cj-1)
cc
cc             rr(ii,isp) = (size(pcles)-nip) /(nxg*nyg*nzg)
cc     .                   *gmetric%grid(1)%jac(ci_ngp,cj,1)
cc           enddo
cc         enddo
cc
cccc         !Collapse to 3-point for 1D
cccc         do cj = 1,nyg
cccc           do ci_ngp=1,nxg
cccc             ci_nm1 = ci_ngp-1
cccc             ci_np1 = ci_ngp+1
cccc            
cccc             ii = ci_ngp + nxg*(cj-1)
cccc             if(ci_ngp==1) then
cccc                cc(nxg   ,1) = cc(nxg   ,1) + MM(ii,1)+MM(ii,4)+MM(ii,7)
cccc                bb(ci_ngp,1) = bb(ci_ngp,1) + MM(ii,2)+MM(ii,5)+MM(ii,8)
cccc                aa(ci_np1,1) = aa(ci_np1,1) + MM(ii,3)+MM(ii,6)+MM(ii,9)
cccc              else if(ci_ngp==nxg) then 
cccc                cc(ci_nm1,1) = cc(ci_nm1,1) + MM(ii,1)+MM(ii,4)+MM(ii,7)
cccc                bb(ci_ngp,1) = bb(ci_ngp,1) + MM(ii,2)+MM(ii,5)+MM(ii,8)
cccc                aa(1     ,1) = aa(1     ,1) + MM(ii,3)+MM(ii,6)+MM(ii,9)
cccc              else
cccc                cc(ci_nm1,1) = cc(ci_nm1,1) + MM(ii,1)+MM(ii,4)+MM(ii,7)
cccc                bb(ci_ngp,1) = bb(ci_ngp,1) + MM(ii,2)+MM(ii,5)+MM(ii,8)
cccc                aa(ci_np1,1) = aa(ci_np1,1) + MM(ii,3)+MM(ii,6)+MM(ii,9)
cccc              end if
cccc            end do
cccc         end do
cccc
cccc         !solve cornered-tridiagonal system
cccc         call BTDS(nxg,neq,nr_bt,aa,bb,cc,r,sol(:,isp))
cc
cc      end do
cc
cccc      call CLEAN_BTDS()
cc
cc      !Solve 9-point stencil problem
cc      do isp=1,size(spcs)
cc        bcs(:,isp) = bcond
cc      enddo
cc
cc      where (bcs == DEF) bcs = NEU
cc      iiout = 0 !No output
cc      call cSolver(size(spcs),nxg*nyg*nzg,rr,sol,bcs
cc     .            ,1,iiout,0,massmat_mtvc,.false.
cc     .            ,tol           = 1d-10
cc     .            ,gm_driver     = .true.
cc     .            ,ks_it         = 100
cccc     .            ,mg_galerkin   = .true.
cccc     .            ,mg_order_res  = 0
cccc     .            ,mg_order_prol = 0
cccc     .            ,mg_coarse_grid_size = mg_coarse_size
cccc     .            ,mg_gm_coarse_solve  = .true.
cc     .            ,mg_vcyc       = 0  !Identity preconditioner
cccc     .            ,mg_smooth     = 'gm'
cccc     .            ,sm_it         = sm_iter
cccccc     .            ,sm_omega      = sm_omega
cccc     .            ,sm_ncolors    = 4
cc     .            ,iters = it_count
cc     .            )  
cc
cc      if (iiout == 0.and.my_rank_pic == 0) then
cc        write (*,*)
cc        write (*,*) 'Mass matrix solve converged in ',it_count,' its'
cc        write (*,*) 
cc      endif
cc
cc      print *,"#", sum(sol(1:nxg*nyg,1)), sum(sol(1:nxg*nyg,2))
cc      print *,"#grid weight"
cc      do ci_ngp = 1,nxg
cc         print *,ci_ngp, sol(ci_ngp,1:2),gmetric%grid(1)%jac(ci_ngp,1,1)
cc      end do
cc      print *
cc      print *
cc      stop
cc
cc      do isp=1, size(spcs) 
cc         pcles => spcs(isp)%pcles
cc         do ip=1,size(pcles)
cc            if(ip.le.nip) then  !reserved inflow particles
cc               pcles(ip)%w = 1d0
cc            else
cc               ci_ngp=pcles(ip)%i_n
cc               pcles(ip)%w = sol(ci_ngp,isp)
cc            end if
cc         end do
cc      end do
