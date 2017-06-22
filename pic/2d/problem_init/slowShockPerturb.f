
cc      !Initialize density profiles
cc      do isp=1,size(spcs)
cc        pcles => spcs(isp)%pcles
cc        ns(:,:,isp) = (size(pcles)-nip) /(nxg*nyg*nzg)
cc      enddo

cc      do isp=1,size(spcs)
cc         aa = 0d0
cc         bb = 0d0
cc         cc = 0d0
cc
cc         pcles => spcs(isp)%pcles
cc         do ip=1,size(pcles)
cc            if(ip>nip) then 
cc               xp = pcles(ip)%x_n(1) / hx
cc               Mx2(:,1) = shape2_offset(xp)
cc               ci_ngp=pcles(ip)%i_n
cc               ci_nm1 = ci_ngp-1
cc               ci_np1 = ci_ngp+1
cc
cc               if(ci_ngp==1) then
cc                  bb(ci_ngp,1) = bb(ci_ngp,1) + Mx2(2,1)+Mx2(1,1)
cc                  aa(ci_np1,1) = aa(ci_np1,1) + Mx2(3,1)
cc               else if(ci_ngp==nxg) then 
cc                  cc(ci_nm1,1) = cc(ci_nm1,1) + Mx2(1,1)
cc                  bb(ci_ngp,1) = bb(ci_ngp,1) + Mx2(2,1)+Mx2(3,1)
cc               else
cc                  cc(ci_nm1,1) = cc(ci_nm1,1) + Mx2(1,1)
cc                  bb(ci_ngp,1) = bb(ci_ngp,1) + Mx2(2,1)
cc                  aa(ci_np1,1) = aa(ci_np1,1) + Mx2(3,1)
cc               end if
cc
ccc$$$               if(ci_ngp==1) then
ccc$$$                  ci_nm1 = nxg
ccc$$$               else if(ci_ngp==nxg) then
ccc$$$                  ci_np1 = 1
ccc$$$               end if
ccc$$$               
ccc$$$               cc(ci_nm1,1) = cc(ci_nm1,1) + Mx2(1,1)
ccc$$$               bb(ci_ngp,1) = bb(ci_ngp,1) + Mx2(2,1)
ccc$$$               aa(ci_np1,1) = aa(ci_np1,1) + Mx2(3,1)
cc            end if
cc         end do
cc!     right hand side
cc         r = (size(pcles)-nip) / (nxg*nyg*nzg)
cc     .       *gmetric%grid(1)%jac(1:nxg,1,1)
cc
cc!     solve cornered-tridiagnoal system
cc         call BTDS(nxg,neq,nr_bt,aa,bb,cc,r,sol(:,isp))
cc
cc      end do
cc
cc
ccc$$$      print *,"#", sum(sol(1:nxg))
ccc$$$      print *,"#grid weight"
ccc$$$      do ci_ngp = 1,nxg
ccc$$$         print *,ci_ngp, sol(ci_ngp)
ccc$$$      end do
ccc$$$      print *
ccc$$$      print *
cc
cc      do isp=1, size(spcs) 
cc         pcles => spcs(isp)%pcles
cc         do ip=1,size(pcles)
cc            if(ip.le.nip) then  !reserved inflow particles
cc               pcles(ip)%w = sol(1,isp)
cc            else
cc               ci_ngp=pcles(ip)%i_n
cc               pcles(ip)%w = sol(ci_ngp,isp)
cc            end if
cc         end do
cc      end do
cc
cc      call CLEAN_BTDS()
cc
