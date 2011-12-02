      module grid_packing

        use grid_petsc

        type :: grid_pack
#if defined(samrai)
          sequence
#endif
          logical    :: pack = .false.
          real(8)    :: xp
          real(8)    :: dx0
        end type grid_pack

        type(grid_pack),dimension(3),save :: g_pack

      contains

c     findGlobalGrid
c     #################################################################
      subroutine findGlobalGrid(lmin,lmax,bcs1,bcs2,xg,dpack,crse_ratio)

c     -----------------------------------------------------------------
c     Finds global grid xg(0:ng+1) in domain (lmin,lmax) with boundary
c     conditions bcs1, bcs2, according to specifications in dpack.
c     The input variable crse_ratio is used to adjust the packing provided
c     in dpack for the fine mesh to coarser grid levels.
c     -----------------------------------------------------------------

        implicit none

c     Call variables

        integer :: crse_ratio,bcs1,bcs2

        real(8) :: xg(0:),lmin,lmax

        type(grid_pack) :: dpack
cc        type(dim_pack) :: dpack

c     Local variables

        integer :: n1,n2,ng,i
        real(8) :: dh,length,lstart,lend,xp,dx0,dxl,dxr
        real(8) :: xi,xil,xir!,alphal=1d0,alphar=1d0
        logical :: pack_l,pack_r,fv_l,fv_r

c     Begin program

        ng = size(xg)-2

        dx0 = dpack%dx0*crse_ratio
        dxl = dpack%dx0*crse_ratio
        dxr = dpack%dx0*crse_ratio

          !Exception for polar coordinates
cc          if ((ig > 1) .and. (.not.dpack%pack) .and. (bcs1 == SP)) then
cc            dpack%pack = .true.
cc            dx0 = (lmax-lmin)/nglobal(1)
cc            dxl = dx0
cc            dxr = dx0
cc            dpack%xp = lmin + dx0
cc          endif

        length = lmax-lmin   !Global length, based on global limits

        fv_l = .not.(abs(bcs1) == SYM)  !FV discr. at opposite boundary except for SYM BC
        fv_r = .not.(abs(bcs2) == SYM)  !FV discr. at opposite boundary except for SYM BC

        if (dpack%pack) then

          dh = length/ng
          lstart = lmin
          lend   = lmax

          xp = dpack%xp

          pack_l = .false.
          pack_r = .false.

          if (xp < lend .and. xp > lstart) then

            n1 = max(int((xp - lstart)/length*ng),1)

            n2 = ng - n1

cc            if (crse_ratio > 1 .and. xil /= 1d0) then
cc              alphal = alphal*(xil**crse_ratio-1d0)/(xil-1d0)
cc            else
cc              alphal = alphal
cc            endif
cc
cc            dxl = alphal*dxl  !Modify initial interval for coarser grids

            if (n1 > 1) then
              call packGrid(n1,lstart,xp-0.5*dxl,-dxl,xg(1:n1),xil,fv_l)
            else
              xg(1) = xp-0.5*dxl
              xil = 1d0
            endif

cc            if (crse_ratio > 1 .and. xir /= 1d0) then
cc              alphar = alphar*(xir**crse_ratio-1d0)/(xir-1d0)
cc            else
cc              alphar = alphar
cc            endif
cc
cc            dxr = alphar*dxr  !Modify initial interval for coarser grids

            if (n2 > 1) then
              call packGrid(n2,xp+0.5*dxr,lend,dxr,xg(n1+1:ng),xir,fv_r)
            else
              xg(ng) = xp+0.5*dxr
              xir = 1d0
            endif

cc            write (*,*) 'alphar,alphal=',alphar,alphal

            if (abs(bcs1) == SP) then
cc              dxl = alphar*dxl
              dxl = crse_ratio*dxl
              dpack%xp = lstart + dxl
            endif

          else

cc            if (crse_ratio > 1 .and. xi /= 1d0) then
cc              alphal = alphal*(xi**crse_ratio-1d0)/(xi-1d0)
cc            else
cc              alphal = alphal
cc            endif
cc
cc            dx0 = alphal*dx0  !Modify initial interval for coarser grids

            if (xp == lstart) then  !Pack to the right
              pack_r = .true.
              if (fv_l) then  !Ensure left boundary is finite-volume type
                call packGrid(ng,lstart+0.5*dx0,lend,dx0,xg(1:ng)
     .                       ,xi,fv_r)
              else
                call packGrid(ng,lstart        ,lend,dx0,xg(1:ng)
     .                       ,xi,fv_r)
              endif
            else                    !Pack to the left
              pack_l = .true.
              if (fv_r) then  !Ensure right boundary is finite-volume type
                call packGrid(ng,lstart,lend-0.5*dx0,-dx0,xg(1:ng)
     .                       ,xi,fv_l)
              else
                call packGrid(ng,lstart,lend        ,-dx0,xg(1:ng)
     .                       ,xi,fv_l)
              endif
            endif

            xil = xi
            xir = xi

          endif

          !X0 ghost cell
          if (bcs1 == PER) then
            xg(0) = xg(1) - (xg(ng)-xg(ng-1))
          elseif (abs(bcs1) == SYM) then
            xg(0) = xg(1) - (xg(2 )-xg(1   ))
          elseif (abs(bcs1) == SP) then
            xg(0) =-xg(1)
          else
            if (pack_r) then
              xg(0) = xg(1) - (xg(2 )-xg(1   ))
            else
              xg(0) = xg(1) - (xg(2 )-xg(1   ))*xil
            endif
          endif

          !Xn ghost cell
          if (bcs2 == PER) then
            xg(ng+1) = xg(ng) + (xg(2) -xg(1)   )
          elseif (abs(bcs2) == SYM) then
            xg(ng+1) = xg(ng) + (xg(ng)-xg(ng-1))
          else
            if (pack_l) then
              xg(ng+1) = xg(ng) + (xg(ng)-xg(ng-1))
            else
              xg(ng+1) = xg(ng) + (xg(ng)-xg(ng-1))*xir
            endif
          endif

        else  !Uniform grid

          dh = length/ng

cc          if (bcs1 == PER .and. bcs2 == PER) then
cc            dh = length/ng
cc            lstart = lmin
cc            lend   = lmax
cc          elseif (bcs1 == SYM .and. bcs2 == SYM) then
          if (bcs1 == SYM .and. bcs2 == SYM) then
            dh = length/(ng-1)
            lstart = lmin
            lend   = lmax
          elseif (bcs1 == SYM .and. bcs2 /= SYM) then
            dh = length/(ng - 0.5)
            lstart = lmin
            lend   = lmax-dh/2.
          elseif (bcs1 /= SYM .and. bcs2 == SYM) then
            dh = length/(ng - 0.5)
            lstart = lmin+dh/2.
            lend   = lmax
          else
            dh = length/ng
            if (bcs1 == PER .and. bcs2 == PER.and.ng == 1) then
              lstart = lmin
              lend   = lmax
            else
              lstart = lmin+dh/2.
              lend   = lmax-dh/2.
            endif
          endif

          xg(1) = lstart
          do i = 2,ng+1
            xg(i) = xg(i-1) + dh
          enddo

          xg(0) = xg(1) - dh

        endif

c     End program

      end subroutine findGlobalGrid

c     packGrid
c     #################################################################
      subroutine packGrid(nn,lmin,lmax,dx0,xx,xi,fv)

c     -----------------------------------------------------------------
c     Finds packed grid, packing to the right if dx0<0 and to the 
c     left if dx0>0. If fv=.true., we pack so that the grid in the
c     unpacked boundary meets the boundary at a face. Otherwise, we
c     pack so that a grid point meets that boundary.
c     -----------------------------------------------------------------

        implicit none

c     Call variables

        integer    :: nn
        real(8)    :: lmin,lmax,dx0,xx(nn),xi
        logical    :: fv

c     Local variables

        integer    :: i
        real(8)    :: rhs

c     Begin program

        rhs = (lmax-lmin)/abs(dx0)

        call find_xi(nn-1,rhs,xi,fv)

        if (dx0 > 0d0) then !Pack to the right

          xx(1) = (1d0+sign(1d-12,lmin))*lmin  !To avoid round-off issues at edges
          do i = 2,nn
cold            xx(i) = xx(i-1) + dx0*xi**(i-2)
            xx(i) = xx(i-1) + dx0*xi**(i-1)
          enddo
            
        else                !Pack to the left

          xx(nn) = (1d0-sign(1d-12,lmax))*lmax  !To avoid round-off issues at edges
          do i = nn,2,-1
cold            xx(i-1) = xx(i) - abs(dx0)*xi**(nn-i)
            xx(i-1) = xx(i) - abs(dx0)*xi**(nn-i+1)
          enddo

        endif

c     End program

      end subroutine packGrid

c     find_xi
c     #################################################################
      subroutine find_xi(nn,rhs,xi,fv)

c     -----------------------------------------------------------------
c     Finds geometric factor xi using Newton's method according to the
c     formula:
c              xi*(xi^nn-1)/(xi-1) = rhs
c     This formula is derived by summing up geometric series starting
c     at xi*dx0.
c
c     The Newton iteration is initialized by using xi=1+eps, and
c     expanding in eps to second order in the formula above to find:
c              eps_0 ~ 2*(rhs-nn)/nn^2
c     This result is used to initialize a single Picard iteration,
c              xi_0 ~ exp(log(1+eps_0*rhs)/nn)
c     which is used as the first guess for Newton.
c     -----------------------------------------------------------------

        implicit none

c     Call variables

        integer    :: nn
        real(8)    :: rhs,xi
        logical    :: fv

c     Local variables

        integer   ,parameter :: nit=30
        integer    :: i
        real(8)    :: ijac,dxi(nit),res(nit)

c     Begin program

c     Initial guess

        xi = 2*(rhs-nn)/nn**2
        xi = exp(log(1+xi*rhs)/nn)

        if (xi == 1d0) return

c     Newton iteration

        do i=1,nit
          if (fv) then  !Do packing in FV fashion at unpacked boundary
cold            res(i) = rhs-0.5*(xi**nn+xi**(nn+1)-2)/(xi-1)
cold            ijac = 2*(xi-1)**2/(nn*(xi**(nn+1)-xi**(nn-1))-2*xi**nn+2)
            res(i) = rhs-0.5*(xi**(nn+1)+xi**(nn+2)-2*xi)/(xi-1)
            ijac=2*(xi-1)**2/((nn+1)*(xi**(nn+2)-xi**nn)-2*xi**(nn+1)+2)
          else          !Position both grid points at both boundaries
cold            res(i) = rhs-(xi**nn-1)/(xi-1)
cold            ijac = (xi-1)**2/((nn-1)*xi**nn-nn*xi**(nn-1)+1)
            res(i) = rhs-xi*(xi**nn-1)/(xi-1)
            ijac = (xi-1)**2/(nn*xi**(nn+1)-(nn+1)*xi**nn+1)
          endif

          dxi(i) = res(i)*ijac
          xi = xi + dxi(i)

          if (abs(res(i)) < 1d-12 .or. (abs(dxi(i)) < 1d-15) ) exit
        enddo

c     Convergence check

        if (i >= nit) then
          write (*,*)
          write (*,*) 'Newton it. in grid-packing routine find_xi'
     .               ,' did not converge'
          write (*,*)
          write (*,*) 'Convergence history:'
          do i=1,nit
            write (*,*) 'Newton it.',i,' residual',abs(res(i))
     .                 ,'update',abs(dxi(i))
          enddo
          write (*,*)
          write (*,*) 'Aborting...'
          stop
        endif

c     End program

      end subroutine find_xi

      end module grid_packing
