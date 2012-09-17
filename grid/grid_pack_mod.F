      module grid_packing

        use math
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

        logical :: map_grid_pack=.false.

        real(8) :: gp_config(6)

      contains

c     scale_coord
c     #################################################################
      function scale_coord(xmin,xmax,x1,dim) result(xlog)

c     -----------------------------------------------------------------
c     Finds global grid xg(0:ng+1) in domain (lmin,lmax) with boundary
c     conditions bcs1, bcs2, according to specifications in dpack.
c     The input variable crse_ratio is used to adjust the packing provided
c     in dpack for the fine mesh to coarser grid levels.
c     -----------------------------------------------------------------

        implicit none

c     Call variables

        integer :: dim
        real(8) :: xlog,x1,xmin,xmax

c     Local variables

        integer :: ibc,exp=3

        real(8) :: gp(2),bbx

cc        real(8),save :: aax=0d0,aay=0d0,aaz=0d0
cc     .                 ,bbx=0d0,bby=0d0,bbz=0d0
cc     .                 ,lbx=0d0,lby=0d0,lbz=0d0

c     Begin program

        ibc = 1+2*(dim-1)

        !Recover grid packing configuration: slope gp(2) at given position gp(1)
        gp = gp_config(ibc:ibc+1)

        if (gp(1) == 0d0 .and. gp(2) == 0d0) then
          xlog = x1
          return
        endif

cc        write (*,*) 'Dim=',dim,'gp=',gp

cc          !Find monotonic atanh profile with given slope (gp(2,4,6)
cc          !  for X, Y, Z) at given position (gp(1,3,5) for X, Y, Z).
cc          if (aax == 0d0) aax = find_a(xmin,xmax,gp(2),gp(1))
cc          if (aay == 0d0) aay = find_a(ymin,ymax,gp(4),gp(3))
cc          if (aaz == 0d0) aaz = find_a(zmin,zmax,gp(6),gp(5))
cc
cc          xx = scl_coord(x1,aax,xmin,xmax,gp(2),gp(1))
cc          yy = scl_coord(x2,aay,ymin,ymax,gp(4),gp(3))
cc          zz = scl_coord(x3,aaz,zmin,zmax,gp(6),gp(5))

c       Packing

        if (bcond(ibc) == PER) then !cosine jacobian
          bbx = (1.-gp(2))  !Slope is in gp(2) (see grid_create_mod.F)
          if (bbx <=-1d0) then
            write (*,*) 'Slope too small along ',dim,' coordinate'
            write (*,*) 'Aborting...'
            stop
          endif
          xlog = x1+(xmax-xmin)/2/pi*bbx*sin(2*pi*(x1-xmin)/(xmax-xmin))
        else  !monotonic cubic polynomial
          bbx= find_eta(xmin,xmax,gp(1),gp(2))
          xlog = scl_coord(x1,bbx,xmin,xmax,gp(1),gp(2))
        endif

      contains

ccc     scl_coord
ccc     #################################################################
cc      function scl_coord(xi,aa,xmin,xmax,x0,lambda) result(x)
cc
ccc     -----------------------------------------------------------------
ccc     Packs coordinate according to x0 (location) and lambda (slope).
ccc     Note:
ccc       * lambda <= 1d0; dx0 ~ dxi*dlambda
ccc       * lambda = 0 --> No packing
ccc     -----------------------------------------------------------------
cc
cc        implicit none
cc
ccc     Global variables
cc
cc        real(8) :: aa,xmin,xmax,x0,lambda,xi,x
cc
ccc     Local variables
cc
cc        real(8) :: bb,xi0
cc
ccc     Begin program
cc
cc        if (lambda == 0d0) then
cc          x = xi
cc        else
cc          bb = aa/lambda
cc          xi0 = xmax - bb*tanh((xmax-x0)/aa)
cc
cc          if (xi > xmax) then
cc            x = 2*xmax -(x0 + aa*atanh((2*xmax-xi-xi0)/bb))
cc          elseif (xi < xmin) then
cc            x = 2*xmin -(x0 + aa*atanh((2*xmin-xi-xi0)/bb))
cc          else
cc            x = x0 + aa*atanh((xi-xi0)/bb)
cc          endif
cccc          write (*,*) xi,x,aa,bb,(xi-xi0)/bb
cc        endif
cc
ccc     End program
cc
cc      end function scl_coord
cc
ccc     find_a
ccc     #################################################################
cc      function find_a(xmin,xmax,x0,cc) result(aa)
cc
cc        implicit none
cc
ccc     Global variables
cc
cc        real(8) :: xmin,xmax,x0,cc,aa
cc
ccc     Local variables
cc
cc        integer,parameter :: nit=20
cc
cc        integer :: i
cc        real(8) :: ijac,dxi(nit),res(nit)
cc        real(8) :: ll,l0m,l0p,dd,ee,ff,bb
cc
ccc     Begin program
cc
cc        if (cc == 0d0) return
cc
cc        ll  = xmax-xmin
cc        l0m = x0  -xmin
cc        l0p = xmax - x0
cc
cc        dd = cc*ll
cc
cc        if (cc > 1d0) then
cccc          call pstop('find_a in grid_mod.F',
cccc     .               'Slope has to <= 1')
cc          write (*,*) 'find_a in grid_anal_map_mod.F: Slope has to <= 1'
cc          write (*,*) 'Aborting...'
cc          stop
cc        endif
cc
ccc     Initial guess
cc
cc        aa = 0.5*ll
cc
ccc     Newton iteration
cc
cc        do i=1,nit
cc          ff = dd/aa
cc
cc          res(i) = dd - aa*(tanh(l0m/aa) + tanh(l0p/aa))
cc          ijac = 1./(-tanh(l0m/aa)-tanh(l0p/aa)
cc     .              +l0m/aa/cosh(l0m/aa)**2
cc     .              +l0p/aa/cosh(l0p/aa)**2)
cc
cc          dxi(i) = -res(i)*ijac
cc          aa = abs(aa + dxi(i))  !Want to find positive root
cc
cccc          write (*,*) abs(res(i)),dxi(i),ijac,aa
cc
cc          if (abs(res(i)) < 1d-10 .or. (abs(dxi(i)) < 1d-15) ) exit
cc        enddo
cc
ccc     Convergence check
cc
cc        if (i >= nit) then
cc          write (*,*)
cc          write (*,*) 'Newton it. in grid-packing routine find_a'
cc     .               ,' did not converge'
cc          write (*,*)
cc          write (*,*) 'Convergence history:'
cc          do i=1,nit
cc            write (*,*) 'Newton it.',i,' residual',abs(res(i))
cc     .                 ,'update',abs(dxi(i))
cc          enddo
cc          write (*,*)
cc          write (*,*) 'Aborting...'
cc          stop
cc        endif
cc
cc      end function find_a

c     scl_coord
c     #################################################################
      function scl_coord(xi,eta,xmin,xmax,x0,lambda) result(x)

c     -----------------------------------------------------------------
c     Packs coordinate according to x0 (location) and lambda (slope).
c     Note:
c       * lambda <> 1d0 --> dx0 ~ dxi*dlambda
c       * lambda == 0d0 --> No packing
c     -----------------------------------------------------------------

        implicit none

c     Global variables

        real(8),INTENT(IN) :: eta,xmin,xmax,x0,lambda,xi
        real(8) :: x

c     Local variables

        real(8) :: xi0,xibar,LL,bb,dx

        integer,parameter :: nit=100

        integer :: i
        real(8) :: ijac,dxi(nit),res(nit)

c     Begin program

        if (eta == 0d0) then

          x = xi

        else

          LL = xmax-xmin
          bb = eta/LL**2

          if (lambda < 1d0) then   !Pack at x0
            xibar = 0.5*(xmax+xmin)

            xi0 = xibar + sign(1d0,x0-xibar)
     .                   *LL*sqrt(4d0*(1d0-lambda)-eta)/sqrt(12.*eta)

            !To ensure periodicity of map
cc            if (xi > xmax) then
cc              dx = xi - xmax
cc              x =  x0 + lambda*(xmax-dx-xi0) + bb*(xmax-dx-xi0)**3
cc              x = 2*xmax - x
cc            elseif (xi < xmin) then
cc              dx = xmin - xi
cc              x =  x0 + lambda*(xmin+dx-xi0) + bb*(xmin+dx-xi0)**3
cc              x = 2*xmin - x
cc            else
              x = x0 + lambda*(xi-xi0) + bb*(xi-xi0)**3
cc            endif

cc            write (*,*) xi0,eta,LL,xmin,xmax,lambda,xi,x
          else                    !Pack at edges

            xi0 = xmin - (xmin-x0)/lambda - bb*(xmin-x0)**exp

            !To ensure periodicity of map
cc            if (xi > xmax) then
cc              dx = xi - xmax
cc            elseif (xi < xmin) then
cc              dx = xi - xmin
cc            else
cc              dx = 0d0
cc            endif

            !Initial guess
            x = xi! - 2*dx

            !Newton iteration to find x(xi)
            do i=1,nit
              res(i) = xi-2*dx-(xi0 + (x-x0)/lambda + bb*(x-x0)**exp)

              ijac = -1d0/lambda -exp*bb*(x-x0)**(exp-1)

              dxi(i) = -res(i)/ijac
              x = x + dxi(i)

cc              write (*,*) abs(res(i)),dxi(i),x

              if (abs(res(i)) < 1d-8 .or. (abs(dxi(i)) < 1d-14) ) exit
            enddo

            !Convergence check
            if (i > nit) then
              write (*,*)
              write (*,*)'Newton it. in grid-packing routine',
     .                   ' scl_coord did not converge'
              write (*,*)
              write (*,*) 'Convergence history:'
              do i=1,nit
                write (*,*) 'Newton it.',i,' residual',abs(res(i))
     .                     ,'update',abs(dxi(i))
              enddo
              write (*,*)
              write (*,*) 'Aborting...'
              stop
            endif

cc            !To ensure periodicity of map
cc            if (xi > xmax) then
cc              x = 2*xmax - x
cc            elseif (xi < xmin) then
cc              x = 2*xmin - x
cc            endif

          endif

        endif

c     End program

      end function scl_coord

c     find_eta
c     #################################################################
      function find_eta(xmin,xmax,x0,lambda) result(eta)

        implicit none

c     Global variables

        real(8) :: xmin,xmax,x0,lambda,eta

c     Local variables

        integer,parameter :: nit=20

        integer :: i
        real(8) :: ijac,dxi(nit),res(nit)
        real(8) :: ll,dd,xbar,dx

c     Begin program

        if (lambda == 0d0 .or. lambda == 1d0) then
          eta = 0d0
          return
        endif

        LL  = xmax-xmin

        if (lambda > 1d0) then  !Pack at edges

cc          eta = LL**3*(1d0-1d0/lambda)/((xmax-x0)**3-(xmin-x0)**3)
          eta = LL**exp*(1d0-1d0/lambda)/((xmax-x0)**exp-(xmin-x0)**exp)

        else                    !Pack at x0

          xbar= 0.5*(xmax+xmin)
          dx  = xbar - x0
          dd  = (dx/LL)**2

c         Initial guess

          eta = 4d0*(1d0-lambda)

c         Newton iteration

          if (dx == 0d0) return

          do i=1,nit
            res(i)=(1+2*(eta+lambda))**2*(1-lambda-0.25*eta) -27.*dd*eta
            ijac = 4*(1+2*(eta+lambda)) *(1-lambda-0.25*eta)
     .         -0.25*(1+2*(eta+lambda))**2                   -27.*dd

            dxi(i) = -res(i)/ijac
            eta = abs(eta + dxi(i)) !Want to find positive root

cc            write (*,*) abs(res(i)),dxi(i),ijac,eta

            if (abs(res(i)) < 1d-10 .or. (abs(dxi(i)) < 1d-15) ) exit
          enddo

c         Convergence check

          if (i >= nit) then
            write (*,*)
            write (*,*) 'Newton it. in grid-packing routine find_a'
     .                 ,' did not converge'
            write (*,*)
            write (*,*) 'Convergence history:'
            do i=1,nit
              write (*,*) 'Newton it.',i,' residual',abs(res(i))
     .                   ,'update',abs(dxi(i))
            enddo
            write (*,*)
            write (*,*) 'Aborting...'
            stop
          endif

        endif

      end function find_eta

      end function scale_coord

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

cc            if (abs(bcs1) == SP) then
cccc              dxl = alphar*dxl
cc              dxl = crse_ratio*dxl
cc              dpack%xp = lstart + dxl
cc            endif

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

cc          if (bcs1 == PER .and. bcs2 == PER) then  !Need to change setGraphicsRange routine as well
cc            dh = length/ng
cc            lstart = lmin
cc            lend   = lmax
cc          elseif (bcs1 == SYM .and. bcs2 == SYM) then
          if (bcs1 == SYM .and. bcs2 == SYM) then
            dh = length/(ng-1)
            lstart = lmin
cc            lend   = lmax
          elseif (bcs1 == SYM .and. bcs2 /= SYM) then
            dh = length/(ng - 0.5)
            lstart = lmin
cc            lend   = lmax-dh/2.
          elseif (bcs1 /= SYM .and. bcs2 == SYM) then
            dh = length/(ng - 0.5)
            lstart = lmin+dh/2.
cc            lend   = lmax
          else
            dh = length/ng
            if (bcs1 == PER .and. bcs2 == PER.and.ng == 1) then
              lstart = lmin
cc              lend   = lmax
            else
              lstart = lmin+dh/2.
cc              lend   = lmax-dh/2.
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
