c module grid_packing
c #####################################################################
      module grid_packing

        use math
        use grid_def_st

        type(grid_pack),dimension(3),save :: g_pack

      contains

c     set_pack_var
c     #################################################################
      subroutine set_pack_slope(gdef)

      implicit none

      type(grid_mg_def),pointer :: gdef

      integer :: nxx,nyy,nzz
      real(8) :: xmin,xmax,ymin,ymax,zmin,zmax

      nxx = gdef%nglx
      nyy = gdef%ngly
      nzz = gdef%nglz

      xmin = gdef%gxmin
      xmax = gdef%gxmax
                         
      ymin = gdef%gymin
      ymax = gdef%gymax
                         
      zmin = gdef%gzmin
      zmax = gdef%gzmax

      if (gdef%g_pack(1)%pack
     .   .and.(.not.gdef%g_pack(1)%processed)) then
        gdef%g_pack(1)%slope = gdef%g_pack(1)%dx0/(xmax-xmin)*nxx
        gdef%g_pack(1)%processed = .true.
        g_pack(1) = gdef%g_pack(1)
      endif

      if (gdef%g_pack(2)%pack
     .   .and.(.not.gdef%g_pack(2)%processed)) then
        gdef%g_pack(2)%slope = gdef%g_pack(2)%dx0/(ymax-ymin)*nyy
        gdef%g_pack(2)%processed = .true.
        g_pack(2) = gdef%g_pack(2)
      endif

      if (gdef%g_pack(3)%pack
     .   .and.(.not.gdef%g_pack(3)%processed)) then
        gdef%g_pack(3)%slope = gdef%g_pack(3)%dx0/(zmax-zmin)*nzz
        gdef%g_pack(3)%processed = .true.
        g_pack(3) = gdef%g_pack(3)
      endif
 
      end subroutine set_pack_slope

c     scale_coord
c     #################################################################
      function scale_coord(gdef,dim,x1) result(xlog)

c     -----------------------------------------------------------------
c     Scales coordinate along dimension dim and evaluates coordinate at
c     logical coordinate x1. Packing configuration (location and slope
c     of packing) is taken from g_def.
c     -----------------------------------------------------------------

        implicit none

c     Call variables

        type(grid_mg_def),pointer :: gdef

        integer :: dim
        real(8) :: xlog,x1

c     Local variables

        integer :: ibc,ibc1,exp=3,nn

        real(8) :: gp(2),lmin,lmax,TT
        real(8), save :: bbx=0d0,bby=0d0,bbz=0d0

        logical :: pack_both_bdrys

c     Begin program

        if (.not.gdef%g_pack(dim)%pack) then
          xlog = x1
          return
        endif

        ibc = 1+2*(dim-1)
        ibc1 = 2+2*(dim-1)

c     Grid packing configuration: slope gp(2) at given position gp(1)

        select case(dim)
        case(1)
          nn   = gdef%nglx
          lmin = gdef%gxmin
          lmax = gdef%gxmax
        case(2)
          nn   = gdef%ngly
          lmin = gdef%gymin
          lmax = gdef%gymax
        case(3)
          nn   = gdef%nglz
          lmin = gdef%gzmin
          lmax = gdef%gzmax
        end select

        gp(1) = gdef%g_pack(dim)%xp
        gp(2) = gdef%g_pack(dim)%slope

        pack_both_bdrys = gdef%g_pack(dim)%pack_both_bdrys
        TT = 2*pi*gdef%g_pack(dim)%modes

        if (.not.gdef%g_pack(dim)%processed) then
          write (*,*) 'Slope not computed in ',dim,' coordinate'
          write (*,*) 'Aborting...'
          stop
        endif

c     Packing

        if (bcond(ibc) == PER) then !cosine jacobian
          bbx = abs(1.-gp(2))
          if (bbx >= 1d0) then
            write (*,*) 'Slope too small along ',dim,' coordinate'
            write (*,*) 'Aborting...'
            stop
          endif

          xlog =x1-(lmax-lmin)/TT*bbx*sin(TT*(x1-gp(1))/(lmax-lmin))
        elseif (pack_both_bdrys) then !sym cosine jacobian
          bbx = abs(1.-gp(2))
          if (bbx >= 1d0) then
            write (*,*) 'Slope too small along ',dim,' coordinate'
            write (*,*) 'Aborting...'
            stop
          endif

          xlog = x1-(lmax-lmin)/TT*bbx*sin(TT*(x1-lmin)/(lmax-lmin))
        else  !monotonic cubic polynomial
          select case(dim)
          case(1)
            if (bbx == 0d0) bbx= find_eta(lmin,lmax,gp(1),gp(2))
            xlog = scl_coord(x1,bbx,lmin,lmax,gp(1),gp(2))
          case(2)
            if (bby == 0d0) bby= find_eta(lmin,lmax,gp(1),gp(2))
            xlog = scl_coord(x1,bby,lmin,lmax,gp(1),gp(2))
          case(3)
            if (bbz == 0d0) bbz= find_eta(lmin,lmax,gp(1),gp(2))
            xlog = scl_coord(x1,bbz,lmin,lmax,gp(1),gp(2))
          end select
        endif

      contains

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

          !Pack at x0
          if (lambda < 1d0) then
            xibar = 0.5*(xmax+xmin)

            xi0 = xibar + sign(1d0,x0-xibar)
     .                   *LL*sqrt(4d0*(1d0-lambda)-eta)/sqrt(12.*eta)

            !Linear extrapolation to ensure equal mesh spacing at ghost points
            if (xi > xmax) then
              dx = xi - xmax
              x =  x0 + lambda*(xmax-dx-xi0) + bb*(xmax-dx-xi0)**3
              x = 2*xmax - x
            elseif (xi < xmin) then
              dx = xmin - xi
              x =  x0 + lambda*(xmin+dx-xi0) + bb*(xmin+dx-xi0)**3
              x = 2*xmin - x
            else
              x = x0 + lambda*(xi-xi0) + bb*(xi-xi0)**3
            endif

          !Pack at edges
          else

            xi0 = xmin - (xmin-x0)/lambda - bb*(xmin-x0)**exp

            !Linear extrapolation to ensure equal mesh spacing at ghost points
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

            if (n1 > 1) then
              call packGrid(n1,lstart,xp-0.5*dxl,-dxl,xg(1:n1),xil,fv_l)
            else
              xg(1) = xp-0.5*dxl
              xil = 1d0
            endif

            if (n2 > 1) then
              call packGrid(n2,xp+0.5*dxr,lend,dxr,xg(n1+1:ng),xir,fv_r)
            else
              xg(ng) = xp+0.5*dxr
              xir = 1d0
            endif

          else

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

          if (bcs1 == PER.and.bcs2 == PER.and.ng == 1) then !Single layer dimension
cc          if (bcs1 == PER.and.bcs2 == PER) then !Sets angle origin at zero
            dh = length/ng
            lstart = lmin
          elseif (bcs1 == SYM .and. bcs2 == SYM) then
            dh = length/(ng-1)
            lstart = lmin
          elseif (bcs1 == SYM .and. bcs2 /= SYM) then
            dh = length/(ng - 0.5d0)
            lstart = lmin
          elseif (bcs1 /= SYM .and. bcs2 == SYM) then
            dh = length/(ng - 0.5d0)
            lstart = lmin+0.5d0*dh
          else
            dh = length/ng
            lstart = lmin+0.5d0*dh
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
