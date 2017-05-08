c module grid_packing
c #####################################################################
      module grid_packing

        use math
        use grid_def_st

        type(grid_pack),dimension(3),private :: g_pack

        real(8), private :: eta(3)=0d0
        integer, private :: pwr=3

      contains

c     init_grid_pack_config
c     #################################################################
      subroutine init_grid_pack_config(gp1,gp2,gp3)

c     -----------------------------------------------------------------
c     Computes slope for grid packing.
c     -----------------------------------------------------------------

      implicit none

      type(grid_pack) :: gp1,gp2,gp3
        
c     Begin program

      g_pack(1) = gp1
      g_pack(2) = gp2
      g_pack(3) = gp3

c     End program

      end subroutine init_grid_pack_config

c     get_grid_pack_config
c     #################################################################
      subroutine get_grid_pack_config(gp)

c     -----------------------------------------------------------------
c     Computes slope for grid packing.
c     -----------------------------------------------------------------

      implicit none

      type(grid_pack),dimension(3) :: gp
        
c     Begin program

      gp = g_pack

c     End program

      end subroutine get_grid_pack_config

c     setup_grid_packing
c     #################################################################
      function setup_grid_packing(map,gdef) result(map_grid_pack)

c     -----------------------------------------------------------------
c     Computes slope for grid packing.
c     -----------------------------------------------------------------

      implicit none

      character(*) :: map

      type(grid_mg_def),pointer :: gdef

      logical :: map_grid_pack

c     Begin program

        if (    gdef%g_pack(1)%pack
     .      .or.gdef%g_pack(2)%pack
     .      .or.gdef%g_pack(3)%pack) then

          map_grid_pack = (map /= 'lpk')  !Do not pack via map if packing logical mesh ("lpk")

          !Save packing configuration
          call find_pack_slope(gdef)
        else
          map_grid_pack = .false.
        endif

c     End program

      end function setup_grid_packing

c     find_pack_slope
c     #################################################################
      subroutine find_pack_slope(gdef)

c     -----------------------------------------------------------------
c     Computes slope for grid packing.
c     -----------------------------------------------------------------

      implicit none

      type(grid_mg_def),pointer :: gdef

      integer :: nn(3),dim,ibc
      real(8) :: lmin(3),lmax(3),idh

      nn(1) = gdef%nglx
      nn(2) = gdef%ngly
      nn(3) = gdef%nglz

      lmin(1) = gdef%gxmin
      lmax(1) = gdef%gxmax
                         
      lmin(2) = gdef%gymin
      lmax(2) = gdef%gymax
                         
      lmin(3) = gdef%gzmin
      lmax(3) = gdef%gzmax

      do dim=1,3
        if (gdef%g_pack(dim)%pack
     .   .and.(.not.gdef%g_pack(dim)%processed)) then
          idh = nn(dim)/(lmax(dim)-lmin(dim))
          gdef%g_pack(dim)%slope = gdef%g_pack(dim)%dx0*idh
          gdef%g_pack(dim)%processed = .true.

          g_pack(dim) = gdef%g_pack(dim) !Save config

          ibc  = 1+2*(dim-1)
          if ((bcond(ibc)/=PER).and.(.not.g_pack(dim)%pack_both_bdrys))
     .       eta(dim)=find_eta(lmin(dim),lmax(dim)
     .                        ,g_pack(dim)%xp,g_pack(dim)%slope)
        endif
      enddo

      contains

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
        real(8) :: ll,dd,xbar,dx,eps

c     Begin program

        eps = 10*findRoundOff()

        if (lambda == 0d0 .or. lambda == 1d0) then
          eta = 0d0
          return
        endif

        LL  = xmax-xmin

        if (lambda > 1d0) then  !Pack at edges

cc          eta = LL**3*(1d0-1d0/lambda)/((xmax-x0)**3-(xmin-x0)**3)
          eta = LL**pwr*(1d0-1d0/lambda)/((xmax-x0)**pwr-(xmin-x0)**pwr)

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

            if (abs(res(i))<1d-10.or.(abs(dxi(i))<eps)) exit
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

      end subroutine find_pack_slope

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

        integer :: ibc,nn

        real(8) :: gp(2),lmin,lmax,TT,bb

        logical :: pack_both_bdrys

c     Begin program

        if (.not.gdef%g_pack(dim)%pack) then
          xlog = x1
          return
        endif

        ibc  = 1+2*(dim-1)

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
          bb = abs(1.-gp(2))
          if (bb >= 1d0) then
            write (*,*) 'Slope too small along ',dim,' coordinate'
            write (*,*) 'Aborting...'
            stop
          endif

          xlog =x1-(lmax-lmin)/TT*bb*sin(TT*(x1-gp(1))/(lmax-lmin))
        elseif (pack_both_bdrys) then !sym cosine jacobian
          bb = abs(1.-gp(2))
          if (bb >= 1d0) then
            write (*,*) 'Slope too small along ',dim,' coordinate'
            write (*,*) 'Aborting...'
            stop
          endif

          xlog = x1-(lmax-lmin)/TT*bb*sin(TT*(x1-lmin)/(lmax-lmin))
        else  !monotonic cubic polynomial
          xlog = scl_coord(x1,eta(dim),lmin,lmax,gp(1),gp(2))
c$$$          select case(dim)
c$$$          case(1)
c$$$            if (bbx == 0d0) bbx= find_eta(lmin,lmax,gp(1),gp(2))
c$$$            xlog = scl_coord(x1,bbx,lmin,lmax,gp(1),gp(2))
c$$$          case(2)
c$$$            if (bby == 0d0) bby= find_eta(lmin,lmax,gp(1),gp(2))
c$$$            xlog = scl_coord(x1,bby,lmin,lmax,gp(1),gp(2))
c$$$          case(3)
c$$$            if (bbz == 0d0) bbz= find_eta(lmin,lmax,gp(1),gp(2))
c$$$            xlog = scl_coord(x1,bbz,lmin,lmax,gp(1),gp(2))
c$$$          end select
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

            xi0 = xmin - (xmin-x0)/lambda - bb*(xmin-x0)**pwr

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
              res(i) = xi-2*dx-(xi0 + (x-x0)/lambda + bb*(x-x0)**pwr)

              ijac = -1d0/lambda -pwr*bb*(x-x0)**(pwr-1)

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
          elseif (bcs1 == SYM .and. bcs2 == SYM.and.ng > 1) then
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
        real(8)    :: ijac,dxi(nit),res(nit),eps

c     Begin program

        eps = 10*findRoundOff()

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

          if (abs(res(i)) < 1d-12.or.(abs(dxi(i))<eps)) exit
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
