
c module grid_anal_map
c #####################################################################
      module grid_anal_map

        use math

        use grid_def

        use grid_operations

        implicit none

        real(8) :: gparams(ngparams)

        logical :: anal_map

      contains

c     checkAnalMapDatabase
c     #################################################################
      function checkAnalMapDatabase() result(anal_map)

c     -----------------------------------------------------------------
c     Checks grid database for analytical mappings
c     -----------------------------------------------------------------

        implicit none

c     Input variables

        logical :: anal_map

c     Begin program

        select case(coords)
        case('car','scl','cyl','hel','hl2','tor','sin')

          anal_map = .true.

        case default

          anal_map = .false.

        end select

      end function checkAnalMapDatabase

c     map
c     #################################################################
      function map(i,j,k,igx,igy,igz,ig,jg,kg) result(car)

c     -----------------------------------------------------------------
c     Give Cartesian coordinates corresponding to node (i,j,k) at grid
c     level (igx,igy,igz).
c     -----------------------------------------------------------------

        implicit none

c     Input variables

        integer    :: i,j,k,igx,igy,igz
        real(8)    :: car(3)
        integer   ,optional,intent(OUT) :: ig,jg,kg

c     Local variables

        integer    :: ig1,jg1,kg1

        real(8)    :: x1,y1,z1

c     External routines

c     Begin program
        
        call getCurvilinearCoordinates(i,j,k,igx,igy,igz,ig1,jg1,kg1
     .                                ,x1,y1,z1)

        car = x_xi(x1,y1,z1)

        if (PRESENT(ig)) ig = ig1
        if (PRESENT(jg)) jg = jg1
        if (PRESENT(kg)) kg = kg1

      end function map

c     x_xi
c     #################################################################
      function x_xi(x1,x2,x3) result(car)

c     -----------------------------------------------------------------
c     Gives Cartesian coordinates from curvilinear coordinates
c     -----------------------------------------------------------------

        implicit none

c     Input variables

        real(8) :: x1,x2,x3,car(3)

c     Local variables

        integer :: inewt,ic,exp=3
        real(8) :: xx,yy,zz,jac_mat(3,3),rhs(3),dx(3),rr,rr0,dd
     .            ,aa,bb,major_r,mm,kk,eps,phi,gp(ngparams)

        real(8),save :: aax=0d0,aay=0d0,aaz=0d0
     .                 ,bbx=0d0,bby=0d0,bbz=0d0
     .                 ,lbx=0d0,lby=0d0,lbz=0d0

c     Begin program

        gp = grid_params%params

        select case (coords)
        case ('car')
          xx = x1
          yy = x2
          zz = x3
        case ('scl')
cc          !Find monotonic atanh profile with given slope (gp(2,4,6)
cc          !  for X, Y, Z) at given position (gp(1,3,5) for X, Y, Z).
cc          if (aax == 0d0) aax = find_a(xmin,xmax,gp(2),gp(1))
cc          if (aay == 0d0) aay = find_a(ymin,ymax,gp(4),gp(3))
cc          if (aaz == 0d0) aaz = find_a(zmin,zmax,gp(6),gp(5))
cc
cc          xx = scl_coord(x1,aax,xmin,xmax,gp(2),gp(1))
cc          yy = scl_coord(x2,aay,ymin,ymax,gp(4),gp(3))
cc          zz = scl_coord(x3,aaz,zmin,zmax,gp(6),gp(5))

          !Find monotonic cubic polynomial with given slope (gp(2,4,6)
          !  for X, Y, Z) at given position (gp(1,3,5) for X, Y, Z).
          if (bbx == 0d0) bbx=find_eta(xmin,xmax,gp(1),gp(2))
          if (bby == 0d0) bby=find_eta(ymin,ymax,gp(3),gp(4))
          if (bbz == 0d0) bbz=find_eta(zmin,zmax,gp(5),gp(6))

          xx = scl_coord(x1,bbx,xmin,xmax,gp(1),gp(2))
          yy = scl_coord(x2,bby,ymin,ymax,gp(3),gp(4))
          zz = scl_coord(x3,bbz,zmin,zmax,gp(5),gp(6))
        case ('cyl')
          xx = x1*cos(x2)
          yy = x1*sin(x2)
          zz = x3
        case ('hel')
          mm = gp(1)
          kk = gp(2)
          aa = kk/mm
          phi = (x2-aa*x3)

          xx = x1*cos(phi)
          yy = x1*sin(phi)
          zz = x3
        case ('hl2')
          mm = gp(1)
          kk = gp(2)
          aa = kk/mm
          phi = (x2-aa*x3)

          xx = sign(sqrt(2*abs(x1)),x1)*cos(phi)
          yy = sign(sqrt(2*abs(x1)),x1)*sin(phi)
          zz = x3
        case ('tor')
          major_r = gp(1)  !Major radius
          aa      = gp(2)  !Horizontal minor radius
          bb      = gp(3)  !Vertical   minor radius

          xx = (major_r + aa*x1*cos(x2))*cos(-x3)
          yy = (major_r + aa*x1*cos(x2))*sin(-x3)
          zz = bb*x1*sin(x2)
        case ('sin')
          pi = acos(-1d0)
          eps = gp(1)

          xx = x1 + eps*sin(2*pi*x1/xmax)*sin(2*pi*x2/ymax)
          yy = x2 + eps*sin(2*pi*x1/xmax)*sin(2*pi*x2/ymax)
          zz = x3
        case default
          write (*,*) 'Grid not implemented in x_xi'
          write (*,*) 'Aborting...'
          stop
        end select

        car = (/ xx,yy,zz /)

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

      end function x_xi

c     jacobian
c     #################################################################
      function jacobian(i,j,k,igx,igy,igz) result(jac)

c     -----------------------------------------------------------------
c     Calculates Jacobian of curvilinear coordinate system
c     -----------------------------------------------------------------

        implicit none

c     Input variables

        integer    :: i,j,k,igx,igy,igz
        real(8)    :: jac

c     Local variables

        integer    :: ig,jg,kg
        real(8)    :: x1,x2,x3,car(3),curv(3),aa,bb,major_r,mm,kk,eps
     .               ,gp(ngparams)

c     Begin program

        gp = grid_params%params

        select case (coords)
        case ('car')
          jac = 1d0
        case ('scl')
          jac = 0d0
cc          call getCurvilinearCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
cc     .                                  ,x1,x2,x3)
cc
cc          curv = (/ x1,x2,x3 /)
cc
cc          lambda = gp(1)
cc
cc          if (lambda /= 0d0) then
cc            cc = 0.5/lambda
cc            cc = 1./tanh(cc)
cc          else
cc            cc = 1d0
cc          endif
cc          ypp = (2*curv(2)/ymax-1.)
cc          jac = cc*lambda/(cc**2-ypp**2)
        case ('cyl')
          call getCurvilinearCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
     .                                  ,x1,x2,x3)
          jac = x1
        case ('hel')
          call getCurvilinearCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
     .                                  ,x1,x2,x3)
          jac = x1
        case ('hl2')
          jac = 1d0
        case ('tor')
          call getCurvilinearCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
     .                                  ,x1,x2,x3)

          major_r = gp(1)
          aa      = gp(2)
          bb      = gp(3)

          jac = aa*bb*x1*(major_r + aa*x1*cos(x2))
        case ('sin')
          pi = acos(-1d0)

          call getCurvilinearCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
     .                                  ,x1,x2,x3)

          eps = gp(1)
          jac = (xmax*ymax + eps*Pi*(xmax - ymax)*
     -          Sin(Pi*((2*x1)/xmax - (2*x2)/ymax)) + 
     -          eps*Pi*(xmax + ymax)*Sin(2*Pi*(x1/xmax + x2/ymax)))
     -          /(xmax*ymax)
        case default
          write (*,*) 'Grid not implemented in jacobian'
          write (*,*) 'Aborting...'
          stop
        end select

      end function jacobian

c     covariantVector
c     #################################################################
      function covariantVector(comp,i,j,k,igx,igy,igz) result (vec)

c     -----------------------------------------------------------------
c     Calculates covariant vectors of curvilinear coordinate system
c     in Cartesian coordinates
c     -----------------------------------------------------------------

        implicit none

c     Input variables

        integer    :: comp,i,j,k,igx,igy,igz

c     Local variables

        integer    :: ig,jg,kg
        real(8)    :: x1,x2,x3,aa,bb,major_r,mm,kk,eps,phi,jac,r
        real(8)    :: car(3),curv(3),vec(3),gp(ngparams)

c     Begin program

        gp = grid_params%params

        select case (coords)
        case ('car')
          select case (comp)
            case (1)
              vec = (/ 1d0,0d0,0d0 /)
            case (2)
              vec = (/ 0d0,1d0,0d0 /)
            case (3)
              vec = (/ 0d0,0d0,1d0 /)
          end select
        case ('scl')
          vec = 0d0
cc          lambda = gp(1)
cc
cc          call getCurvilinearCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
cc     .                                  ,x1,x2,x3)
cc
cc          curv = (/ x1,x2,x3 /)
cc
cc          cc = 0.5/lambda
cc          cc = 1./tanh(cc)
cc          ypp = (2*curv(2)/ymax-1.)
cc          jac = cc*lambda/(cc**2-ypp**2)
cc          jac = 1./jac
cc          select case (comp)
cc            case (1)
cc              vec = (/ 1d0,0d0,0d0 /)
cc            case (2)
cc              vec = (/ 0d0,jac,0d0 /)
cc            case (3)
cc              vec = (/ 0d0,0d0,1d0 /)
cc          end select
        case ('cyl')
          call getCurvilinearCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
     .                                  ,x1,x2,x3)
          curv = (/ x1,x2,x3 /)

          select case (comp)
            case (1)
              vec = (/ cos(curv(2)),sin(curv(2)),0d0 /)
            case (2)
              vec = (/-sin(curv(2)),cos(curv(2)),0d0 /)/curv(1)
            case (3)
              vec = (/ 0d0,0d0,1d0 /)
          end select
        case ('hel')
          call getCurvilinearCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
     .                                  ,x1,x2,x3)
          curv = (/ x1,x2,x3 /)

          mm = gp(1)
          kk = gp(2)
          aa = kk/mm
          phi = (curv(2)-aa*curv(3))

          select case (comp)
            case (1)
              vec = (/ cos(phi),sin(phi),0d0 /)
            case (2)
              vec = (/-sin(phi)/curv(1),cos(phi)/curv(1),aa /)
            case (3)
              vec = (/ 0d0,0d0,1d0 /)
          end select
        case ('hl2')
          call getCurvilinearCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
     .                                  ,x1,x2,x3)
          curv = (/ sign(sqrt(2*abs(x1)),x1),x2,x3 /)

          mm = gp(1)
          kk = gp(2)
          aa = kk/mm
          phi = (curv(2)-aa*curv(3))

          select case (comp)
            case (1)
              vec = (/ cos(phi),sin(phi),0d0 /)*curv(1)
            case (2)
              vec = (/-sin(phi)/curv(1),cos(phi)/curv(1),aa /)
            case (3)
              vec = (/ 0d0,0d0,1d0 /)
          end select
        case ('tor')
          major_r = gp(1)
          aa      = gp(2)
          bb      = gp(3)

          call getCurvilinearCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
     .                                  ,x1,x2,x3)
          curv = (/ x1,x2,x3 /)

          select case (comp)
            case (1)
              vec = (/  cos(curv(2))*cos(curv(3))/aa
     .                ,-cos(curv(2))*sin(curv(3))/aa
     .                , sin(curv(2))/bb             /)
            case (2)
              vec = (/-sin(curv(2))*cos(curv(3))/aa
     .                ,sin(curv(2))*sin(curv(3))/aa
     .                ,cos(curv(2))/bb             /)/curv(1)
            case (3)
              vec = (/ -sin(curv(3)),-cos(curv(3)),0d0 /)
     .              /(major_r + aa*curv(1)*cos(curv(2)))
          end select
        case ('sin')

          pi = acos(-1d0)

          call getCurvilinearCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
     .                                  ,x1,x2,x3)

          eps = gp(1)

          select case (comp)
            case (1)
              vec = (/ (xmax*(ymax + 2*eps*Pi*Cos((2*Pi*x2)/ymax)*
     -                                        Sin((2*Pi*x1)/xmax)))/
     -                 (xmax*ymax + eps*Pi*(xmax - ymax)
     -                         *Sin(Pi*((2*x1)/xmax - (2*x2)/ymax))
     -                            + eps*Pi*(xmax + ymax)
     .                         *Sin(2*Pi*(x1/xmax + x2/ymax)))
     .              , (-2*eps*Pi*xmax*Cos((2*Pi*x2)/ymax)
     .                               *Sin((2*Pi*x1)/xmax))/
     -                (xmax*ymax + eps*Pi*(xmax - ymax)
     -                         *Sin(Pi*((2*x1)/xmax - (2*x2)/ymax))  
     -                           + eps*Pi*(xmax + ymax)
     .                         *Sin(2*Pi*(x1/xmax + x2/ymax)))
     .              , 0d0 /)
            case (2)
              vec = (/ (-2*eps*Pi*ymax*Cos((2*Pi*x1)/xmax)
     .                                *Sin((2*Pi*x2)/ymax))/
     -                 (xmax*ymax + eps*Pi*(xmax - ymax)
     -                         *Sin(Pi*((2*x1)/xmax - (2*x2)/ymax))
     -                            + eps*Pi*(xmax + ymax)
     .                         *Sin(2*Pi*(x1/xmax + x2/ymax)))
     .              , (ymax*(xmax + 2*eps*Pi*Cos((2*Pi*x1)/xmax)
     .                                      *Sin((2*Pi*x2)/ymax)))/
     -                (xmax*ymax + eps*Pi*(xmax - ymax)
     -                         *Sin(Pi*((2*x1)/xmax - (2*x2)/ymax))
     -                           + eps*Pi*(xmax + ymax)
     .                         *Sin(2*Pi*(x1/xmax + x2/ymax)))
     .              ,0d0 /)
            case (3)
              vec = (/ 0d0,0d0,1d0 /)
          end select
        case default
          write (*,*) 'Grid not implemented in covariantVector'
          write (*,*) 'Aborting...'
          stop
        end select

      end function covariantVector

c     contravariantVector
c     #################################################################
      function contravariantVector(comp,i,j,k,igx,igy,igz) result (vec)

c     -----------------------------------------------------------------
c     Calculates contravariant vectors of curvilinear coordinate system
c     in Cartesian coordinates
c     -----------------------------------------------------------------

        implicit none

c     Input variables

        integer    :: comp,i,j,k,igx,igy,igz
        real(8)    :: vec(3)

c     Local variables

        integer    :: ig,jg,kg
        real(8)    :: car(3),curv(3),jac
        real(8)    :: x1,x2,x3,aa,bb,major_r,mm,kk,eps,phi,gp(ngparams)

c     Begin program

        gp = grid_params%params

        select case (coords)
        case ('car')
          select case (comp)
            case (1)
              vec = (/ 1d0,0d0,0d0 /)
            case (2)
              vec = (/ 0d0,1d0,0d0 /)
            case (3)
              vec = (/ 0d0,0d0,1d0 /)
          end select
       case ('scl')
         vec = 0d0
cc          call getCurvilinearCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
cc     .                                  ,x1,x2,x3)
cc          curv = (/ x1,x2,x3 /)
cc
cc          lambda = gp(1)
cc          cc = 0.5/lambda
cc          cc = 1./tanh(cc)
cc          ypp = (2*curv(2)/ymax-1.)
cc          jac = cc*lambda/(cc**2-ypp**2)
cc          jac = 1./jac
cc          select case (comp)
cc            case (1)
cc              vec = (/ jac,0d0,0d0 /)
cc            case (2)
cc              vec = (/ 0d0,1d0,0d0 /)
cc            case (3)
cc              vec = (/ 0d0,0d0,jac /)
cc          end select
        case ('cyl')
          call getCurvilinearCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
     .                                  ,x1,x2,x3)
          curv = (/ x1,x2,x3 /)

          select case (comp)
            case (1)
              vec = (/ cos(curv(2)),sin(curv(2)),0d0 /)/curv(1)
            case (2)
              vec = (/-sin(curv(2)),cos(curv(2)),0d0 /)
            case (3)
              vec = (/ 0d0,0d0,1d0 /)/curv(1)
          end select
        case ('hel')
          call getCurvilinearCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
     .                                  ,x1,x2,x3)
          curv = (/ x1,x2,x3 /)

          mm = gp(1)
          kk = gp(2)
          aa = kk/mm
          phi = (curv(2)-aa*curv(3))
          select case (comp)
            case (1)
              vec = (/ cos(phi),sin(phi),0d0 /)/curv(1)
            case (2)
              vec = (/-sin(phi),cos(phi),0d0 /)
            case (3)
              vec = (/ aa*sin(phi),-aa*cos(phi),1d0/curv(1) /)
          end select
        case ('hl2')
          call getCurvilinearCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
     .                                  ,x1,x2,x3)
          curv = (/ sign(sqrt(2*abs(x1)),x1),x2,x3 /)

          mm = gp(1)
          kk = gp(2)
          aa = kk/mm
          phi = (curv(2)-aa*curv(3))
          select case (comp)
            case (1)
              vec = (/ cos(phi),sin(phi),0d0 /)/curv(1)
            case (2)
              vec = (/-sin(phi),cos(phi),0d0 /)*curv(1)
            case (3)
              vec = (/ aa*sin(phi)*curv(1),-aa*cos(phi)*curv(1),1d0 /)
          end select
        case ('tor')
          major_r = gp(1)
          aa      = gp(2)
          bb      = gp(3)

          call getCurvilinearCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
     .                                  ,x1,x2,x3)
          curv = (/ x1,x2,x3 /)

          select case (comp)
            case (1)
              vec = (/  cos(curv(2))*cos(curv(3))/bb
     .                ,-cos(curv(2))*sin(curv(3))/bb
     .                , sin(curv(2))/aa             /)
     .              /curv(1)/(major_r + curv(1)*cos(curv(2)))
            case (2)
              vec = (/-sin(curv(2))*cos(curv(3))/bb
     .                ,sin(curv(2))*sin(curv(3))/bb
     .                ,cos(curv(2))/aa              /)
     .              /(major_r + curv(1)*cos(curv(2)))
            case (3)
              vec = (/ -sin(curv(3)),-cos(curv(3)),0d0 /)/curv(1)/aa/bb
          end select
        case ('sin')

          pi = acos(-1d0)

          call getCurvilinearCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
     .                                  ,x1,x2,x3)

          eps = gp(1)

          select case (comp)
            case (1)
              vec = (/ (ymax*(xmax + 2*eps*Pi*Cos((2*Pi*x1)/xmax)*
     -                                        Sin((2*Pi*x2)/ymax)))/
     -                 (xmax*ymax + eps*Pi*(xmax - ymax)
     -                         *Sin(Pi*((2*x1)/xmax - (2*x2)/ymax))
     -                            + eps*Pi*(xmax + ymax)
     .                         *Sin(2*Pi*(x1/xmax + x2/ymax)))
     .              , (2*eps*Pi*ymax*Cos((2*Pi*x1)/xmax)
     .                              *Sin((2*Pi*x2)/ymax))/
     -                (xmax*ymax + eps*Pi*(xmax - ymax)
     -                         *Sin(Pi*((2*x1)/xmax - (2*x2)/ymax))
     -                           + eps*Pi*(xmax + ymax)
     .                         *Sin(2*Pi*(x1/xmax + x2/ymax)))
     .              ,0d0 /)
            case (2)
              vec = (/ (2*eps*Pi*xmax*Cos((2*Pi*x2)/ymax)
     .                               *Sin((2*Pi*x1)/xmax))/
     -                 (xmax*ymax + eps*Pi*(xmax - ymax)
     -                         *Sin(Pi*((2*x1)/xmax - (2*x2)/ymax))
     -                            + eps*Pi*(xmax + ymax)
     .                         *Sin(2*Pi*(x1/xmax + x2/ymax)))
     .              ,(xmax*(ymax + 2*eps*Pi*Cos((2*Pi*x2)/ymax)
     .                                     *Sin((2*Pi*x1)/xmax)))/
     -               (xmax*ymax + eps*Pi*(xmax - ymax)
     -                         *Sin(Pi*((2*x1)/xmax - (2*x2)/ymax))
     -                          + eps*Pi*(xmax + ymax)
     .                         *Sin(2*Pi*(x1/xmax + x2/ymax)))
     .              ,0d0 /)
            case (3)
              vec = (/ 0d0
     .                ,0d0
     .                ,(xmax*ymax)/(xmax*ymax
     .                  + eps*Pi*(xmax - ymax)
     -                       *Sin(Pi*((2*x1)/xmax - (2*x2)/ymax))
     -                  + eps*Pi*(xmax + ymax)
     .                       *Sin(2*Pi*(x1/xmax + x2/ymax))) /)
          end select
        case default
          write (*,*) 'Grid not implemented in contravariantVector'
          write (*,*) 'Aborting...'
          stop
        end select

      end function contravariantVector

c     hessian22
c     #################################################################
      function hessian22(l,i,j,k,igx,igy,igz) result (tensor)

c     -----------------------------------------------------------------
c     Calculates hessian elements of curvilinear coordinate system in
c     the covariant basis, i.e., 
c            hessian[l](i,j) = J^2 <cnv(i)|grad(cov[l])|cnv(j)>
c     -----------------------------------------------------------------

        implicit none

c     Input variables

        integer    :: l,i,j,k,igx,igy,igz

c     Local variables

        integer    :: ig,jg,kg
        real(8)    :: car(3),curv(3),vec(9)
        real(8)    :: x1,x2,x3,tensor(3,3),aa,bb,major_r,mm,kk,eps,phi
     .               ,gp(ngparams)

c     Begin program

        gp = grid_params%params

        select case (coords)
        case ('car')
          vec = (/ 0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0 /)
        case ('scl')
          vec = 0d0
cc          call getCurvilinearCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
cc     .                                  ,x1,x2,x3)
cc          curv = (/ x1,x2,x3 /)
cc
cc          lambda = gp(1)
cc          cc = 0.5/lambda
cc          cc = 1./tanh(cc)
cc          ypp = (2*curv(2)/ymax-1.)
cc
cc          select case (l)
cc            case (1)
cc              vec(1) = 0d0
cc              vec(2) = 0d0
cc              vec(3) = 0d0
cc              vec(4) = vec(2)
cc              vec(5) = 0d0
cc              vec(6) = 0d0
cc              vec(7) = vec(3)
cc              vec(8) = vec(6)
cc              vec(9) = 0d0
cc            case (2)
cc              vec(1) = 0d0
cc              vec(2) = 0d0
cc              vec(3) = 0d0
cc              vec(4) = vec(2)
cc              vec(5) = 2.*ypp/(ypp**2-cc**2)
cc              vec(6) = 0d0
cc              vec(7) = vec(3)
cc              vec(8) = vec(6)
cc              vec(9) = 0d0
cc            case (3)
cc              vec(1) = 0d0
cc              vec(2) = 0d0
cc              vec(3) = 0d0
cc              vec(4) = vec(2)
cc              vec(5) = 0d0
cc              vec(6) = 0d0
cc              vec(7) = vec(3)
cc              vec(8) = vec(6)
cc              vec(9) = 0d0
cc          end select
cc
cc          vec = -vec

        case ('cyl')
          call getCurvilinearCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
     .                                  ,x1,x2,x3)
          curv = (/ x1,x2,x3 /)

          select case (l)
            case (1)
              vec(1) = 0d0
              vec(2) = 0d0
              vec(3) = 0d0
              vec(4) = vec(2)
              vec(5) = curv(1)
              vec(6) = 0d0
              vec(7) = vec(3)
              vec(8) = vec(6)
              vec(9) = 0d0
            case (2)
              vec(1) = 0d0
              vec(2) = -1./curv(1)
              vec(3) = 0d0
              vec(4) = vec(2)
              vec(5) = 0d0
              vec(6) = 0d0
              vec(7) = vec(3)
              vec(8) = vec(6)
              vec(9) = 0d0
            case (3)
              vec(1) = 0d0
              vec(2) = 0d0
              vec(3) = 0d0
              vec(4) = vec(2)
              vec(5) = 0d0
              vec(6) = 0d0
              vec(7) = vec(3)
              vec(8) = vec(6)
              vec(9) = 0d0
          end select

          vec = -vec

        case ('hel')
          call getCurvilinearCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
     .                                  ,x1,x2,x3)
          curv = (/ x1,x2,x3 /)

          mm = gp(1)
          kk = gp(2)
          aa = kk/mm
          phi = (curv(2)-aa*curv(3))
          select case (l)
            case (1)
              vec(1) = 0d0
              vec(2) = 0d0
              vec(3) = 0d0
              vec(4) = vec(2)
              vec(5) = curv(1)
              vec(6) = -aa*curv(1)
              vec(7) = vec(3)
              vec(8) = vec(6)
              vec(9) = aa**2*curv(1)
            case (2)
              vec(1) = 0d0
              vec(2) = -1./curv(1)
              vec(3) = aa/curv(1)
              vec(4) = vec(2)
              vec(5) = 0d0
              vec(6) = 0d0
              vec(7) = vec(3)
              vec(8) = vec(6)
              vec(9) = 0d0
            case (3)
              vec(1) = 0d0
              vec(2) = 0d0
              vec(3) = 0d0
              vec(4) = vec(2)
              vec(5) = 0d0
              vec(6) = 0d0
              vec(7) = vec(3)
              vec(8) = vec(6)
              vec(9) = 0d0
          end select

          vec = -vec

        case ('hl2')
          call getCurvilinearCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
     .                                  ,x1,x2,x3)
          curv = (/ abs(2*x1),x2,x3 /)

          mm = gp(1)
          kk = gp(2)
          aa = kk/mm
          phi = (curv(2)-aa*curv(3))
          select case (l)
            case (1)
              vec(1) = -1d0/curv(1)
              vec(2) = 0d0
              vec(3) = 0d0
              vec(4) = vec(2)
              vec(5) = -curv(1)
              vec(6) = aa*curv(1)
              vec(7) = vec(3)
              vec(8) = vec(6)
              vec(9) = -aa**2*curv(1)
            case (2)
              vec(1) = 0d0
              vec(2) = 1./curv(1)
              vec(3) = -aa/curv(1)
              vec(4) = vec(2)
              vec(5) = 0d0
              vec(6) = 0d0
              vec(7) = vec(3)
              vec(8) = vec(6)
              vec(9) = 0d0
            case (3)
              vec(1) = 0d0
              vec(2) = 0d0
              vec(3) = 0d0
              vec(4) = vec(2)
              vec(5) = 0d0
              vec(6) = 0d0
              vec(7) = vec(3)
              vec(8) = vec(6)
              vec(9) = 0d0
          end select

        case ('tor')
          call getCurvilinearCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
     .                                  ,x1,x2,x3)
          curv = (/ x1,x2,x3 /)

          major_r = gp(1)
          aa      = gp(2)
          bb      = gp(3)

          select case (l)
            case (1)
              vec(1) = 0d0
              vec(2) = 0d0
              vec(3) = 0d0
              vec(4) = vec(2)
              vec(5) =-curv(1)
              vec(6) = 0d0
              vec(7) = vec(3)
              vec(8) = vec(6)
              vec(9) =-cos(curv(2))*(major_r+aa*curv(1)*cos(curv(2)))/aa
            case (2)
              vec(1) = 0d0
              vec(2) = 1./curv(1)
              vec(3) = 0d0
              vec(4) = vec(2)
              vec(5) = 0d0
              vec(6) = 0d0
              vec(7) = vec(3)
              vec(8) = vec(6)
              vec(9) = (major_r + aa*curv(1)*cos(curv(2)))*sin(curv(2))
     .                 /curv(1)/aa
            case (3)
              vec(1) = 0d0
              vec(2) = 0d0
              vec(3) = aa*cos(curv(2))/(major_r+aa*curv(1)*cos(curv(2)))
              vec(4) = vec(2)
              vec(5) = 0d0
              vec(6) = -aa*curv(1)*sin(curv(2))
     .                 /(major_r + aa*curv(1)*cos(curv(2)))
              vec(7) = vec(3)
              vec(8) = vec(6)
              vec(9) = 0d0
          end select

        case ('sin')

          pi = acos(-1d0)

          call getCurvilinearCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
     .                                  ,x1,x2,x3)
          curv = (/ x1,x2,x3 /)

          eps = gp(1)

          select case (l)
            case (1)
              vec(1) = (-4*eps*Pi**2*ymax*Sin((2*Pi*x1)/xmax)*
     -                                    Sin((2*Pi*x2)/ymax))/
     -          (xmax*(xmax*ymax + eps*Pi*(xmax - ymax)*
     -           Sin(Pi*((2*x1)/xmax - (2*x2)/ymax)) + 
     -           eps*Pi*(xmax + ymax)*Sin(2*Pi*(x1/xmax + x2/ymax))))
              vec(2) = (4*eps*Pi**2*Cos((2*Pi*x1)/xmax)
     .                             *Cos((2*Pi*x2)/ymax))/
     -            (xmax*ymax + eps*Pi*(xmax - ymax)*
     -              Sin(Pi*((2*x1)/xmax - (2*x2)/ymax)) + 
     -            eps*Pi*(xmax + ymax)*Sin(2*Pi*(x1/xmax + x2/ymax)))
              vec(3) = 0d0
              vec(4) = vec(2)
              vec(5) = (-4*eps*Pi**2*xmax*Sin((2*Pi*x1)/xmax)
     .                                   *Sin((2*Pi*x2)/ymax))/
     -            (ymax*(xmax*ymax + eps*Pi*(xmax - ymax)*
     -            Sin(Pi*((2*x1)/xmax - (2*x2)/ymax)) + 
     -            eps*Pi*(xmax + ymax)*Sin(2*Pi*(x1/xmax + x2/ymax))))
              vec(6) = 0d0
              vec(7) = vec(3)
              vec(8) = vec(6)
              vec(9) = 0d0
            case (2)
              vec(1) = (-4*eps*Pi**2*ymax*Sin((2*Pi*x1)/xmax)
     -                                   *Sin((2*Pi*x2)/ymax))/
     -                 (xmax*(xmax*ymax
     .                 + eps*Pi*(xmax - ymax)
     -                     *Sin(Pi*((2*x1)/xmax - (2*x2)/ymax))
     -                 + eps*Pi*(xmax + ymax)
     .                     *Sin(2*Pi*(x1/xmax + x2/ymax))))
              vec(2) = (4*eps*Pi**2*Cos((2*Pi*x1)/xmax)
     .                             *Cos((2*Pi*x2)/ymax))/
     -                 (xmax*ymax
     .                 + eps*Pi*(xmax - ymax)
     -                     *Sin(Pi*((2*x1)/xmax - (2*x2)/ymax))
     -                 + eps*Pi*(xmax + ymax)
     .                     *Sin(2*Pi*(x1/xmax + x2/ymax)))
              vec(3) = 0d0
              vec(4) = vec(2)
              vec(5) = (-4*eps*Pi**2*xmax*Sin((2*Pi*x1)/xmax)
     .                                   *Sin((2*Pi*x2)/ymax))/
     -                  (ymax*(xmax*ymax
     .                  + eps*Pi*(xmax - ymax)
     -                      *Sin(Pi*((2*x1)/xmax - (2*x2)/ymax))
     -                  + eps*Pi*(xmax + ymax)
     .                      *Sin(2*Pi*(x1/xmax + x2/ymax))))
              vec(6) = 0d0
              vec(7) = vec(3)
              vec(8) = vec(6)
              vec(9) = 0d0
            case (3)
              vec(1) = 0d0
              vec(2) = 0d0
              vec(3) = 0d0
              vec(4) = vec(2)
              vec(5) = 0d0
              vec(6) = 0d0
              vec(7) = vec(3)
              vec(8) = vec(6)
              vec(9) = 0d0
          end select

        case default
          write (*,*) 'Grid not implemented in hessian'
          write (*,*) 'Aborting...'
          stop
        end select

        tensor = reshape(vec,(/3,3/))

      end function hessian22

c     christ_2knd
c     #################################################################
      function christ_2knd(i,j,k,igx,igy,igz) result (tensor)

c     -----------------------------------------------------------------
c     Calculates elements of Christoffel symbol of the second kind,
c     Gamma[l](i,j) defined as: 
c            gamma[l](i,j) = -J^2 <cnv(i)|grad(cov[l])|cnv(j)>
c     -----------------------------------------------------------------

        implicit none

c     Input variables

        integer    :: l,i,j,k,igx,igy,igz
        real(8)    :: tensor(3,3,3)

c     Local variables

c     Begin program

        tensor(1,:,:) = hessian22(1,i,j,k,igx,igy,igz)
        tensor(2,:,:) = hessian22(2,i,j,k,igx,igy,igz)
        tensor(3,:,:) = hessian22(3,i,j,k,igx,igy,igz)

c     End program

      end function christ_2knd

c     hessian_cnv
c     #################################################################
      function hessian_cnv(k,x1,x2,x3) result (tensor)

c     -----------------------------------------------------------------
c     Calculates elements of tensor grad(cnv) in a mixed coordinate 
c     system:
c              hessian_cnv[k](i,j) = J^2 <cnv(i)|grad(cnv[k])|cov(j)>
c     -----------------------------------------------------------------

        implicit none

c     Input variables

        integer    :: k
        real(8)    :: x1,x2,x3,tensor(3,3)

c     Local variables

        real(8)    :: vec(9),car(3),aa,bb,major_r,mm,kk,eps,gp(ngparams)

c     Begin program

        gp = grid_params%params

        select case (coords)
        case ('car')
          vec = (/ 0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0 /)
        case ('scl')
          vec = 0d0
cc          lambda = gp(1)
cc          cc = 0.5/lambda
cc          cc = 1./tanh(cc)
cc          ypp = (2*x2/ymax-1.)
cc          select case (k)
cc            case (1)
cc              vec(1) = 0d0
cc              vec(2) = 0d0
cc              vec(3) = 0d0
cc              vec(4) = 2.*ypp/(ypp**2-cc**2)
cc              vec(5) = 0d0
cc              vec(6) = 0d0
cc              vec(7) = 0d0
cc              vec(8) = 0d0
cc              vec(9) = 0d0
cc            case (2)
cc              vec(1) = 0d0
cc              vec(2) = 0d0
cc              vec(3) = 0d0
cc              vec(4) = 0d0
cc              vec(5) = 0d0
cc              vec(6) = 0d0
cc              vec(7) = 0d0
cc              vec(8) = 0d0
cc              vec(9) = 0d0
cc            case (3)
cc              vec(1) = 0d0
cc              vec(2) = 0d0
cc              vec(3) = 0d0
cc              vec(4) = 0d0
cc              vec(5) = 0d0
cc              vec(6) = 2.*ypp/(ypp**2-cc**2)
cc              vec(7) = 0d0
cc              vec(8) = 0d0
cc              vec(9) = 0d0
cc          end select
        case ('cyl')
          select case (k)
            case (1)
              vec(1) = -1./x1
              vec(2) = 0d0
              vec(3) = 0d0
              vec(4) = 0d0
              vec(5) = 1./x1
              vec(6) = 0d0
              vec(7) = 0d0
              vec(8) = 0d0
              vec(9) = 0d0
            case (2)
              vec(1) = 0d0
              vec(2) = 0d0
              vec(3) = 0d0
              vec(4) = -x1
              vec(5) = 0d0
              vec(6) = 0d0
              vec(7) = 0d0
              vec(8) = 0d0
              vec(9) = 0d0
            case (3)
              vec(1) = 0d0
              vec(2) = 0d0
              vec(3) = -1./x1
              vec(4) = 0d0
              vec(5) = 0d0
              vec(6) = 0d0
              vec(7) = 0d0
              vec(8) = 0d0
              vec(9) = 0d0
          end select
        case ('hel')
          mm = gp(1)
          kk = gp(2)
          aa = kk/mm 
          select case (k)
          case (1)
            vec(1) = -1./x1
            vec(2) = 0d0
            vec(3) = 0d0
            vec(4) = 0d0
            vec(5) = 1./x1
            vec(6) = 0d0
            vec(7) = 0d0
            vec(8) = -aa/x1
            vec(9) = 0d0
          case (2)
            vec(1) = 0d0
            vec(2) = 0d0
            vec(3) = 0d0
            vec(4) = -x1
            vec(5) = 0d0
            vec(6) = 0d0
            vec(7) = aa*x1
            vec(8) = 0d0
            vec(9) = 0d0
          case (3)
            vec(1) = 0d0
            vec(2) = -aa/x1
            vec(3) = -1./x1
            vec(4) = aa*x1
            vec(5) = 0d0
            vec(6) = 0d0
            vec(7) = -aa**2*x1
            vec(8) = 0d0
            vec(9) = 0d0
          end select
        case ('tor')
          major_r = gp(1)
          select case (k)
            case (1)
              vec(1) = -(major_r + 2*x1*sin(x2))
     .                  /x1/(major_r + x1*sin(x2))
              vec(2) = 0d0
              vec(3) = 0d0
              vec(4) = -x1*cos(x2)/(major_r + x1*sin(x2))
              vec(5) = 1./x1
              vec(6) = 0d0
              vec(7) = 0d0
              vec(8) = 0d0
              vec(9) = sin(x2)/(major_r + x1*sin(x2))
            case (2)
              vec(1) = 0d0
              vec(2) = -sin(x2)/(major_r + x1*Sin(x2))
              vec(3) = 0d0
              vec(4) = -x1
              vec(5) = -x1*cos(x2)/(major_r + x1*Sin(x2))
              vec(6) = 0d0
              vec(7) = 0d0
              vec(8) = 0d0
              vec(9) =  x1*cos(x2)/(major_r + x1*Sin(x2))
            case (3)
              vec(1) =  0d0
              vec(2) =  0d0
              vec(3) =  -1./x1
              vec(4) = 0d0
              vec(5) = 0d0
              vec(6) = 0d0
              vec(7) = -sin(x2)*(major_r + x1*Sin(x2))
              vec(8) = -cos(x2)*(major_r + x1*Sin(x2))/x1
              vec(9) = 0d0
          end select

        case default
          write (*,*) 'Grid not implemented in hessian_cnv'
          write (*,*) 'Aborting...'
          stop
        end select

        tensor = transpose(reshape(vec, (/3,3/)))

      end function hessian_cnv

c     g_sub
c     #################################################################
      function g_sub(i,j,k,igx,igy,igz) result (tensor)

c     -----------------------------------------------------------------
c     Calculates contravariant metric tensor of curvilinear coordinate 
c     system
c     -----------------------------------------------------------------

        implicit none

c     Input variables

        integer    :: i,j,k,igx,igy,igz
        real(8)    :: tensor(3,3)

c     Local variables

        integer    :: ig,jg,kg
        real(8)    :: vec(9),jac,car(3),curv(3)
        real(8)    :: x1,x2,x3,aa,bb,major_r,mm,kk,eps,gp(ngparams)

c     Begin program

        gp = grid_params%params

        select case (coords)
        case ('car')
          vec = (/ 1d0, 0d0, 0d0
     .            ,0d0, 1d0, 0d0
     .            ,0d0, 0d0, 1d0 /)
        case ('scl')
          vec = 0d0
cc          call getCurvilinearCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
cc     .                                  ,curv(1),curv(2),curv(3))
cccc          curv = (/ x1,x2,x3 /)
cc
cc          lambda = gp(1)
cc          cc = 0.5/lambda
cc          cc = 1./tanh(cc)
cc          ypp = (2*curv(2)/ymax-1.)
cc          jac = cc*lambda/(cc**2-ypp**2)
cc          vec = (/ 1d0/jac, 0d0, 0d0
cc     .            ,0d0    , jac, 0d0
cc     .            ,0d0    , 0d0, 1d0/jac /)
        case ('cyl')
          call getCurvilinearCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
     .                                  ,x1,x2,x3)

          vec = (/ 1d0/x1, 0d0  , 0d0
     .            ,0d0   , x1   , 0d0
     .            ,0d0   , 0d0  , 1d0/x1 /)
        case ('hel')
          call getCurvilinearCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
     .                                  ,x1,x2,x3)

          mm = gp(1)
          kk = gp(2)
          aa = kk/mm
          vec = (/ 1./x1, 0d0  ,      0d0
     .            ,0d0  ,  x1  ,    -aa*x1
     .            ,0d0  ,-aa*x1, 1./x1 + aa**2*x1 /)
        case ('hl2')
          call getCurvilinearCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
     .                                  ,x1,x2,x3)
          curv = (/ abs(2*x1),x2,x3 /)

          mm = gp(1)
          kk = gp(2)
          aa = kk/mm
          vec = (/ 1./curv(1),0d0, 0d0
     .            ,0d0, curv(1) ,-aa*curv(1)
     .            ,0d0,-aa*curv(1), 1 + aa**2*curv(1) /)
        case ('tor')
          call getCurvilinearCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
     .                                  ,x1,x2,x3)

          major_r = gp(1)
          aa      = gp(2)
          bb      = gp(3)

cc          curv = (/ x1,x2,x3 /)
cc          vec = (/ 1d0/curv(1)/(major_r + curv(1)*cos(curv(2))),0d0,0d0
cc     .            ,0d0, curv(1)/(major_r + curv(1)*cos(curv(2))), 0d0
cc     .            ,0d0, 0d0, (major_r/curv(1) + cos(curv(2))) /)

          vec(1) = (aa**2 + bb**2 + (aa**2 - bb**2)*Cos(2*x2))/
     -             (2.*aa*bb*x1*(major_r + aa*x1*Cos(x2)))

          vec(2) = ((-aa**2 + bb**2)*Cos(x2)*Sin(x2))
     .             /(aa*bb*(major_r + aa*x1*Cos(x2)))

          vec(3) = 0d0

          vec(4) = vec(2)

          vec(5) = (x1*(aa**2 + bb**2 + (-aa**2 + bb**2)*Cos(2*x2)))/
     -             (2.*aa*bb*(major_r + aa*x1*Cos(x2)))

          vec(6) = 0d0

          vec(7) = vec(3)
          vec(8) = vec(6)
          vec(9) = (major_r + aa*x1*Cos(x2))/(aa*bb*x1)

        case ('sin')

          pi = acos(-1d0)

          call getCurvilinearCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
     .                                  ,x1,x2,x3)

          eps = gp(1)

          vec(1) =
     .         (ymax*(2*eps**2*Pi**2 + xmax**2 + 
     -        2*eps**2*Pi**2*Cos((4*Pi*x1)/xmax) - 
     -        eps**2*Pi**2*Cos(Pi*((4*x1)/xmax - (4*x2)/ymax)) - 
     -        eps**2*Pi**2*Cos(4*Pi*(x1/xmax + x2/ymax)) - 
     -        2*eps**2*Pi**2*Cos((4*Pi*x2)/ymax) - 
     -        2*eps*Pi*xmax*Sin(Pi*((2*x1)/xmax - (2*x2)/ymax)) + 
     -        2*eps*Pi*xmax*Sin(2*Pi*(x1/xmax + x2/ymax))))/
     -    (xmax*(xmax*ymax + eps*Pi*(xmax - ymax)*
     -         Sin(Pi*((2*x1)/xmax - (2*x2)/ymax)) + 
     -        eps*Pi*(xmax + ymax)*Sin(2*Pi*(x1/xmax + x2/ymax))))

          vec(2) =
     .         (eps*Pi*(eps*Pi*Cos(Pi*((4*x1)/xmax - (4*x2)/ymax)) - 
     -        eps*Pi*Cos(4*Pi*(x1/xmax + x2/ymax)) + 
     -        xmax*Sin(Pi*((2*x1)/xmax - (2*x2)/ymax)) - 
     -        ymax*Sin(Pi*((2*x1)/xmax - (2*x2)/ymax)) + 
     -        xmax*Sin(2*Pi*(x1/xmax + x2/ymax)) + 
     -        ymax*Sin(2*Pi*(x1/xmax + x2/ymax))))/
     -    (xmax*ymax + eps*Pi*(xmax - ymax)*
     -       Sin(Pi*((2*x1)/xmax - (2*x2)/ymax)) + 
     -      eps*Pi*(xmax + ymax)*Sin(2*Pi*(x1/xmax + x2/ymax)))

          vec(3) = 0d0
          vec(4) = vec(2)

          vec(5) =
     .         (xmax*(2*eps**2*Pi**2 + ymax**2 - 
     -        2*eps**2*Pi**2*Cos((4*Pi*x1)/xmax) - 
     -        eps**2*Pi**2*Cos(Pi*((4*x1)/xmax - (4*x2)/ymax)) - 
     -        eps**2*Pi**2*Cos(4*Pi*(x1/xmax + x2/ymax)) + 
     -        2*eps**2*Pi**2*Cos((4*Pi*x2)/ymax) + 
     -        2*eps*Pi*ymax*Sin(Pi*((2*x1)/xmax - (2*x2)/ymax)) + 
     -        2*eps*Pi*ymax*Sin(2*Pi*(x1/xmax + x2/ymax))))/
     -    (ymax*(xmax*ymax + eps*Pi*(xmax - ymax)*
     -         Sin(Pi*((2*x1)/xmax - (2*x2)/ymax)) + 
     -        eps*Pi*(xmax + ymax)*Sin(2*Pi*(x1/xmax + x2/ymax))))!

          vec(6) = 0d0
          vec(7) = vec(3)
          vec(8) = vec(6)

          vec(9) = (xmax*ymax)/
     -    (xmax*ymax + eps*Pi*(xmax - ymax)*
     -       Sin(Pi*((2*x1)/xmax - (2*x2)/ymax)) + 
     -      eps*Pi*(xmax + ymax)*Sin(2*Pi*(x1/xmax + x2/ymax)))

        case default
          write (*,*) 'Grid not implemented in g_sub'
          write (*,*) 'Aborting...'
          stop
        end select

        tensor = reshape(vec, (/3,3/))

      end function g_sub

c     g_sup
c     #################################################################
      function g_sup(i,j,k,igx,igy,igz) result (tensor)

c     -----------------------------------------------------------------
c     Calculates covariant metric tensor of curvilinear coordinate 
c     system
c     -----------------------------------------------------------------

        implicit none

c     Input variables

        integer    :: i,j,k,igx,igy,igz
        real(8)    :: tensor(3,3)

c     Local variables

        integer    :: ig,jg,kg
        real(8)    :: vec(9),jac,car(3),curv(3)
        real(8)    :: x1,x2,x3,aa,bb,major_r,mm,kk,eps,gp(ngparams)

c     Begin program

        gp = grid_params%params

        select case (coords)
        case ('car')
          vec = (/ 1d0, 0d0, 0d0
     .            ,0d0, 1d0, 0d0
     .            ,0d0, 0d0, 1d0 /)
       case ('scl')
          vec = 0d0
cc          call getCurvilinearCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
cc     .                                  ,x1,x2,x3)
cc          curv = (/ x1,x2,x3 /)
cc
cc          lambda = gp(1)
cc          cc = 0.5/lambda
cc          cc = 1./tanh(cc)
cc          ypp = (2*curv(2)/ymax-1.)
cc          jac = cc*lambda/(cc**2-ypp**2)
cc          vec = (/ jac, 0d0    , 0d0
cc     .            ,0d0, 1d0/jac, 0d0
cc     .            ,0d0, 0d0    , jac /)
        case ('cyl')
          call getCurvilinearCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
     .                                  ,x1,x2,x3)
          curv = (/ x1,x2,x3 /)

          vec = (/ curv(1) , 0d0   , 0d0
     .            ,0d0, 1d0/curv(1), 0d0
     .            ,0d0, 0d0   , curv(1) /)
        case ('hel')
          call getCurvilinearCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
     .                                  ,x1,x2,x3)
          curv = (/ x1,x2,x3 /)

          mm = gp(1)
          kk = gp(2)
          aa = kk/mm
          vec = (/ curv(1),0d0,0d0
     .            ,0d0,1./curv(1) + aa**2*curv(1), aa*curv(1)
     .            ,0d0, aa*curv(1), curv(1) /)
        case ('hl2')
          call getCurvilinearCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
     .                                  ,x1,x2,x3)
          curv = (/ abs(2*x1),x2,x3 /)

          mm = gp(1)
          kk = gp(2)
          aa = kk/mm
          vec = (/ curv(1),0d0,0d0
     .            ,0d0,1./curv(1) + aa**2, aa
     .            ,0d0, aa, 1d0 /)
        case ('tor')
          call getCurvilinearCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
     .                                  ,x1,x2,x3)
          curv = (/ x1,x2,x3 /)

          major_r = gp(1)
          aa      = gp(2)
          bb      = gp(3)

cc          vec = (/ curv(1)*(major_r + curv(1)*cos(curv(2))), 0d0, 0d0
cc     .            ,0d0, (major_r/curv(1) + cos(curv(2))), 0d0
cc     .            ,0d0, 0d0, curv(1)/(major_r + curv(1)*cos(curv(2))) /)

          vec(1) = -(x1*(major_r + aa*x1*Cos(x2))*
     .             (-aa**2 - bb**2 + (aa**2 - bb**2)*Cos(2*x2)))
     .             /(2.*aa*bb)

          vec(2) = ((aa**2 - bb**2)*Cos(x2)
     .             *(major_r + aa*x1*Cos(x2))*Sin(x2))
     .             /(aa*bb)

          vec(3) = 0d0

          vec(4) = vec(2)

          vec(5) = ((major_r + aa*x1*Cos(x2))
     .            *(aa**2 + bb**2 + (aa**2 - bb**2)*Cos(2*x2)))
     .            /(2.*aa*bb*x1)

          vec(6) = 0d0

          vec(7) = vec(3)
          vec(8) = vec(6)
          vec(9) = (aa*bb*x1)/(major_r + aa*x1*Cos(x2))

        case ('sin')

          pi = acos(-1d0)

          call getCurvilinearCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
     .                                  ,x1,x2,x3)

          eps = gp(1)

          vec(1) =
     .         (xmax*(2*eps**2*Pi**2 + ymax**2 - 
     -        2*eps**2*Pi**2*Cos((4*Pi*x1)/xmax) - 
     -        eps**2*Pi**2*Cos(Pi*((4*x1)/xmax - (4*x2)/ymax)) - 
     -        eps**2*Pi**2*Cos(4*Pi*(x1/xmax + x2/ymax)) + 
     -        2*eps**2*Pi**2*Cos((4*Pi*x2)/ymax) + 
     -        2*eps*Pi*ymax*Sin(Pi*((2*x1)/xmax - (2*x2)/ymax)) + 
     -        2*eps*Pi*ymax*Sin(2*Pi*(x1/xmax + x2/ymax))))/
     -    (ymax*(xmax*ymax + eps*Pi*(xmax - ymax)*
     -         Sin(Pi*((2*x1)/xmax - (2*x2)/ymax)) + 
     -        eps*Pi*(xmax + ymax)*Sin(2*Pi*(x1/xmax + x2/ymax))))

          vec(2) =
     .         -((eps*Pi*(eps*Pi*Cos(Pi*((4*x1)/xmax - (4*x2)/ymax)) - 
     -          eps*Pi*Cos(4*Pi*(x1/xmax + x2/ymax)) + 
     -          xmax*Sin(Pi*((2*x1)/xmax - (2*x2)/ymax)) - 
     -          ymax*Sin(Pi*((2*x1)/xmax - (2*x2)/ymax)) + 
     -          xmax*Sin(2*Pi*(x1/xmax + x2/ymax)) + 
     -          ymax*Sin(2*Pi*(x1/xmax + x2/ymax))))/
     -      (xmax*ymax + eps*Pi*(xmax - ymax)*
     -         Sin(Pi*((2*x1)/xmax - (2*x2)/ymax)) + 
     -        eps*Pi*(xmax + ymax)*Sin(2*Pi*(x1/xmax + x2/ymax))))

          vec(3) = 0d0
          vec(4) = vec(2)

          vec(5) =
     .         (ymax*(2*eps**2*Pi**2 + xmax**2 + 
     -        2*eps**2*Pi**2*Cos((4*Pi*x1)/xmax) - 
     -        eps**2*Pi**2*Cos(Pi*((4*x1)/xmax - (4*x2)/ymax)) - 
     -        eps**2*Pi**2*Cos(4*Pi*(x1/xmax + x2/ymax)) - 
     -        2*eps**2*Pi**2*Cos((4*Pi*x2)/ymax) - 
     -        2*eps*Pi*xmax*Sin(Pi*((2*x1)/xmax - (2*x2)/ymax)) + 
     -        2*eps*Pi*xmax*Sin(2*Pi*(x1/xmax + x2/ymax))))/
     -    (xmax*(xmax*ymax + eps*Pi*(xmax - ymax)*
     -         Sin(Pi*((2*x1)/xmax - (2*x2)/ymax)) + 
     -        eps*Pi*(xmax + ymax)*Sin(2*Pi*(x1/xmax + x2/ymax))))

          vec(6) = 0d0
          vec(7) = vec(3)
          vec(8) = vec(6)

          vec(9) = (xmax*ymax + eps*Pi*(xmax - ymax)*
     -       Sin(Pi*((2*x1)/xmax - (2*x2)/ymax)) + 
     -      eps*Pi*(xmax + ymax)*Sin(2*Pi*(x1/xmax + x2/ymax)))
     .         /(xmax*ymax)

        case default
          write (*,*) 'Grid not implemented in g_sup'
          write (*,*) 'Aborting...'
          stop
        end select

        tensor = reshape(vec, (/3,3/))

      end function g_sup

      end module grid_anal_map
