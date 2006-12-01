
c module grid_definition
c #####################################################################
      module grid_definition

        use grid_structures

        use math

        implicit none

        real(8)         :: gparams(ngparams)

        character*(3)   :: coords

        real(8)         :: xmax,ymax,zmax,xmin,ymin,zmin  !3D domain dimensions

        real(8),private :: pi,lambda,cc,xpp,ypp,zpp,eps,mm,kk,aa,phi
     .                    ,major_r,xi0,yi0,zi0,lambdax,lambday,lambdaz
     .                    ,aax=0d0,bbx=0d0,xx0
     .                    ,aay=0d0,bby=0d0,yy0
     .                    ,aaz=0d0,bbz=0d0,zz0

        logical         :: numerical_grid,anal_map

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

c     checkMapDatabase
c     #################################################################
      function checkMapDatabase() result(is_map)

c     -----------------------------------------------------------------
c     Checks grid database for analytical mappings
c     -----------------------------------------------------------------

        implicit none

c     Input variables

        logical :: is_map

c     Begin program

        select case(coords)
        case('car','scl','cyl','hel','hl2','tor','sin','vmc')

          is_map = .true.

        case default

          is_map = .false.

        end select

      end function checkMapDatabase

c     map
c     #################################################################
      function map(i,j,k,igx,igy,igz,ig,jg,kg) result(car)

c     -----------------------------------------------------------------
c     Give Cartesian coordinates corresponding to node (i,j,k) at grid
c     level (igx,igy,igz).
c     -----------------------------------------------------------------

        implicit none

c     Input variables

        integer(4) :: i,j,k,igx,igy,igz,ig,jg,kg
        real(8)    :: car(3)

c     Local variables

        real(8)    :: x1,y1,z1

c     External routines

c     Begin program

        call getCurvilinearCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
     .                                  ,x1,y1,z1)

        car = x_xi(x1,y1,z1)

      end function map

c     getCartesianCoordinates
c     #################################################################
      subroutine getCartesianCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
     .                                  ,x1,y1,z1)

c     -----------------------------------------------------------------
c     Inverts curvilinear coordinates to give Cartesian coordinates.
c     Requires external routine 'map', with call sequence:
c
c             map(i,j,k,igx,igy,igz,ig,jg,kg,x1,y1,z1)
c
c     -----------------------------------------------------------------

        implicit none

c     Input variables

        integer(4) :: i,j,k,igx,igy,igz,ig,jg,kg
        real(8)    :: x1,y1,z1

c     Local variables

        real(8)    :: car(3)

c     Begin program

        car = map(i,j,k,igx,igy,igz,ig,jg,kg)

        x1 = car(1)
        y1 = car(2)
        z1 = car(3)

      end subroutine getCartesianCoordinates

c     x_xi
c     #################################################################
      function x_xi(x1,x2,x3) result(car)

c     -----------------------------------------------------------------
c     Gives Cartesian coordinates from curvilinear coordinates
c     -----------------------------------------------------------------

        implicit none

c     Input variables

        real(8)    :: x1,x2,x3,car(3)

c     Local variables

        integer(4) :: inewt,ic
        real(8)    :: xx,yy,zz,jac_mat(3,3),rhs(3),dx(3),rr,rr0,dd

c     Begin program

        select case (coords)
        case ('car')
          xx = x1
          yy = x2
          zz = x3
        case ('scl')
          lambdax = gparams(1)
          xx0     = gparams(2)
          if (lambdax == 0d0) then
            xx = x1
          else
            if (aax == 0d0 .and. bbx == 0d0) then
              dd = 0.01*(xmax-xmin)
              aax = find_a(xmin,xmax,xx0,lambdax)
              bbx = aax/lambdax
              xi0 = xmax - bbx*tanh((xmax-xx0)/aax)
            endif

            if (x1 > xmax) then
              xx = 2*xmax -(xx0 + aax*atanh((2*xmax-x1-xi0)/bbx))
            elseif (x1 < xmin) then
              xx = 2*xmin -(xx0 + aax*atanh((2*xmin-x1-xi0)/bbx))
            else
              xx = xx0 + aax*atanh((x1-xi0)/bbx)
            endif
cc            write (*,*) x1,xx,aax,bbx,(x1-xi0)/bbx
          endif

          lambday = gparams(3)
          yy0     = gparams(4)
          if (lambday == 0d0) then
            yy = x2
          else
            if (aay == 0d0 .and. bby == 0d0) then
              aay = find_a(ymin,ymax,yy0,lambday)
              bby = aay/lambday
              yi0 = ymax - bby*tanh((ymax-yy0)/aay)
            endif

            if (x2 > ymax) then
              yy = 2*ymax -(yy0 + aay*atanh((2*ymax-x2-yi0)/bby))
            elseif (x2 < ymin) then
              yy = 2*ymin -(yy0 + aay*atanh((2*ymin-x2-yi0)/bby))
            else
              yy = yy0 + aay*atanh((x2-yi0)/bby)
            endif
          endif

          lambdaz = gparams(5)
          zz0     = gparams(6)
          if (lambdaz == 0d0) then
            zz = x3
          else
            if (aaz == 0d0 .and. bbz == 0d0) then
              aaz = find_a(zmin,zmax,zz0,lambdaz)
              bbz = aaz/lambdaz
              zi0 = zmax - bbz*tanh((zmax-zz0)/aaz)
            endif

            if (x3 > zmax) then
              zz = 2*zmax -(zz0 + aaz*atanh((2*zmax-x3-zi0)/bbz))
            elseif (x3 < zmin) then
              zz = 2*zmin -(zz0 + aaz*atanh((2*zmin-x3-zi0)/bbz))
            else
              zz = zz0 + aaz*atanh((x3-zi0)/bbz)
            endif
          endif
        case ('cyl')
          xx = x1*cos(x2)
          yy = x1*sin(x2)
          zz = x3
        case ('hel')
          mm = gparams(1)
          kk = gparams(2)
          aa = kk/mm
          phi = (x2-aa*x3)

          xx = x1*cos(phi)
          yy = x1*sin(phi)
          zz = x3
        case ('hl2')
          mm = gparams(1)
          kk = gparams(2)
          aa = kk/mm
          phi = (x2-aa*x3)

          xx = sign(sqrt(2*abs(x1)),x1)*cos(phi)
          yy = sign(sqrt(2*abs(x1)),x1)*sin(phi)
          zz = x3
        case ('tor')
          major_r = gparams(1)

          xx = (major_r + x1*sin(x2))*cos(x3)
          yy = (major_r + x1*sin(x2))*sin(x3)
          zz = x1*cos(x2)
        case ('sin')
          pi = acos(-1d0)
          eps = gparams(1)

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

c     find_a
c     #################################################################
      function find_a(xmin,xmax,x0,cc) result(aa)

c     Global variables

        real(8) :: xmin,xmax,x0,cc,aa

c     Local variables

        real(8) :: ll,l0m,l0p,dd,ee,ff,bb

        integer(4),parameter :: nit=20
        integer(4) :: i
        real(8)    :: ijac,dxi(nit),res(nit)

c     Begin program

        ll  = xmax-xmin
        l0m = x0  -xmin
        l0p = xmax - x0

        dd = cc*ll

        if (cc > 1d0) then
cc          call pstop('find_a in grid_mod.F',
cc     .               'Slope has to <= 1')
          write (*,*) 'Slope has to <= 1'
          write (*,*) 'Aborting...'
          stop
        endif

c     Initial guess

        aa = 0.5*ll

c     Newton iteration

        do i=1,nit
          ff = dd/aa

          res(i) = dd - aa*(tanh(l0m/aa) + tanh(l0p/aa))
          ijac = 1./(-tanh(l0m/aa)-tanh(l0p/aa)
     .              +l0m/aa/cosh(l0m/aa)**2
     .              +l0p/aa/cosh(l0p/aa)**2)

          dxi(i) = -res(i)*ijac
          aa = abs(aa + dxi(i))  !Want to find positive root

cc          write (*,*) abs(res(i)),dxi(i),ijac,aa

          if (abs(res(i)) < 1d-10 .or. (abs(dxi(i)) < 1d-15) ) exit
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

      end function find_a

      end function x_xi

c     jacobian
c     #################################################################
      function jacobian(i,j,k,igx,igy,igz) result(jac)

c     -----------------------------------------------------------------
c     Calculates Jacobian of curvilinear coordinate system
c     -----------------------------------------------------------------

        implicit none

c     Input variables

        integer(4) :: i,j,k,igx,igy,igz
        real(8)    :: jac

c     Local variables

        integer(4) :: ig,jg,kg
        real(8)    :: x1,x2,x3,car(3),curv(3)

c     Begin program

        select case (coords)
        case ('car')
          jac = 1d0
        case ('scl')
          call getCurvilinearCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
     .                                  ,x1,x2,x3)

          curv = (/ x1,x2,x3 /)

          lambda = gparams(1)

          cc = 0.5/lambda
          cc = 1./tanh(cc)
          ypp = (2*curv(2)/ymax-1.)
          jac = cc*lambda/(cc**2-ypp**2)
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

          major_r = gparams(1)

          jac = x1*(major_r + x1*sin(x2))
        case ('sin')
          pi = acos(-1d0)

          call getCurvilinearCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
     .                                  ,x1,x2,x3)

          eps = gparams(1)
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

        integer(4) :: comp,i,j,k,igx,igy,igz

c     Local variables

        integer(4) :: ig,jg,kg
        real(8)    :: x1,x2,x3,vec(3)

        real(8)    :: car(3),curv(3),jac,r

c     Begin program

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
          lambda = gparams(1)

          call getCurvilinearCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
     .                                  ,x1,x2,x3)

          curv = (/ x1,x2,x3 /)

          cc = 0.5/lambda
          cc = 1./tanh(cc)
          ypp = (2*curv(2)/ymax-1.)
          jac = cc*lambda/(cc**2-ypp**2)
          jac = 1./jac
          select case (comp)
            case (1)
              vec = (/ 1d0,0d0,0d0 /)
            case (2)
              vec = (/ 0d0,jac,0d0 /)
            case (3)
              vec = (/ 0d0,0d0,1d0 /)
          end select
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

          mm = gparams(1)
          kk = gparams(2)
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

          mm = gparams(1)
          kk = gparams(2)
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
          major_r = gparams(1)

          call getCurvilinearCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
     .                                  ,x1,x2,x3)
          curv = (/ x1,x2,x3 /)

          select case (comp)
            case (1)
              vec = (/ sin(curv(2))*cos(curv(3))
     .                ,sin(curv(2))*sin(curv(3))
     .                ,cos(curv(2)) /)
            case (2)
              vec = (/ cos(curv(2))*cos(curv(3))
     .                ,cos(curv(2))*sin(curv(3))
     .                ,-sin(curv(2))/)/curv(1)
            case (3)
              vec = (/ -sin(curv(3)),cos(curv(3)),0d0 /)
     .              /(major_r + curv(1)*sin(curv(2)))
          end select
        case ('sin')

          pi = acos(-1d0)

          call getCurvilinearCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
     .                                  ,x1,x2,x3)

          eps = gparams(1)

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

        integer(4) :: comp,i,j,k,igx,igy,igz

c     Local variables

        integer(4) :: ig,jg,kg
        real(8)    :: car(3),curv(3),jac
        real(8)    :: x1,x2,x3,vec(3)

c     Begin program

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
          call getCurvilinearCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
     .                                  ,x1,x2,x3)
          curv = (/ x1,x2,x3 /)

          lambda = gparams(1)
          cc = 0.5/lambda
          cc = 1./tanh(cc)
          ypp = (2*curv(2)/ymax-1.)
          jac = cc*lambda/(cc**2-ypp**2)
          jac = 1./jac
          select case (comp)
            case (1)
              vec = (/ jac,0d0,0d0 /)
            case (2)
              vec = (/ 0d0,1d0,0d0 /)
            case (3)
              vec = (/ 0d0,0d0,jac /)
          end select
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

          mm = gparams(1)
          kk = gparams(2)
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

          mm = gparams(1)
          kk = gparams(2)
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
          major_r = gparams(1)

          call getCurvilinearCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
     .                                  ,x1,x2,x3)
          curv = (/ x1,x2,x3 /)

          select case (comp)
            case (1)
              vec = (/ sin(curv(2))*cos(curv(3))
     .                ,sin(curv(2))*sin(curv(3))
     .                ,cos(curv(2)) /)
     .              /curv(1)/(major_r + curv(1)*sin(curv(2)))
            case (2)
              vec = (/ cos(curv(2))*cos(curv(3))
     .                ,cos(curv(2))*sin(curv(3))
     .                ,-sin(curv(2))/)
     .              /(major_r + curv(1)*sin(curv(2)))
            case (3)
              vec = (/ -sin(curv(3)),cos(curv(3)),0d0 /)/curv(1)
          end select
        case ('sin')

          pi = acos(-1d0)

          call getCurvilinearCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
     .                                  ,x1,x2,x3)

          eps = gparams(1)

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

        integer(4) :: l,i,j,k,igx,igy,igz

c     Local variables

        integer(4) :: ig,jg,kg
        real(8)    :: car(3),curv(3),vec(9)
        real(8)    :: x1,x2,x3,tensor(3,3)

c     Begin program

        select case (coords)
        case ('car')
          vec = (/ 0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0 /)
        case ('scl')
          call getCurvilinearCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
     .                                  ,x1,x2,x3)
          curv = (/ x1,x2,x3 /)

          lambda = gparams(1)
          cc = 0.5/lambda
          cc = 1./tanh(cc)
          ypp = (2*curv(2)/ymax-1.)

          select case (l)
            case (1)
              vec(1) = 0d0
              vec(2) = 0d0
              vec(3) = 0d0
              vec(4) = vec(2)
              vec(5) = 0d0
              vec(6) = 0d0
              vec(7) = vec(3)
              vec(8) = vec(6)
              vec(9) = 0d0
            case (2)
              vec(1) = 0d0
              vec(2) = 0d0
              vec(3) = 0d0
              vec(4) = vec(2)
              vec(5) = 2.*ypp/(ypp**2-cc**2)
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

          mm = gparams(1)
          kk = gparams(2)
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

          mm = gparams(1)
          kk = gparams(2)
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

          major_r = gparams(1)
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
              vec(9) = sin(curv(2))*(major_r + curv(1)*sin(curv(2)))
            case (2)
              vec(1) = 0d0
              vec(2) = -1./curv(1)
              vec(3) = 0d0
              vec(4) = vec(2)
              vec(5) = 0d0
              vec(6) = 0d0
              vec(7) = vec(3)
              vec(8) = vec(6)
              vec(9) = (major_r + curv(1)*sin(curv(2)))*cos(curv(2))
     .                 /curv(1)
            case (3)
              vec(1) =  0d0
              vec(2) =  0d0
              vec(3) =  -sin(curv(2))/(major_r + curv(1)*Sin(curv(2)))
              vec(4) = vec(2)
              vec(5) = 0d0
              vec(6) = -curv(1)*cos(curv(2))
     .                 /(major_r + curv(1)*Sin(curv(2)))
              vec(7) = vec(3)
              vec(8) = vec(6)
              vec(9) = 0d0
          end select

          vec = -vec

        case ('sin')

          pi = acos(-1d0)

          call getCurvilinearCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
     .                                  ,x1,x2,x3)
          curv = (/ x1,x2,x3 /)

          eps = gparams(1)

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

        integer(4) :: l,i,j,k,igx,igy,igz
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

        integer(4) :: k
        real(8)    :: x1,x2,x3,tensor(3,3)

c     Local variables

        real(8)    :: vec(9),car(3),eps

c     Begin program

        select case (coords)
        case ('car')
          vec = (/ 0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0 /)
        case ('scl')
          lambda = gparams(1)
          cc = 0.5/lambda
          cc = 1./tanh(cc)
          ypp = (2*x2/ymax-1.)
          select case (k)
            case (1)
              vec(1) = 0d0
              vec(2) = 0d0
              vec(3) = 0d0
              vec(4) = 2.*ypp/(ypp**2-cc**2)
              vec(5) = 0d0
              vec(6) = 0d0
              vec(7) = 0d0
              vec(8) = 0d0
              vec(9) = 0d0
            case (2)
              vec(1) = 0d0
              vec(2) = 0d0
              vec(3) = 0d0
              vec(4) = 0d0
              vec(5) = 0d0
              vec(6) = 0d0
              vec(7) = 0d0
              vec(8) = 0d0
              vec(9) = 0d0
            case (3)
              vec(1) = 0d0
              vec(2) = 0d0
              vec(3) = 0d0
              vec(4) = 0d0
              vec(5) = 0d0
              vec(6) = 2.*ypp/(ypp**2-cc**2)
              vec(7) = 0d0
              vec(8) = 0d0
              vec(9) = 0d0
          end select
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
          mm = gparams(1)
          kk = gparams(2)
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
          major_r = gparams(1)
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

        integer(4) :: i,j,k,igx,igy,igz
        real(8)    :: tensor(3,3)

c     Local variables

        integer(4) :: ig,jg,kg
        real(8)    :: vec(9),jac,car(3),curv(3)
        real(8)    :: x1,x2,x3

c     Begin program

        select case (coords)
        case ('car')
          vec = (/ 1d0, 0d0, 0d0
     .            ,0d0, 1d0, 0d0
     .            ,0d0, 0d0, 1d0 /)
        case ('scl')
          call getCurvilinearCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
     .                                  ,x1,x2,x3)
          curv = (/ x1,x2,x3 /)

          lambda = gparams(1)
          cc = 0.5/lambda
          cc = 1./tanh(cc)
          ypp = (2*curv(2)/ymax-1.)
          jac = cc*lambda/(cc**2-ypp**2)
          vec = (/ 1d0/jac, 0d0, 0d0
     .            ,0d0    , jac, 0d0
     .            ,0d0    , 0d0, 1d0/jac /)
        case ('cyl')
          call getCurvilinearCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
     .                                  ,x1,x2,x3)
          curv = (/ x1,x2,x3 /)

          vec = (/ 1d0/curv(1), 0d0     , 0d0
     .            ,0d0        , curv(1) , 0d0
     .            ,0d0        , 0d0     , 1d0/curv(1) /)
        case ('hel')
          call getCurvilinearCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
     .                                  ,x1,x2,x3)
          curv = (/ x1,x2,x3 /)

          mm = gparams(1)
          kk = gparams(2)
          aa = kk/mm
          vec = (/ 1./curv(1),0d0, 0d0
     .            ,0d0, curv(1) ,-aa*curv(1)
     .            ,0d0,-aa*curv(1), 1./curv(1) + aa**2*curv(1) /)
        case ('hl2')
          call getCurvilinearCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
     .                                  ,x1,x2,x3)
          curv = (/ abs(2*x1),x2,x3 /)

          mm = gparams(1)
          kk = gparams(2)
          aa = kk/mm
          vec = (/ 1./curv(1),0d0, 0d0
     .            ,0d0, curv(1) ,-aa*curv(1)
     .            ,0d0,-aa*curv(1), 1 + aa**2*curv(1) /)
        case ('tor')
          call getCurvilinearCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
     .                                  ,x1,x2,x3)
          curv = (/ x1,x2,x3 /)

          major_r = gparams(1)
          vec = (/ 1d0/curv(1)/(major_r + curv(1)*sin(curv(2))),0d0,0d0
     .            ,0d0, curv(1)/(major_r + curv(1)*sin(curv(2))), 0d0
     .            ,0d0, 0d0, (major_r/curv(1) + sin(curv(2))) /)
        case ('sin')

          pi = acos(-1d0)

          call getCurvilinearCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
     .                                  ,x1,x2,x3)

          eps = gparams(1)

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

        integer(4) :: i,j,k,igx,igy,igz
        real(8)    :: tensor(3,3)

c     Local variables

        integer(4) :: ig,jg,kg
        real(8)    :: vec(9),jac,car(3),curv(3)
        real(8)    :: x1,x2,x3

c     Begin program

        select case (coords)
        case ('car')
          vec = (/ 1d0, 0d0, 0d0
     .            ,0d0, 1d0, 0d0
     .            ,0d0, 0d0, 1d0 /)
       case ('scl')
          call getCurvilinearCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
     .                                  ,x1,x2,x3)
          curv = (/ x1,x2,x3 /)

          lambda = gparams(1)
          cc = 0.5/lambda
          cc = 1./tanh(cc)
          ypp = (2*curv(2)/ymax-1.)
          jac = cc*lambda/(cc**2-ypp**2)
          vec = (/ jac, 0d0    , 0d0
     .            ,0d0, 1d0/jac, 0d0
     .            ,0d0, 0d0    , jac /)
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

          mm = gparams(1)
          kk = gparams(2)
          aa = kk/mm
          vec = (/ curv(1),0d0,0d0
     .            ,0d0,1./curv(1) + aa**2*curv(1), aa*curv(1)
     .            ,0d0, aa*curv(1), curv(1) /)
        case ('hl2')
          call getCurvilinearCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
     .                                  ,x1,x2,x3)
          curv = (/ abs(2*x1),x2,x3 /)

          mm = gparams(1)
          kk = gparams(2)
          aa = kk/mm
          vec = (/ curv(1),0d0,0d0
     .            ,0d0,1./curv(1) + aa**2, aa
     .            ,0d0, aa, 1d0 /)
        case ('tor')
          call getCurvilinearCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
     .                                  ,x1,x2,x3)
          curv = (/ x1,x2,x3 /)

          major_r = gparams(1)
          vec = (/ curv(1)*(major_r + curv(1)*sin(curv(2))), 0d0, 0d0
     .            ,0d0, (major_r/curv(1) + sin(curv(2))), 0d0
     .            ,0d0, 0d0, curv(1)/(major_r + curv(1)*sin(curv(2))) /)
        case ('sin')

          pi = acos(-1d0)

          call getCurvilinearCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
     .                                  ,x1,x2,x3)

          eps = gparams(1)

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

      end module grid_definition
