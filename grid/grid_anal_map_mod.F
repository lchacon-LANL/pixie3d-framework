
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
        case('car','cyl','hel','hl2','tor','sin','lpk')

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

        integer :: inewt,ic
        real(8) :: xx,yy,zz,jac_mat(3,3),rhs(3),dx(3),rr,rr0,dd
     .            ,aa,bb,major_r,mm,kk,eps,phi,gp(ngparams)
     .            ,xlog(3)

c     Begin program

        if (map_grid_pack) then
          xlog(1) = scale_coord(xmin,xmax,x1,1)
          xlog(2) = scale_coord(ymin,ymax,x2,2)
          xlog(3) = scale_coord(zmin,zmax,x3,3)
        else
          xlog = (/ x1, x2, x3 /)
        endif

        gp = grid_params%params

        select case (coords)
        case ('car','lpk')
          car = xlog
        case ('cyl')
          car(1) = xlog(1)*cos(xlog(2))
          car(2) = xlog(1)*sin(xlog(2))
          car(3) = xlog(3)
        case ('hel')
          mm = gp(1)
          kk = gp(2)
          aa = kk/mm
cc          phi = (x2-aa*x3)
cc          xx = x1*cos(phi)
cc          yy = x1*sin(phi)
cc          zz = x3
          phi = xlog(2)-aa*xlog(3)
          car(1) = xlog(1)*cos(phi)
          car(2) = xlog(1)*sin(phi)
          car(3) = xlog(3)
        case ('hl2')
          mm = gp(1)
          kk = gp(2)
          aa = kk/mm
cc          phi = (x2-aa*x3)
cc
cc          xx = sign(sqrt(2*abs(x1)),x1)*cos(phi)
cc          yy = sign(sqrt(2*abs(x1)),x1)*sin(phi)
cc          zz = x3
          phi = xlog(2)-aa*xlog(3)
          car(1) = sign(sqrt(2*abs(xlog(1))),xlog(1))*cos(phi)
          car(2) = sign(sqrt(2*abs(xlog(1))),xlog(1))*sin(phi)
          car(3) = xlog(3)
        case ('tor')
          major_r = gp(1)  !Major radius
          aa      = gp(2)  !Horizontal minor radius
          bb      = gp(3)  !Vertical   minor radius

cc          xx = (major_r + aa*x1*cos(x2))*cos(-x3)
cc          yy = (major_r + aa*x1*cos(x2))*sin(-x3)
cc          zz = bb*x1*sin(x2)
          car(1) = (major_r + aa*xlog(1)*cos(xlog(2)))*cos(-xlog(3))
          car(2) = (major_r + aa*xlog(1)*cos(xlog(2)))*sin(-xlog(3))
          car(3) = bb*xlog(1)*sin(xlog(2))
        case ('sin')
          pi = acos(-1d0)
          eps = gp(1)

cc          xx = x1 + eps*sin(2*pi*x1/xmax)*sin(2*pi*x2/ymax)
cc          yy = x2 + eps*sin(2*pi*x1/xmax)*sin(2*pi*x2/ymax)
cc          zz = x3
          car(1) = xlog(1) + eps*sin(2*pi*xlog(1)/xmax)
     .                          *sin(2*pi*xlog(2)/ymax)
          car(2) = xlog(2) + eps*sin(2*pi*xlog(1)/xmax)
     .                          *sin(2*pi*xlog(2)/ymax)
          car(3) = xlog(3)
        case default
          write (*,*) 'Grid not implemented in x_xi'
          write (*,*) 'Aborting...'
          stop
        end select

      end function x_xi

c     jacobian
c     #################################################################
      function jacobian(i,j,k,igx,igy,igz,vertex) result(jac)

c     -----------------------------------------------------------------
c     Calculates Jacobian of curvilinear coordinate system
c     -----------------------------------------------------------------

        implicit none

c     Input variables

        integer :: i,j,k,igx,igy,igz
        real(8) :: jac
        logical,optional :: vertex

c     Local variables

        integer :: ig,jg,kg
        real(8) :: x1,x2,x3,car(3),curv(3),aa,bb,major_r,mm,kk,eps
     .            ,gp(ngparams),x1p,x2p,x3p
        logical :: vrtx

c     Begin program

        if (PRESENT(vertex)) then
          vrtx = vertex
        else
          vrtx = .false.
        endif

        gp = grid_params%params

        call getCurvilinearCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
     .                                ,x1,x2,x3)

        if (vrtx) then
          call getCurvilinearCoordinates(i+1,j+1,k+1,igx,igy,igz
     .                                  ,ig,jg,kg,x1p,x2p,x3p)
          x1 = 0.5*(x1+x1p)
          x2 = 0.5*(x2+x2p)
          x3 = 0.5*(x3+x3p)
        endif

        select case (coords)
        case ('car','lpk')
          jac = 1d0
        case ('cyl')
          jac = x1
        case ('hel')
          jac = x1
        case ('hl2')
          jac = 1d0
        case ('tor')
          major_r = gp(1)
          aa      = gp(2)
          bb      = gp(3)

          jac = aa*bb*x1*(major_r + aa*x1*cos(x2))
        case ('sin')
          pi = acos(-1d0)

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
      function covariantVector(comp,i,j,k,igx,igy,igz,vertex)
     .         result(vec)

c     -----------------------------------------------------------------
c     Calculates covariant vectors of curvilinear coordinate system
c     in Cartesian coordinates
c     -----------------------------------------------------------------

        implicit none

c     Input variables

        integer :: comp,i,j,k,igx,igy,igz
        logical,optional :: vertex

c     Local variables

        integer :: ig,jg,kg
        real(8) :: x1,x2,x3,aa,bb,major_r,mm,kk,eps,phi,jac,r
     .            ,x1p,x2p,x3p
        real(8) :: car(3),curv(3),vec(3),gp(ngparams)
        logical :: vrtx

c     Begin program

        if (PRESENT(vertex)) then
          vrtx = vertex
        else
          vrtx = .false.
        endif

        gp = grid_params%params

        call getCurvilinearCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
     .                                ,x1,x2,x3)

        if (vrtx) then
          call getCurvilinearCoordinates(i+1,j+1,k+1,igx,igy,igz
     .                                  ,ig,jg,kg,x1p,x2p,x3p)
          x1 = 0.5*(x1+x1p)
          x2 = 0.5*(x2+x2p)
          x3 = 0.5*(x3+x3p)
        endif

        select case (coords)
        case ('car','lpk')
          select case (comp)
            case (1)
              vec = (/ 1d0,0d0,0d0 /)
            case (2)
              vec = (/ 0d0,1d0,0d0 /)
            case (3)
              vec = (/ 0d0,0d0,1d0 /)
          end select
        case ('cyl')
          if (x1 == 0d0) then
            vec = 0d0
            return
          endif

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
          if (x1 == 0d0) then
            vec = 0d0
            return
          endif

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
          if (x1 == 0d0) then
            vec = 0d0
            return
          endif

          major_r = gp(1)
          aa      = gp(2)
          bb      = gp(3)

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
      function contravariantVector(comp,i,j,k,igx,igy,igz,vertex)
     .         result (vec)

c     -----------------------------------------------------------------
c     Calculates contravariant vectors of curvilinear coordinate system
c     in Cartesian coordinates
c     -----------------------------------------------------------------

        implicit none

c     Input variables

        integer :: comp,i,j,k,igx,igy,igz
        real(8) :: vec(3)
        logical,optional :: vertex

c     Local variables

        integer :: ig,jg,kg
        real(8) :: car(3),curv(3),jac
        real(8) :: x1,x2,x3,aa,bb,major_r,mm,kk,eps,phi,gp(ngparams)
     .            ,x1p,x2p,x3p
        logical :: vrtx

c     Begin program

        if (PRESENT(vertex)) then
          vrtx = vertex
        else
          vrtx = .false.
        endif

        gp = grid_params%params

        call getCurvilinearCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
     .                                ,x1,x2,x3)

        if (vrtx) then
          call getCurvilinearCoordinates(i+1,j+1,k+1,igx,igy,igz
     .                                  ,ig,jg,kg,x1p,x2p,x3p)
          x1 = 0.5*(x1+x1p)
          x2 = 0.5*(x2+x2p)
          x3 = 0.5*(x3+x3p)
        endif

        select case (coords)
        case ('car','lpk')
          select case (comp)
            case (1)
              vec = (/ 1d0,0d0,0d0 /)
            case (2)
              vec = (/ 0d0,1d0,0d0 /)
            case (3)
              vec = (/ 0d0,0d0,1d0 /)
          end select
        case ('cyl')
          if (x1 == 0d0) then
            vec = 0d0
            return
          endif

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
          if (x1 == 0d0) then
            vec = 0d0
            return
          endif

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
          if (x1 == 0d0) then
            vec = 0d0
            return
          endif

          major_r = gp(1)
          aa      = gp(2)
          bb      = gp(3)

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
        case ('car','lpk')
          vec = (/ 0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0 /)
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
        case ('car','lpk')
          vec = (/ 0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0 /)
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
      function g_sub(i,j,k,igx,igy,igz,vertex) result (tensor)

c     -----------------------------------------------------------------
c     Calculates contravariant metric tensor of curvilinear coordinate 
c     system
c     -----------------------------------------------------------------

        implicit none

c     Input variables

        integer :: i,j,k,igx,igy,igz
        real(8) :: tensor(3,3)
        logical,optional :: vertex

c     Local variables

        integer :: ig,jg,kg
        real(8) :: vec(9),jac,car(3),curv(3)
        real(8) :: x1,x2,x3,aa,bb,major_r,mm,kk,eps,gp(ngparams)
     .            ,x1p,x2p,x3p,i_x1
        logical :: vrtx

c     Begin program

        if (PRESENT(vertex)) then
          vrtx = vertex
        else
          vrtx = .false.
        endif

        gp = grid_params%params

        call getCurvilinearCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
     .                                ,x1,x2,x3)

        if (vrtx) then
          call getCurvilinearCoordinates(i+1,j+1,k+1,igx,igy,igz
     .                                  ,ig,jg,kg,x1p,x2p,x3p)
          x1 = 0.5*(x1+x1p)
          x2 = 0.5*(x2+x2p)
          x3 = 0.5*(x3+x3p)
        endif

        select case (coords)
        case ('car','lpk')
          vec = (/ 1d0, 0d0, 0d0
     .            ,0d0, 1d0, 0d0
     .            ,0d0, 0d0, 1d0 /)
        case ('cyl')
          if (x1 == 0d0) then
            vec = 0d0
          else
            vec = (/ 1d0/x1, 0d0  , 0d0
     .              ,0d0   , x1   , 0d0
     .              ,0d0   , 0d0  , 1d0/x1 /)
          endif
        case ('hel')
          mm = gp(1)
          kk = gp(2)
          aa = kk/mm
          if (x1 == 0d0) then
            vec = 0d0
          else
            vec = (/ 1./x1, 0d0  ,      0d0
     .              ,0d0  ,  x1  ,    -aa*x1
     .              ,0d0  ,-aa*x1, 1./x1 + aa**2*x1 /)
          endif
        case ('hl2')
          curv = (/ abs(2*x1),x2,x3 /)

          mm = gp(1)
          kk = gp(2)
          aa = kk/mm
          vec = (/ 1./curv(1),0d0, 0d0
     .            ,0d0, curv(1) ,-aa*curv(1)
     .            ,0d0,-aa*curv(1), 1 + aa**2*curv(1) /)
        case ('tor')
          if (x1 == 0d0) then
            vec = 0d0
          else
            major_r = gp(1)
            aa      = gp(2)
            bb      = gp(3)

            vec(1) = (aa**2 + bb**2 + (aa**2 - bb**2)*Cos(2*x2))/
     -               (2.*aa*bb*x1*(major_r + aa*x1*Cos(x2)))

            vec(2) = ((-aa**2 + bb**2)*Cos(x2)*Sin(x2))
     .               /(aa*bb*(major_r + aa*x1*Cos(x2)))

            vec(3) = 0d0

            vec(4) = vec(2)

            vec(5) = (x1*(aa**2 + bb**2 + (-aa**2 + bb**2)*Cos(2*x2)))/
     -               (2.*aa*bb*(major_r + aa*x1*Cos(x2)))

            vec(6) = 0d0

            vec(7) = vec(3)
            vec(8) = vec(6)
            vec(9) = (major_r + aa*x1*Cos(x2))/(aa*bb*x1)
          endif
        case ('sin')

          pi = acos(-1d0)

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

        tensor = reshape(vec,(/3,3/))

      end function g_sub

c     g_sup
c     #################################################################
      function g_sup(i,j,k,igx,igy,igz,vertex) result (tensor)

c     -----------------------------------------------------------------
c     Calculates covariant metric tensor of curvilinear coordinate 
c     system
c     -----------------------------------------------------------------

        implicit none

c     Input variables

        integer :: i,j,k,igx,igy,igz
        real(8) :: tensor(3,3)
        logical,optional :: vertex

c     Local variables

        integer :: ig,jg,kg
        real(8) :: vec(9),jac,car(3),curv(3)
        real(8) :: x1,x2,x3,aa,bb,major_r,mm,kk,eps,gp(ngparams)
     .            ,x1p,x2p,x3p
        logical :: vrtx

c     Begin program

        if (PRESENT(vertex)) then
          vrtx = vertex
        else
          vrtx = .false.
        endif

        gp = grid_params%params

        call getCurvilinearCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
     .                                ,x1,x2,x3)

        if (vrtx) then
          call getCurvilinearCoordinates(i+1,j+1,k+1,igx,igy,igz
     .                                  ,ig,jg,kg,x1p,x2p,x3p)
          x1 = 0.5*(x1+x1p)
          x2 = 0.5*(x2+x2p)
          x3 = 0.5*(x3+x3p)
        endif

        select case (coords)
        case ('car','lpk')
          vec = (/ 1d0, 0d0, 0d0
     .            ,0d0, 1d0, 0d0
     .            ,0d0, 0d0, 1d0 /)
        case ('cyl')
          if (x1 == 0d0) then
            vec = 0d0
          else
            vec = (/ x1 , 0d0   , 0d0
     .              ,0d0, 1d0/x1, 0d0
     .              ,0d0, 0d0   , x1 /)
          endif
        case ('hel')
          mm = gp(1)
          kk = gp(2)
          aa = kk/mm
          if (x1 == 0d0) then
            vec = 0d0
          else
            vec = (/  x1  , 0d0              , 0d0
     .              , 0d0 , 1./x1 + aa**2*x1 , aa*x1
     .              , 0d0 , aa*x1            , x1 /)
          endif
        case ('hl2')
          mm = gp(1)
          kk = gp(2)
          aa = kk/mm

          x1 = abs(2*x1)

          if (x1 == 0d0) then
            vec = 0d0
          else
            vec = (/ x1 , 0d0          , 0d0
     .              ,0d0, 1./x1 + aa**2, aa
     .              ,0d0, aa           , 1d0 /)
          endif
        case ('tor')
          if (x1 == 0d0) then
            vec = 0d0
          else
            major_r = gp(1)
            aa      = gp(2)
            bb      = gp(3)

            vec(1) = -(x1*(major_r + aa*x1*Cos(x2))*
     .               (-aa**2 - bb**2 + (aa**2 - bb**2)*Cos(2*x2)))
     .               /(2.*aa*bb)

            vec(2) = ((aa**2 - bb**2)*Cos(x2)
     .               *(major_r + aa*x1*Cos(x2))*Sin(x2))
     .               /(aa*bb)

            vec(3) = 0d0

            vec(4) = vec(2)

            vec(5) = ((major_r + aa*x1*Cos(x2))
     .              *(aa**2 + bb**2 + (aa**2 - bb**2)*Cos(2*x2)))
     .              /(2.*aa*bb*x1)

            vec(6) = 0d0

            vec(7) = vec(3)
            vec(8) = vec(6)
            vec(9) = (aa*bb*x1)/(major_r + aa*x1*Cos(x2))
          endif
        case ('sin')

          pi = acos(-1d0)

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

        tensor = reshape(vec,(/3,3/))

      end function g_sup

      end module grid_anal_map
