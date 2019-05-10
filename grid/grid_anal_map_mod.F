
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
        case('car','cyl','hel','hl2','tor','sin','lpk','s3d','hsn'
     .      ,'sha','sh2','stl','tsq')

          anal_map = .true.

        case default

          anal_map = .false.

        end select

      end function checkAnalMapDatabase

c     map
c     #################################################################
      function map(g_def,i,j,k,igx,igy,igz,ig,jg,kg,vertex) result(car)

c     -----------------------------------------------------------------
c     Give Cartesian coordinates corresponding to node (i,j,k) at grid
c     level (igx,igy,igz).
c     -----------------------------------------------------------------

        implicit none

c     Input variables

        type(grid_mg_def),pointer :: g_def

        integer :: i,j,k,igx,igy,igz
        real(8) :: car(3)
        integer,optional,intent(OUT) :: ig,jg,kg
        logical,optional :: vertex

c     Local variables

        integer :: ig1,jg1,kg1

        real(8) :: x1,y1,z1

c     External routines

c     Begin program

        call getCurvilinearCoordinates(g_def,i,j,k,igx,igy,igz
     .                                ,ig1,jg1,kg1,x1,y1,z1
     .                                ,vertex=vertex)

        car = x_xi(g_def,x1,y1,z1)

        if (PRESENT(ig)) ig = ig1
        if (PRESENT(jg)) jg = jg1
        if (PRESENT(kg)) kg = kg1

      end function map

c     x_xi
c     #################################################################
      function x_xi(gl_def,x1,x2,x3) result(car)

c     -----------------------------------------------------------------
c     Gives Cartesian coordinates from curvilinear coordinates
c     -----------------------------------------------------------------

        implicit none

c     Input variables

        real(8) :: x1,x2,x3,car(3)

        type(grid_mg_def),pointer :: gl_def

c     Local variables

        integer :: inewt,ic
        real(8) :: xx,yy,zz,jac_mat(3,3),rhs(3),dx(3),rr,rr0,dd
     .            ,aa,bb,major_r,mm,kk,eps,phi,gp(ngparams)
     .            ,xlog(3),period1,period2
        real(8) :: xmax,ymax,zmax,xmin,ymin,zmin
! Daniele, 11/07/2017. Local variables needed for shaped toroidal coordinates
        integer :: i
        real(8) :: ai(10),bi(10)
! Daniele, 11/08/2017. Local variables needed for stellarator (helias) coordinates
        real(8) :: nn,uu,vv,d10,d0,d20,d1m1,d2m1,d2m2,rc,zc

c     Begin program

        xmin = gl_def%gxmin
        xmax = gl_def%gxmax

        ymin = gl_def%gymin
        ymax = gl_def%gymax

        zmin = gl_def%gzmin
        zmax = gl_def%gzmax

        if (coords /= 'lpk') then
          xlog(1) = scale_coord(gl_def,1,x1)
          xlog(2) = scale_coord(gl_def,2,x2)
          xlog(3) = scale_coord(gl_def,3,x3)
        else
          xlog = (/ x1, x2, x3 /)
        endif

        gp = gl_def%params

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

          phi = xlog(2)-aa*xlog(3)
          car(1) = xlog(1)*cos(phi)
          car(2) = xlog(1)*sin(phi)
          car(3) = xlog(3)
        case ('hl2')
          mm = gp(1)
          kk = gp(2)
          aa = kk/mm

          phi = xlog(2)-aa*xlog(3)
          car(1) = sign(sqrt(2*abs(xlog(1))),xlog(1))*cos(phi)
          car(2) = sign(sqrt(2*abs(xlog(1))),xlog(1))*sin(phi)
          car(3) = xlog(3)
        case ('tor')
          major_r = gp(1)  !Major radius
          aa      = gp(2)  !Horizontal minor radius
          bb      = gp(3)  !Vertical   minor radius

          car(1) = (major_r + aa*xlog(1)*cos(xlog(2)))*cos(-xlog(3))
          car(2) = (major_r + aa*xlog(1)*cos(xlog(2)))*sin(-xlog(3))
          car(3) = bb*xlog(1)*sin(xlog(2))
        case ('sha')
! Daniele, 11/07/2017
! Shaped toroidal coordinates system
! x = (R0 + r*sum_i(ai*cos(i*theta))*cos(-phi)
! y = (R0 + r*sum_i(ai*cos(i*theta))*sin(-phi)
! z = r*sum_i(bi*sin(i*theta)

          major_r = gp(1)  !Major radius

          ai(1) = gp(2)  !Horizontal minor radius
          bi(1) = gp(3)  !Vertical   minor radius

          do i=2,10
             ai(i) = gp(7+i) ! Higher order cosine coefficients are stored in gp(9) to gp(17)
             bi(i) = gp(16+i) ! Higher order sine coefficients are stored in gp(18) to gp(26)
          enddo

          aa = 0d0 ; bb = 0d0
          do i=1,10
             aa = aa + ai(i)*cos(i*xlog(2))
             bb = bb + bi(i)*sin(i*xlog(2))
          enddo

          car(1) = (major_r + aa*xlog(1))*cos(-xlog(3))
          car(2) = (major_r + aa*xlog(1))*sin(-xlog(3))
          car(3) = bb*xlog(1)

        case ('sh2')
! Daniele, 11/07/2017
! Shaped toroidal coordinates system
! x = (R0 + r*a1*cos(theta)+r^2*sum_i(ai*cos(i*theta))*cos(-phi)
! y = (R0 + r*a1*cos(theta)+r^2*sum_i(ai*cos(i*theta))*sin(-phi)
!! x = (R0 + r*sum_i(r^(i-1)*ai*cos(i*theta))*cos(-phi)
!! y = (R0 + r*sum_i(r^(i-1)*ai*cos(i*theta))*sin(-phi)
! z = r*sum_i(bi*sin(i*theta) <<< INCORRECT???

          major_r = gp(1)  !Major radius

          ai(1) = gp(2)  !Horizontal minor radius
          bi(1) = gp(3)  !Vertical   minor radius

          do i=2,10
             ai(i) = gp(7+i) ! Higher order cosine coefficients are stored in gp(9) to gp(17)
             bi(i) = gp(16+i) ! Higher order sine coefficients are stored in gp(18) to gp(26)
          enddo

          aa = ai(1)*cos(1*xlog(2))
          bb = bi(1)*sin(1*xlog(2))
          do i=2,10
             aa = aa + xlog(1)*ai(i)*cos(i*xlog(2))
             bb = bb + xlog(1)*bi(i)*sin(i*xlog(2))
          enddo
c$$$          aa = 0d0 ; bb = 0d0
c$$$          do i=1,10
c$$$             aa = aa + xlog(1)**(i-1)*ai(i)*cos(i*xlog(2))
c$$$             bb = bb + xlog(1)**(i-1)*bi(i)*sin(i*xlog(2))
c$$$          enddo

          car(1) = (major_r + aa*xlog(1))*cos(-xlog(3))
          car(2) = (major_r + aa*xlog(1))*sin(-xlog(3))
          car(3) = bb*xlog(1)

        case ('tsq')  !Square torus
          major_r = gp(1)  !Major radius
          aa      = gp(2)  !Horizontal dimension
          bb      = gp(3)  !Vertical   dimension

          car(1) = (major_r + aa*(xlog(1)-0.5))*cos(-xlog(3))
          car(2) = (major_r + aa*(xlog(1)-0.5))*sin(-xlog(3))
          car(3) =            bb*(xlog(2)-0.5)

        case ('sin')
          eps = gp(1)

          car(1) = xlog(1) + eps*sin(2*pi*xlog(1)/xmax)
     .                          *sin(2*pi*xlog(2)/ymax)
          car(2) = xlog(2) + eps*sin(2*pi*xlog(1)/xmax)
     .                          *sin(2*pi*xlog(2)/ymax)
          car(3) = xlog(3)
        case ('hsn')
          eps = gp(1)
          period1 = 2*pi
          if (bcond(1) /= PER) period1 = gp(2)
          period2 = 2*pi
          if (bcond(3) /= PER) period2 = gp(3)

          car(1) = xlog(1) + eps*sin(period1*xlog(1)/xmax)
     .                          *sin(period2*xlog(2)/ymax)
          car(2) = xlog(2) + eps*sin(period1*xlog(1)/xmax)
     .                          *sin(period2*xlog(2)/ymax)
          car(3) = xlog(3)
        case ('s3d')
          eps = gp(1)

          car(1) = xlog(1) + eps*sin(2*pi*xlog(1)/xmax)
     .                          *sin(2*pi*xlog(2)/ymax)
     .                          *sin(2*pi*xlog(3)/zmax)
          car(2) = xlog(2) + eps*sin(2*pi*xlog(1)/xmax)
     .                          *sin(2*pi*xlog(2)/ymax)
     .                          *sin(2*pi*xlog(3)/zmax)
          car(3) = xlog(3) + eps*sin(2*pi*xlog(1)/xmax)
     .                          *sin(2*pi*xlog(2)/ymax)
     .                          *sin(2*pi*xlog(3)/zmax)

        case ('stl')
! Daniele, 11/08/2017
! Helias stellarator coordinate system
! see Nuhrenberg-Zille1986PL

          major_r = gp(1)  !Major radius

          aa = gp(2)       !Horizontal minor radius
          bb = gp(3)       !Vertical   minor radius

          nn   = gp(9)
          d10  = gp(10)
          d0   = gp(11)
          d20  = gp(12)
          d1m1 = gp(13)
          d2m1 = gp(14)
          d2m2 = gp(15)

          rc = gp(16)
          zc = gp(17)

          vv = nn*xlog(3)
          uu = xlog(2)

          rr=major_r
     .       +aa*(1d0+rc*(1d0-xlog(1)))*cos(vv)
     .       +xlog(1)*((1d0-d10-d0*cos(vv))*cos(uu)
     .       +sqrt(abs(xlog(1)))
     .           *(d20*cos(2*uu)
     .            -d1m1*cos(uu-vv)
     .            +d2m1*cos(2*uu-vv)
     .            +d2m2*cos(2*uu-2*vv)))
          zz=bb*(1d0+zc*(1d0-xlog(1)))*sin(vv)
     .       +xlog(1)*((1d0+d10-d0*cos(vv))*sin(uu)
     .       +sqrt(abs(xlog(1)))
     .           *(d20*sin(2*uu)
     .            +d1m1*sin(uu-vv)
     .            +d2m1*sin(2*uu-vv)
     .            -d2m2*sin(2*uu-2*vv)))

          car(1) = rr*cos(-xlog(3))
          car(2) = rr*sin(-xlog(3))
          car(3) = zz

        case default
          write (*,*) 'Grid not implemented in x_xi'
          write (*,*) 'Aborting...'
          stop
        end select

      end function x_xi

c     jacobian
c     #################################################################
      function jacobian(gl_def,i,j,k,igx,igy,igz,vertex) result(jac)

c     -----------------------------------------------------------------
c     Calculates Jacobian of curvilinear coordinate system
c     -----------------------------------------------------------------

        implicit none

c     Input variables

        integer :: i,j,k,igx,igy,igz
        real(8) :: jac
        logical,optional :: vertex

        type(grid_mg_def),pointer :: gl_def

c     Local variables

        integer :: ig,jg,kg
        real(8) :: x1,x2,x3,car(3),curv(3),aa,bb,major_r,mm,kk,eps
     .            ,gp(ngparams),x1p,x2p,x3p
        logical :: vrtx

        real(8) :: xmax,ymax,zmax,xmin,ymin,zmin
! Daniele, 11/08/2017. Local variables needed for shaped toroidal coordinates
        integer :: ii
        real(8) :: ai(10),bi(10)
        real(8) :: ap,bp

c     Begin program

        xmin = gl_def%gxmin
        xmax = gl_def%gxmax

        ymin = gl_def%gymin
        ymax = gl_def%gymax

        zmin = gl_def%gzmin
        zmax = gl_def%gzmax

        if (PRESENT(vertex)) then
          vrtx = vertex
        else
          vrtx = .false.
        endif

        gp = gl_def%params

        call getCurvilinearCoordinates(gl_def,i,j,k,igx,igy,igz,ig,jg,kg
     .                                ,x1,x2,x3)

        if (vrtx) then
          call getCurvilinearCoordinates(gl_def,i+1,j+1,k+1,igx,igy,igz
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
        case ('sha')
! Daniele, 11/08/2017
! Shaped toroidal coordinates system
! x = (R0 + r*sum_i(ai*cos(i*theta))*cos(-phi)
! y = (R0 + r*sum_i(ai*cos(i*theta))*sin(-phi)
! z = r*sum_i(bi*sin(i*theta)

          major_r = gp(1)  !Major radius

          ai(1) = gp(2)  !Horizontal minor radius
          bi(1) = gp(3)  !Vertical   minor radius

          do ii=2,10
             ai(ii) = gp(7+ii) ! Higher order cosine coefficients are stored in gp(9) to gp(17)
             bi(ii) = gp(16+ii) ! Higher order sine coefficients are stored in gp(18) to gp(26)
          enddo

          aa = 0d0
          bb = 0d0
          ap = 0d0
          bp = 0d0
          do ii=1,10
             aa = aa + ai(ii)*cos(ii*x2) ! A(theta)=sum_i(ai*cos(i*theta))
             bb = bb + bi(ii)*sin(ii*x2) ! B(theta)=sum_i(bi*sin(i*theta))
             bp = bp + bi(ii)*ii*cos(ii*x2) ! B'(theta)
             ap = ap - ai(ii)*ii*sin(ii*x2) ! A'(theta)
          enddo

          jac = (aa*bp-bb*ap)*x1*(major_r + aa*x1)

        case ('sin')
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
      function covariantVector(gl_def,comp,i,j,k,igx,igy,igz,vertex)
     .         result(vec)

c     -----------------------------------------------------------------
c     Calculates covariant vectors of curvilinear coordinate system
c     in Cartesian coordinates
c     -----------------------------------------------------------------

        implicit none

c     Input variables

        integer :: comp,i,j,k,igx,igy,igz
        logical,optional :: vertex

        type(grid_mg_def),pointer :: gl_def

c     Local variables

        integer :: ig,jg,kg
        real(8) :: x1,x2,x3,aa,bb,major_r,mm,kk,eps,phi,jac,r
     .            ,x1p,x2p,x3p
        real(8) :: car(3),curv(3),vec(3),gp(ngparams)
        logical :: vrtx

        real(8) :: xmax,ymax,zmax,xmin,ymin,zmin
! Daniele, 11/08/2017. Local variables needed for shaped toroidal coordinates
        integer :: ii
        real(8) :: ai(10),bi(10)
        real(8) :: ap,bp

c     Begin program

        xmin = gl_def%gxmin
        xmax = gl_def%gxmax

        ymin = gl_def%gymin
        ymax = gl_def%gymax

        zmin = gl_def%gzmin
        zmax = gl_def%gzmax

        if (PRESENT(vertex)) then
          vrtx = vertex
        else
          vrtx = .false.
        endif

        gp = gl_def%params

        call getCurvilinearCoordinates(gl_def,i,j,k,igx,igy,igz,ig,jg,kg
     .                                ,x1,x2,x3)

        if (vrtx) then
          call getCurvilinearCoordinates(gl_def,i+1,j+1,k+1,igx,igy,igz
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
        case ('sha')
! Daniele, 11/08/2017
! Shaped toroidal coordinates system
! x = (R0 + r*sum_i(ai*cos(i*theta))*cos(-phi)
! y = (R0 + r*sum_i(ai*cos(i*theta))*sin(-phi)
! z = r*sum_i(bi*sin(i*theta)

          if (x1 == 0d0) then
            vec = 0d0
            return
          endif

          major_r = gp(1)  !Major radius

          ai(1) = gp(2)  !Horizontal minor radius
          bi(1) = gp(3)  !Vertical   minor radius

          do ii=2,10
             ai(ii) = gp(7+ii) ! Higher order cosine coefficients are stored in gp(9) to gp(17)
             bi(ii) = gp(16+ii) ! Higher order sine coefficients are stored in gp(18) to gp(26)
          enddo

          curv = (/ x1,x2,x3 /)

          aa = 0d0
          bb = 0d0
          ap = 0d0
          bp = 0d0
          do ii=1,10
             aa = aa + ai(ii)*cos(ii*curv(2)) ! A(theta)=sum_i(ai*cos(i*theta))
             bb = bb + bi(ii)*sin(ii*curv(2)) ! B(theta)=sum_i(bi*sin(i*theta))
             bp = bp + bi(ii)*ii*cos(ii*curv(2)) ! B'(theta)
             ap = ap - ai(ii)*ii*sin(ii*curv(2)) ! A'(theta)
          enddo

          select case (comp)
            case (1)
              vec = (/  bp*cos(curv(3))
     .                ,-bp*sin(curv(3))
     .                ,-ap             /)/(aa*bp-bb*ap)
            case (2)
              vec = (/ -bb*cos(curv(3))
     .                , bb*sin(curv(3))
     .                , aa             /)/curv(1)/(aa*bp-bb*ap)
            case (3)
              vec = (/ -sin(curv(3)),-cos(curv(3)),0d0 /)
     .              /(major_r + aa*curv(1))
          end select

        case ('sin')

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

c     cov_vector_anal
c     #################################################################
      function cov_vector_anal(gl_def,curv) result(vec)

c     -----------------------------------------------------------------
c     Calculates covariant vectors of curvilinear coordinate system
c     in Cartesian coordinates
c     -----------------------------------------------------------------

        implicit none

c     Input variables

        real(8) :: curv(3),vec(3,3)

        type(grid_mg_def),pointer :: gl_def

c     Local variables

        real(8) :: aa,bb,major_r,mm,kk,eps,phi,jac,r
     .            ,x1,x2,x3
        real(8) :: gp(ngparams)

        real(8) :: xmax,ymax,zmax,xmin,ymin,zmin
! Daniele, 11/08/2017. Local variables needed for shaped toroidal coordinates
        integer :: ii
        real(8) :: ai(10),bi(10)
        real(8) :: ap,bp

c     Begin program

        xmin = gl_def%gxmin
        xmax = gl_def%gxmax

        ymin = gl_def%gymin
        ymax = gl_def%gymax

        zmin = gl_def%gzmin
        zmax = gl_def%gzmax

        gp = gl_def%params

        select case (coords)
        case ('car','lpk')
          vec(1,:) = (/ 1d0,0d0,0d0 /)
          vec(2,:) = (/ 0d0,1d0,0d0 /)
          vec(3,:) = (/ 0d0,0d0,1d0 /)
        case ('cyl')
          if (curv(1) == 0d0) then
            vec = 0d0
            return
          endif

          vec(1,:) = (/ cos(curv(2)),sin(curv(2)),0d0 /)
          vec(2,:) = (/-sin(curv(2)),cos(curv(2)),0d0 /)/curv(1)
          vec(3,:) = (/ 0d0,0d0,1d0 /)
        case ('hel')
          if (curv(1) == 0d0) then
            vec = 0d0
            return
          endif

          mm = gp(1)
          kk = gp(2)
          aa = kk/mm
          phi = (curv(2)-aa*curv(3))

          vec(1,:) = (/ cos(phi),sin(phi),0d0 /)
          vec(2,:) = (/-sin(phi)/curv(1),cos(phi)/curv(1),aa /)
          vec(3,:) = (/ 0d0,0d0,1d0 /)
cc        case ('hl2')
cc          curv = (/ sign(sqrt(2*abs(x1)),x1),x2,x3 /)
cc
cc          mm = gp(1)
cc          kk = gp(2)
cc          aa = kk/mm
cc          phi = (curv(2)-aa*curv(3))
cc
cc          vec(1,:) = (/ cos(phi),sin(phi),0d0 /)*curv(1)
cc          vec(2,:) = (/-sin(phi)/curv(1),cos(phi)/curv(1),aa /)
cc          vec(3,:) = (/ 0d0,0d0,1d0 /)
        case ('tor')
          if (curv(1) == 0d0) then
            vec = 0d0
            return
          endif

          major_r = gp(1)
          aa      = gp(2)
          bb      = gp(3)

          vec(1,:) = (/  cos(curv(2))*cos(curv(3))/aa
     .                 ,-cos(curv(2))*sin(curv(3))/aa
     .                 , sin(curv(2))/bb             /)
          vec(2,:) = (/-sin(curv(2))*cos(curv(3))/aa
     .                 ,sin(curv(2))*sin(curv(3))/aa
     .                 ,cos(curv(2))/bb             /)/curv(1)
          vec(3,:) = (/ -sin(curv(3)),-cos(curv(3)),0d0 /)
     .                /(major_r + aa*curv(1)*cos(curv(2)))
        case ('sha')
! Daniele, 11/08/2017
! Shaped toroidal coordinates system
! x = (R0 + r*sum_i(ai*cos(i*theta))*cos(-phi)
! y = (R0 + r*sum_i(ai*cos(i*theta))*sin(-phi)
! z = r*sum_i(bi*sin(i*theta)

          if (curv(1) == 0d0) then
            vec = 0d0
            return
          endif

          major_r = gp(1)  !Major radius

          ai(1) = gp(2)  !Horizontal minor radius
          bi(1) = gp(3)  !Vertical   minor radius

          do ii=2,10
             ai(ii) = gp(7+ii) ! Higher order cosine coefficients are stored in gp(9) to gp(17)
             bi(ii) = gp(16+ii) ! Higher order sine coefficients are stored in gp(18) to gp(26)
          enddo

          aa = 0d0
          bb = 0d0
          ap = 0d0
          bp = 0d0
          do ii=1,10
             aa = aa + ai(ii)*cos(ii*curv(2)) ! A(theta)=sum_i(ai*cos(i*theta))
             bb = bb + bi(ii)*sin(ii*curv(2)) ! B(theta)=sum_i(bi*sin(i*theta))
             bp = bp + bi(ii)*ii*cos(ii*curv(2)) ! B'(theta)
             ap = ap - ai(ii)*ii*sin(ii*curv(2)) ! A'(theta)
          enddo

          vec(1,:) = (/  bp*cos(curv(3))
     .                 ,-bp*sin(curv(3))
     .                 ,-ap             /)/(aa*bp-bb*ap)
          vec(2,:) = (/ -bb*cos(curv(3))
     .                 , bb*sin(curv(3))
     .                 , aa             /)/curv(1)/(aa*bp-bb*ap)
          vec(3,:) = (/ -sin(curv(3)),-cos(curv(3)),0d0 /)
     .                /(major_r + aa*curv(1))

        case ('sin')

          eps = gp(1)

          x1 = curv(1)
          x2 = curv(2)
          x3 = curv(3)
          vec(1,:) = (/(xmax*(ymax + 2*eps*Pi*Cos((2*Pi*x2)/ymax)*
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
          vec(2,:) = (/(-2*eps*Pi*ymax*Cos((2*Pi*x1)/xmax)
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
          vec(3,:) = (/ 0d0,0d0,1d0 /)
        case default
          write (*,*) 'Grid not implemented in covariantVector'
          write (*,*) 'Aborting...'
          stop
        end select

      end function cov_vector_anal

c     contravariantVector
c     #################################################################
      function contravariantVector(gl_def,comp,i,j,k,igx,igy,igz,vertex)
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

        type(grid_mg_def),pointer :: gl_def

c     Local variables

        integer :: ig,jg,kg
        real(8) :: car(3),curv(3),jac
        real(8) :: x1,x2,x3,aa,bb,major_r,mm,kk,eps,phi,gp(ngparams)
     .            ,x1p,x2p,x3p
        logical :: vrtx

        real(8) :: xmax,ymax,zmax,xmin,ymin,zmin
! Daniele, 11/08/2017. Local variables needed for shaped toroidal coordinates
        integer :: ii
        real(8) :: ai(10),bi(10)
        real(8) :: ap,bp

c     Begin program

        xmin = gl_def%gxmin
        xmax = gl_def%gxmax

        ymin = gl_def%gymin
        ymax = gl_def%gymax

        zmin = gl_def%gzmin
        zmax = gl_def%gzmax

        if (PRESENT(vertex)) then
          vrtx = vertex
        else
          vrtx = .false.
        endif

        gp = gl_def%params

        call getCurvilinearCoordinates(gl_def,i,j,k,igx,igy,igz,ig,jg,kg
     .                                ,x1,x2,x3)

        if (vrtx) then
          call getCurvilinearCoordinates(gl_def,i+1,j+1,k+1,igx,igy,igz
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
     .              /curv(1)/(major_r + aa*curv(1)*cos(curv(2)))
            case (2)
              vec = (/-sin(curv(2))*cos(curv(3))/bb
     .                ,sin(curv(2))*sin(curv(3))/bb
     .                ,cos(curv(2))/aa              /)
     .              /(major_r + aa*curv(1)*cos(curv(2)))
            case (3)
              vec = (/ -sin(curv(3)),-cos(curv(3)),0d0 /)/curv(1)/aa/bb
          end select

        case ('sha')
! Daniele, 11/08/2017
! Shaped toroidal coordinates system
! x = (R0 + r*sum_i(ai*cos(i*theta))*cos(-phi)
! y = (R0 + r*sum_i(ai*cos(i*theta))*sin(-phi)
! z = r*sum_i(bi*sin(i*theta)

          if (x1 == 0d0) then
            vec = 0d0
            return
          endif

          major_r = gp(1)  !Major radius

          ai(1) = gp(2)  !Horizontal minor radius
          bi(1) = gp(3)  !Vertical   minor radius

          do ii=2,10
             ai(ii) = gp(7+ii) ! Higher order cosine coefficients are stored in gp(9) to gp(17)
             bi(ii) = gp(16+ii) ! Higher order sine coefficients are stored in gp(18) to gp(26)
          enddo

          curv = (/ x1,x2,x3 /)

          aa = 0d0
          bb = 0d0
          ap = 0d0
          bp = 0d0
          do ii=1,10
             aa = aa + ai(ii)*cos(ii*curv(2)) ! A(theta)=sum_i(ai*cos(i*theta))
             bb = bb + bi(ii)*sin(ii*curv(2)) ! B(theta)=sum_i(bi*sin(i*theta))
             bp = bp + bi(ii)*ii*cos(ii*curv(2)) ! B'(theta)
             ap = ap - ai(ii)*ii*sin(ii*curv(2)) ! A'(theta)
          enddo

          select case (comp)
            case (1)
              vec = (/  aa*cos(curv(3))
     .                ,-aa*sin(curv(3))
     .                , bb             /)
     .              /curv(1)/(major_r + aa*curv(1))/(aa*bp-bb*ap)
            case (2)
              vec = (/  ap*cos(curv(3))
     .                ,-ap*sin(curv(3))
     .                , bp             /)
     .              /(major_r + aa*curv(1))/(aa*bp-bb*ap)
            case (3)
              vec = (/ -sin(curv(3)),-cos(curv(3)),0d0 /)
     .              /curv(1)/(aa*bp-bb*ap)
          end select

        case ('sin')

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

c     cnv_vector_anal
c     #################################################################
      function cnv_vector_anal(gl_def,curv) result(vec)

c     -----------------------------------------------------------------
c     Calculates covariant vectors of curvilinear coordinate system
c     in Cartesian coordinates
c     -----------------------------------------------------------------

        implicit none

c     Input variables

        real(8) :: curv(3),vec(3,3)

        type(grid_mg_def),pointer :: gl_def

c     Local variables

        real(8) :: aa,bb,major_r,mm,kk,eps,phi,jac,r
     .            ,x1,x2,x3
        real(8) :: gp(ngparams)

        real(8) :: xmax,ymax,zmax,xmin,ymin,zmin
! Daniele, 11/08/2017. Local variables needed for shaped toroidal coordinates
        integer :: ii
        real(8) :: ai(10),bi(10)
        real(8) :: ap,bp

c     Begin program

        xmin = gl_def%gxmin
        xmax = gl_def%gxmax

        ymin = gl_def%gymin
        ymax = gl_def%gymax

        zmin = gl_def%gzmin
        zmax = gl_def%gzmax

        gp = gl_def%params

        select case (coords)
        case ('car','lpk')
          vec(1,:) = (/ 1d0,0d0,0d0 /)
          vec(2,:) = (/ 0d0,1d0,0d0 /)
          vec(3,:) = (/ 0d0,0d0,1d0 /)
        case ('cyl')
          if (curv(1) == 0d0) then
            vec = 0d0
            return
          endif

          vec(1,:) = (/ cos(curv(2)),sin(curv(2)),0d0 /)/curv(1)
          vec(2,:) = (/-sin(curv(2)),cos(curv(2)),0d0 /)
          vec(3,:) = (/ 0d0,0d0,1d0 /)/curv(1)
        case ('hel')
          if (curv(1) == 0d0) then
            vec = 0d0
            return
          endif

          mm = gp(1)
          kk = gp(2)
          aa = kk/mm
          phi = (curv(2)-aa*curv(3))

          vec(1,:) = (/ cos(phi),sin(phi),0d0 /)/curv(1)
          vec(2,:) = (/-sin(phi),cos(phi),0d0 /)
          vec(3,:) = (/ aa*sin(phi),-aa*cos(phi),1d0/curv(1) /)
cc        case ('hl2')
cc          curv = (/ sign(sqrt(2*abs(x1)),x1),x2,x3 /)
cc
cc          mm = gp(1)
cc          kk = gp(2)
cc          aa = kk/mm
cc          phi = (curv(2)-aa*curv(3))
cc
cc          vec(1,:) = (/ cos(phi),sin(phi),0d0 /)/curv(1)
cc          vec(2,:) = (/-sin(phi),cos(phi),0d0 /)*curv(1)
cc          vec(3,:) = (/ aa*sin(phi)*curv(1),-aa*cos(phi)*curv(1),1d0 /)
        case ('tor')
          if (curv(1) == 0d0) then
            vec = 0d0
            return
          endif

          major_r = gp(1)
          aa      = gp(2)
          bb      = gp(3)

          vec(1,:) = (/  cos(curv(2))*cos(curv(3))/bb
     .                ,-cos(curv(2))*sin(curv(3))/bb
     .                , sin(curv(2))/aa             /)
     .              /curv(1)/(major_r + aa*curv(1)*cos(curv(2)))
          vec(2,:) = (/-sin(curv(2))*cos(curv(3))/bb
     .                ,sin(curv(2))*sin(curv(3))/bb
     .                ,cos(curv(2))/aa              /)
     .              /(major_r + aa*curv(1)*cos(curv(2)))
          vec(3,:) = (/ -sin(curv(3)),-cos(curv(3)),0d0 /)/curv(1)/aa/bb

        case ('sha')
! Daniele, 11/08/2017
! Shaped toroidal coordinates system
! x = (R0 + r*sum_i(ai*cos(i*theta))*cos(-phi)
! y = (R0 + r*sum_i(ai*cos(i*theta))*sin(-phi)
! z = r*sum_i(bi*sin(i*theta)

          if (curv(1) == 0d0) then
            vec = 0d0
            return
          endif

          major_r = gp(1)  !Major radius

          ai(1) = gp(2)  !Horizontal minor radius
          bi(1) = gp(3)  !Vertical   minor radius

          do ii=2,10
             ai(ii) = gp(7+ii) ! Higher order cosine coefficients are stored in gp(9) to gp(17)
             bi(ii) = gp(16+ii) ! Higher order sine coefficients are stored in gp(18) to gp(26)
          enddo

          aa = 0d0
          bb = 0d0
          ap = 0d0
          bp = 0d0
          do ii=1,10
             aa = aa + ai(ii)*cos(ii*curv(2)) ! A(theta)=sum_i(ai*cos(i*theta))
             bb = bb + bi(ii)*sin(ii*curv(2)) ! B(theta)=sum_i(bi*sin(i*theta))
             bp = bp + bi(ii)*ii*cos(ii*curv(2)) ! B'(theta)
             ap = ap - ai(ii)*ii*sin(ii*curv(2)) ! A'(theta)
          enddo

          vec(1,:) = (/ aa*cos(curv(3))
     .                ,-aa*sin(curv(3))
     .                , bb             /)
     .              /curv(1)/(major_r + aa*curv(1))/(aa*bp-bb*ap)
          vec(2,:) = (/ ap*cos(curv(3))
     .                ,-ap*sin(curv(3))
     .                , bp             /)
     .              /(major_r + aa*curv(1))/(aa*bp-bb*ap)
          vec(3,:) = (/ -sin(curv(3)),-cos(curv(3)),0d0 /)
     .              /curv(1)/(aa*bp-bb*ap)

        case ('sin')

          eps = gp(1)

          x1 = curv(1)
          x2 = curv(2)
          x3 = curv(3)

          vec(1,:) =(/ (ymax*(xmax + 2*eps*Pi*Cos((2*Pi*x1)/xmax)*
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
          vec(2,:) =(/ (2*eps*Pi*xmax*Cos((2*Pi*x2)/ymax)
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
          vec(3,:) =(/ 0d0
     .                ,0d0
     .                ,(xmax*ymax)/(xmax*ymax
     .                  + eps*Pi*(xmax - ymax)
     -                       *Sin(Pi*((2*x1)/xmax - (2*x2)/ymax))
     -                  + eps*Pi*(xmax + ymax)
     .                       *Sin(2*Pi*(x1/xmax + x2/ymax))) /)
        case default
          write (*,*) 'Grid not implemented in covariantVector'
          write (*,*) 'Aborting...'
          stop
        end select

      end function cnv_vector_anal

c     hessian22
c     #################################################################
      function hessian22(gl_def,l,i,j,k,igx,igy,igz) result (tensor)

c     -----------------------------------------------------------------
c     Calculates hessian elements of curvilinear coordinate system in
c     the covariant basis, i.e., 
c            hessian[l](i,j) = J^2 <cnv(i)|grad(cov[l])|cnv(j)>
c     -----------------------------------------------------------------

        implicit none

c     Input variables

        integer    :: l,i,j,k,igx,igy,igz

        type(grid_mg_def),pointer :: gl_def

c     Local variables

        integer    :: ig,jg,kg
        real(8)    :: car(3),curv(3),vec(9)
        real(8)    :: x1,x2,x3,tensor(3,3),aa,bb,major_r,mm,kk,eps,phi
     .               ,gp(ngparams)

        real(8) :: xmax,ymax,zmax,xmin,ymin,zmin
! Daniele, 11/08/2017. Local variables needed for shaped toroidal coordinates
        integer :: ii
        real(8) :: ai(10),bi(10)
        real(8) :: ap,bp
        real(8) :: as,bs

c     Begin program

        xmin = gl_def%gxmin
        xmax = gl_def%gxmax

        ymin = gl_def%gymin
        ymax = gl_def%gymax

        zmin = gl_def%gzmin
        zmax = gl_def%gzmax

        gp = gl_def%params

        select case (coords)
        case ('car','lpk')
          vec = (/ 0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0 /)
        case ('cyl')
          call getCurvilinearCoordinates(gl_def,i,j,k,igx,igy,igz
     .                                  ,ig,jg,kg,x1,x2,x3)
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
          call getCurvilinearCoordinates(gl_def,i,j,k,igx,igy,igz
     .                                  ,ig,jg,kg,x1,x2,x3)
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
          call getCurvilinearCoordinates(gl_def,i,j,k,igx,igy,igz
     .                                  ,ig,jg,kg,x1,x2,x3)
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
          call getCurvilinearCoordinates(gl_def,i,j,k,igx,igy,igz
     .                                  ,ig,jg,kg,x1,x2,x3)
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

        case ('sha')
          call getCurvilinearCoordinates(gl_def,i,j,k,igx,igy,igz
     .                                  ,ig,jg,kg,x1,x2,x3)
          curv = (/ x1,x2,x3 /)

          major_r = gp(1)  !Major radius

          ai(1) = gp(2)  !Horizontal minor radius
          bi(1) = gp(3)  !Vertical   minor radius

          do ii=2,10
             ai(ii) = gp(7+ii) ! Higher order cosine coefficients are stored in gp(9) to gp(17)
             bi(ii) = gp(16+ii) ! Higher order sine coefficients are stored in gp(18) to gp(26)
          enddo

          aa = 0d0
          bb = 0d0
          ap = 0d0
          bp = 0d0
          as = 0d0
          bs = 0d0
          do ii=1,10
             aa = aa + ai(ii)*cos(ii*curv(2)) ! A(theta)=sum_i(ai*cos(i*theta))
             bb = bb + bi(ii)*sin(ii*curv(2)) ! B(theta)=sum_i(bi*sin(i*theta))
             bp = bp + bi(ii)*ii*cos(ii*curv(2)) ! B'(theta)
             ap = ap - ai(ii)*ii*sin(ii*curv(2)) ! A'(theta)
             as = as - ai(ii)*ii*ii*cos(ii*curv(2)) ! A''(theta)
             bs = bs - bi(ii)*ii*ii*sin(ii*curv(2)) ! B''(theta)
          enddo

          select case (l)
            case (1)
              vec(1) = 0d0
              vec(2) = 0d0
              vec(3) = 0d0
              vec(4) = vec(2)
              vec(5) = curv(1)*(bp*as-ap*bs)/(aa*bp-bb*ap)
              vec(6) = 0d0
              vec(7) = vec(3)
              vec(8) = vec(6)
              vec(9) = -bp*(major_r+aa*curv(1))/(aa*bp-bb*ap)
            case (2)
              vec(1) = 0d0
              vec(2) = 1./curv(1)
              vec(3) = 0d0
              vec(4) = vec(2)
              vec(5) = (aa*bs-bb*as)/(aa*bp-bb*ap)
              vec(6) = 0d0
              vec(7) = vec(3)
              vec(8) = vec(6)
              vec(9) = bb*(major_r+aa*curv(1))/curv(1)/(aa*bp-bb*ap)
            case (3)
              vec(1) = 0d0
              vec(2) = 0d0
              vec(3) = aa/(major_r+aa*curv(1))
              vec(4) = vec(2)
              vec(5) = 0d0
              vec(6) = ap*curv(1)/(major_r + aa*curv(1))
              vec(7) = vec(3)
              vec(8) = vec(6)
              vec(9) = 0d0
          end select

        case ('sin')

          call getCurvilinearCoordinates(gl_def,i,j,k,igx,igy,igz
     .                                  ,ig,jg,kg,x1,x2,x3)
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

c$$$c     gamma_anal
c$$$c     #################################################################
c$$$      function gamma_anal(gl_def,curv) result(gamma)
c$$$
c$$$c     -----------------------------------------------------------------
c$$$c     Calculates covariant vectors of curvilinear coordinate system
c$$$c     in Cartesian coordinates
c$$$c     -----------------------------------------------------------------
c$$$
c$$$        implicit none
c$$$
c$$$c     Input variables
c$$$
c$$$        real(8) :: curv(3),gamma(3,3,3)
c$$$
c$$$        type(grid_mg_def),pointer :: gl_def
c$$$
c$$$c     Local variables
c$$$
c$$$        real(8) :: aa,bb,major_r,mm,kk,eps,phi,jac,r
c$$$     .            ,vec(3,3,3),x1,x2,x3
c$$$        real(8) :: gp(ngparams)
c$$$
c$$$        real(8) :: xmax,ymax,zmax,xmin,ymin,zmin
c$$$
c$$$c     Begin program
c$$$
c$$$        xmin = gl_def%gxmin
c$$$        xmax = gl_def%gxmax
c$$$
c$$$        ymin = gl_def%gymin
c$$$        ymax = gl_def%gymax
c$$$
c$$$        zmin = gl_def%gzmin
c$$$        zmax = gl_def%gzmax
c$$$
c$$$        gp = gl_def%params
c$$$
c$$$        select case (coords)
c$$$        case ('car','lpk')
c$$$          vec = 0d0
c$$$        case ('cyl')
c$$$
c$$$          vec(1,1,1) = 0d0
c$$$          vec(1,1,2) = 0d0
c$$$          vec(1,1,3) = 0d0
c$$$          vec(1,2,1) = vec(1,1,2)
c$$$          vec(1,2,2) = curv(1)
c$$$          vec(1,2,3) = 0d0
c$$$          vec(1,3,1) = vec(1,1,3)
c$$$          vec(1,3,2) = vec(1,2,3)
c$$$          vec(1,3,3) = 0d0
c$$$
c$$$          vec(2,1,1) = 0d0
c$$$          vec(2,1,2) = -1./curv(1)
c$$$          vec(2,1,3) = 0d0
c$$$          vec(2,2,1) = vec(2,1,2)
c$$$          vec(2,2,2) = 0d0
c$$$          vec(2,2,3) = 0d0
c$$$          vec(2,3,1) = vec(2,1,3)
c$$$          vec(2,3,2) = vec(2,2,3)
c$$$          vec(2,3,3) = 0d0
c$$$
c$$$          vec(3,1,1) = 0d0
c$$$          vec(3,1,2) = 0d0
c$$$          vec(3,1,3) = 0d0
c$$$          vec(3,2,1) = vec(3,1,2)
c$$$          vec(3,2,2) = 0d0
c$$$          vec(3,2,3) = 0d0
c$$$          vec(3,3,1) = vec(3,1,3)
c$$$          vec(3,3,2) = vec(3,2,3)
c$$$          vec(3,3,3) = 0d0
c$$$
c$$$          vec = -vec
c$$$
c$$$        case ('hel')
c$$$
c$$$          mm = gp(1)
c$$$          kk = gp(2)
c$$$          aa = kk/mm
c$$$          phi = (curv(2)-aa*curv(3))
c$$$
c$$$          vec(1,1,1) = 0d0
c$$$          vec(1,1,2) = 0d0
c$$$          vec(1,1,3) = 0d0
c$$$          vec(1,2,1) = vec(1,1,2)
c$$$          vec(1,2,2) = curv(1)
c$$$          vec(1,2,3) = -aa*curv(1)
c$$$          vec(1,3,1) = vec(1,1,3)
c$$$          vec(1,3,2) = vec(1,2,3)
c$$$          vec(1,3,3) = aa**2*curv(1)
c$$$                    
c$$$          vec(2,1,1) = 0d0
c$$$          vec(2,1,2) = -1./curv(1)
c$$$          vec(2,1,3) = aa/curv(1)
c$$$          vec(2,2,1) = vec(2,1,2)
c$$$          vec(2,2,2) = 0d0
c$$$          vec(2,2,3) = 0d0
c$$$          vec(2,3,1) = vec(2,1,3)
c$$$          vec(2,3,2) = vec(2,2,3)
c$$$          vec(2,3,3) = 0d0
c$$$                    
c$$$          vec(3,1,1) = 0d0
c$$$          vec(3,1,2) = 0d0
c$$$          vec(3,1,3) = 0d0
c$$$          vec(3,2,1) = vec(3,1,2)
c$$$          vec(3,2,2) = 0d0
c$$$          vec(3,2,3) = 0d0
c$$$          vec(3,3,1) = vec(3,1,3)
c$$$          vec(3,3,2) = vec(3,2,3)
c$$$          vec(3,3,3) = 0d0
c$$$
c$$$          vec = -vec
c$$$
c$$$cc        case ('hl2')
c$$$cc
c$$$cc          curv = (/ abs(2*x1),x2,x3 /)
c$$$cc
c$$$cc          mm = gp(1)
c$$$cc          kk = gp(2)
c$$$cc          aa = kk/mm
c$$$cc          phi = (curv(2)-aa*curv(3))
c$$$cc          select case (l)
c$$$cc            case (1)
c$$$cc              vec(1) = -1d0/curv(1)
c$$$cc              vec(2) = 0d0
c$$$cc              vec(3) = 0d0
c$$$cc              vec(4) = vec(2)
c$$$cc              vec(5) = -curv(1)
c$$$cc              vec(6) = aa*curv(1)
c$$$cc              vec(7) = vec(3)
c$$$cc              vec(8) = vec(6)
c$$$cc              vec(9) = -aa**2*curv(1)
c$$$cc            case (2)
c$$$cc              vec(1) = 0d0
c$$$cc              vec(2) = 1./curv(1)
c$$$cc              vec(3) = -aa/curv(1)
c$$$cc              vec(4) = vec(2)
c$$$cc              vec(5) = 0d0
c$$$cc              vec(6) = 0d0
c$$$cc              vec(7) = vec(3)
c$$$cc              vec(8) = vec(6)
c$$$cc              vec(9) = 0d0
c$$$cc            case (3)
c$$$cc              vec(1) = 0d0
c$$$cc              vec(2) = 0d0
c$$$cc              vec(3) = 0d0
c$$$cc              vec(4) = vec(2)
c$$$cc              vec(5) = 0d0
c$$$cc              vec(6) = 0d0
c$$$cc              vec(7) = vec(3)
c$$$cc              vec(8) = vec(6)
c$$$cc              vec(9) = 0d0
c$$$cc          end select
c$$$
c$$$        case ('tor')
c$$$
c$$$          major_r = gp(1)
c$$$          aa      = gp(2)
c$$$          bb      = gp(3)
c$$$
c$$$          vec(1,1,1) = 0d0
c$$$          vec(1,1,2) = 0d0
c$$$          vec(1,1,3) = 0d0
c$$$          vec(1,2,1) = vec(1,1,2)
c$$$          vec(1,2,2) =-curv(1)
c$$$          vec(1,2,3) = 0d0
c$$$          vec(1,3,1) = vec(1,1,3)
c$$$          vec(1,3,2) = vec(1,2,3)
c$$$          vec(1,3,3) =-cos(curv(2))*(major_r+aa*curv(1)*cos(curv(2)))/aa
c$$$                    
c$$$          vec(2,1,1) = 0d0
c$$$          vec(2,1,2) = 1./curv(1)
c$$$          vec(2,1,3) = 0d0
c$$$          vec(2,2,1) = vec(2,1,2)
c$$$          vec(2,2,2) = 0d0
c$$$          vec(2,2,3) = 0d0
c$$$          vec(2,3,1) = vec(2,1,3)
c$$$          vec(2,3,2) = vec(2,2,3)
c$$$          vec(2,3,3) = (major_r + aa*curv(1)*cos(curv(2)))*sin(curv(2))
c$$$     .                 /curv(1)/aa
c$$$                    
c$$$          vec(3,1,1) = 0d0
c$$$          vec(3,1,2) = 0d0
c$$$          vec(3,1,3) = aa*cos(curv(2))/(major_r+aa*curv(1)*cos(curv(2)))
c$$$          vec(3,2,1) = vec(3,2,1)
c$$$          vec(3,2,2) = 0d0
c$$$          vec(3,2,3) = -aa*curv(1)*sin(curv(2))
c$$$     .                 /(major_r + aa*curv(1)*cos(curv(2)))
c$$$          vec(3,3,1) = vec(3,3,1)
c$$$          vec(3,3,2) = vec(3,3,2)
c$$$          vec(3,3,3) = 0d0
c$$$
c$$$        case ('sin')
c$$$
c$$$          eps = gp(1)
c$$$
c$$$          x1 = curv(1)
c$$$          x2 = curv(2)
c$$$          x3 = curv(3)
c$$$
c$$$          vec(1,1,1) = (-4*eps*Pi**2*ymax*Sin((2*Pi*x1)/xmax)*
c$$$     -                                    Sin((2*Pi*x2)/ymax))/
c$$$     -          (xmax*(xmax*ymax + eps*Pi*(xmax - ymax)*
c$$$     -           Sin(Pi*((2*x1)/xmax - (2*x2)/ymax)) + 
c$$$     -           eps*Pi*(xmax + ymax)*Sin(2*Pi*(x1/xmax + x2/ymax))))
c$$$          vec(1,1,2) = (4*eps*Pi**2*Cos((2*Pi*x1)/xmax)
c$$$     .                             *Cos((2*Pi*x2)/ymax))/
c$$$     -            (xmax*ymax + eps*Pi*(xmax - ymax)*
c$$$     -              Sin(Pi*((2*x1)/xmax - (2*x2)/ymax)) + 
c$$$     -            eps*Pi*(xmax + ymax)*Sin(2*Pi*(x1/xmax + x2/ymax)))
c$$$          vec(1,1,3) = 0d0
c$$$          vec(1,2,1) = vec(1,1,2)
c$$$          vec(1,2,2) = (-4*eps*Pi**2*xmax*Sin((2*Pi*x1)/xmax)
c$$$     .                                   *Sin((2*Pi*x2)/ymax))/
c$$$     -            (ymax*(xmax*ymax + eps*Pi*(xmax - ymax)*
c$$$     -            Sin(Pi*((2*x1)/xmax - (2*x2)/ymax)) + 
c$$$     -            eps*Pi*(xmax + ymax)*Sin(2*Pi*(x1/xmax + x2/ymax))))
c$$$          vec(1,2,3) = 0d0
c$$$          vec(1,3,1) = vec(1,1,3)
c$$$          vec(1,3,2) = vec(1,2,3)
c$$$          vec(1,3,3) = 0d0
c$$$
c$$$          vec(2,1,1) = (-4*eps*Pi**2*ymax*Sin((2*Pi*x1)/xmax)
c$$$     -                                   *Sin((2*Pi*x2)/ymax))/
c$$$     -                 (xmax*(xmax*ymax
c$$$     .                 + eps*Pi*(xmax - ymax)
c$$$     -                     *Sin(Pi*((2*x1)/xmax - (2*x2)/ymax))
c$$$     -                 + eps*Pi*(xmax + ymax)
c$$$     .                     *Sin(2*Pi*(x1/xmax + x2/ymax))))
c$$$          vec(2,1,2) = (4*eps*Pi**2*Cos((2*Pi*x1)/xmax)
c$$$     .                             *Cos((2*Pi*x2)/ymax))/
c$$$     -                 (xmax*ymax
c$$$     .                 + eps*Pi*(xmax - ymax)
c$$$     -                     *Sin(Pi*((2*x1)/xmax - (2*x2)/ymax))
c$$$     -                 + eps*Pi*(xmax + ymax)
c$$$     .                     *Sin(2*Pi*(x1/xmax + x2/ymax)))
c$$$          vec(2,1,3) = 0d0
c$$$          vec(2,2,1) = vec(2,1,2)
c$$$          vec(2,2,2) = (-4*eps*Pi**2*xmax*Sin((2*Pi*x1)/xmax)
c$$$     .                                   *Sin((2*Pi*x2)/ymax))/
c$$$     -                  (ymax*(xmax*ymax
c$$$     .                  + eps*Pi*(xmax - ymax)
c$$$     -                      *Sin(Pi*((2*x1)/xmax - (2*x2)/ymax))
c$$$     -                  + eps*Pi*(xmax + ymax)
c$$$     .                      *Sin(2*Pi*(x1/xmax + x2/ymax))))
c$$$          vec(2,2,3) = 0d0
c$$$          vec(2,3,1) = vec(2,1,3)
c$$$          vec(2,3,2) = vec(2,2,3)
c$$$          vec(2,3,3) = 0d0
c$$$
c$$$          vec(3,1,1) = 0d0
c$$$          vec(3,1,2) = 0d0
c$$$          vec(3,1,3) = 0d0
c$$$          vec(3,2,1) = vec(3,1,2)
c$$$          vec(3,2,2) = 0d0
c$$$          vec(3,2,3) = 0d0
c$$$          vec(3,3,1) = vec(3,1,3)
c$$$          vec(3,3,2) = vec(3,2,3)
c$$$          vec(3,3,3) = 0d0
c$$$
c$$$        case default
c$$$          write (*,*) 'Grid not implemented in hessian'
c$$$          write (*,*) 'Aborting...'
c$$$          stop
c$$$        end select
c$$$
c$$$        gamma = vec
c$$$
c$$$      end function gamma_anal

c     christ_2knd
c     #################################################################
      function christ_2knd(gl_def,i,j,k,igx,igy,igz) result (tensor)

c     -----------------------------------------------------------------
c     Calculates elements of Christoffel symbol of the second kind,
c     Gamma[l](i,j) defined as: 
c            gamma[l](i,j) = -J^2 <cnv(i)|grad(cov[l])|cnv(j)>
c     -----------------------------------------------------------------

        implicit none

c     Input variables

        integer    :: l,i,j,k,igx,igy,igz
        real(8)    :: tensor(3,3,3)

        type(grid_mg_def),pointer :: gl_def

c     Local variables

c     Begin program

        tensor(1,:,:) = hessian22(gl_def,1,i,j,k,igx,igy,igz)
        tensor(2,:,:) = hessian22(gl_def,2,i,j,k,igx,igy,igz)
        tensor(3,:,:) = hessian22(gl_def,3,i,j,k,igx,igy,igz)

c     End program

      end function christ_2knd

c     hessian_cnv
c     #################################################################
      function hessian_cnv(gl_def,k,x1,x2,x3) result (tensor)

c     -----------------------------------------------------------------
c     Calculates elements of tensor grad(cnv) in a mixed coordinate 
c     system:
c              hessian_cnv[k](i,j) = J^2 <cnv(i)|grad(cnv[k])|cov(j)>
c     -----------------------------------------------------------------

        implicit none

c     Input variables

        integer    :: k
        real(8)    :: x1,x2,x3,tensor(3,3)

        type(grid_mg_def),pointer :: gl_def

c     Local variables

        real(8)    :: vec(9),car(3),aa,bb,major_r,mm,kk,eps,gp(ngparams)

c     Begin program

        gp = gl_def%params

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
      function g_sub(gl_def,i,j,k,igx,igy,igz,vertex) result (tensor)

c     -----------------------------------------------------------------
c     Calculates contravariant metric tensor of curvilinear coordinate 
c     system
c     -----------------------------------------------------------------

        implicit none

c     Input variables

        integer :: i,j,k,igx,igy,igz
        real(8) :: tensor(3,3)
        logical,optional :: vertex

        type(grid_mg_def),pointer :: gl_def

c     Local variables

        integer :: ig,jg,kg
        real(8) :: vec(9),jac,car(3),curv(3)
        real(8) :: x1,x2,x3,aa,bb,major_r,mm,kk,eps,gp(ngparams)
     .            ,x1p,x2p,x3p,i_x1
        logical :: vrtx

        real(8) :: xmax,ymax,zmax,xmin,ymin,zmin
! Daniele, 11/08/2017. Local variables needed for shaped toroidal coordinates
        integer :: ii
        real(8) :: ai(10),bi(10)
        real(8) :: ap,bp

c     Begin program

        xmin = gl_def%gxmin
        xmax = gl_def%gxmax

        ymin = gl_def%gymin
        ymax = gl_def%gymax

        zmin = gl_def%gzmin
        zmax = gl_def%gzmax

        if (PRESENT(vertex)) then
          vrtx = vertex
        else
          vrtx = .false.
        endif

        gp = gl_def%params

        call getCurvilinearCoordinates(gl_def,i,j,k,igx,igy,igz,ig,jg,kg
     .                                ,x1,x2,x3)

        if (vrtx) then
          call getCurvilinearCoordinates(gl_def,i+1,j+1,k+1,igx,igy,igz
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
        case ('sha')
! Daniele, 11/08/2017
! Shaped toroidal coordinates system
! x = (R0 + r*sum_i(ai*cos(i*theta))*cos(-phi)
! y = (R0 + r*sum_i(ai*cos(i*theta))*sin(-phi)
! z = r*sum_i(bi*sin(i*theta)
          if (x1 == 0d0) then
            vec = 0d0
          else
          major_r = gp(1)  !Major radius

          ai(1) = gp(2)  !Horizontal minor radius
          bi(1) = gp(3)  !Vertical   minor radius

          do ii=2,10
             ai(ii) = gp(7+ii) ! Higher order cosine coefficients are stored in gp(9) to gp(17)
             bi(ii) = gp(16+ii) ! Higher order sine coefficients are stored in gp(18) to gp(26)
          enddo

          aa = 0d0
          bb = 0d0
          ap = 0d0
          bp = 0d0
          do ii=1,10
             aa = aa + ai(ii)*cos(ii*x2) ! A(theta)=sum_i(ai*cos(i*theta))
             bb = bb + bi(ii)*sin(ii*x2) ! B(theta)=sum_i(bi*sin(i*theta))
             bp = bp + bi(ii)*ii*cos(ii*x2) ! B'(theta)
             ap = ap - ai(ii)*ii*sin(ii*x2) ! A'(theta)
          enddo

            vec(1) = (aa**2 + bb**2)/
     .               (x1*(major_r + aa*x1)*(aa*bp-bb*ap))

            vec(2) = (aa*ap + bb*bp)/
     .               ((major_r + aa*x1)*(aa*bp-bb*ap))

            vec(3) = 0d0

            vec(4) = vec(2)

            vec(5) = x1*(ap**2 + bp**2)/
     .               ((major_r + aa*x1)*(aa*bp-bb*ap))

            vec(6) = 0d0

            vec(7) = vec(3)
            vec(8) = vec(6)
            vec(9) = (major_r + aa*x1)/x1/(aa*bp-bb*ap)

          endif
        case ('sin')

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
      function g_sup(gl_def,i,j,k,igx,igy,igz,vertex) result (tensor)

c     -----------------------------------------------------------------
c     Calculates covariant metric tensor of curvilinear coordinate 
c     system
c     -----------------------------------------------------------------

        implicit none

c     Input variables

        integer :: i,j,k,igx,igy,igz
        real(8) :: tensor(3,3)
        logical,optional :: vertex

        type(grid_mg_def),pointer :: gl_def

c     Local variables

        integer :: ig,jg,kg
        real(8) :: vec(9),jac,car(3),curv(3)
        real(8) :: x1,x2,x3,aa,bb,major_r,mm,kk,eps,gp(ngparams)
     .            ,x1p,x2p,x3p
        logical :: vrtx

        real(8) :: xmax,ymax,zmax,xmin,ymin,zmin
! Daniele, 11/08/2017. Local variables needed for shaped toroidal coordinates
        integer :: ii
        real(8) :: ai(10),bi(10)
        real(8) :: ap,bp

c     Begin program

        xmin = gl_def%gxmin
        xmax = gl_def%gxmax

        ymin = gl_def%gymin
        ymax = gl_def%gymax

        zmin = gl_def%gzmin
        zmax = gl_def%gzmax

        if (PRESENT(vertex)) then
          vrtx = vertex
        else
          vrtx = .false.
        endif

        gp = gl_def%params

        call getCurvilinearCoordinates(gl_def,i,j,k,igx,igy,igz,ig,jg,kg
     .                                ,x1,x2,x3)

        if (vrtx) then
          call getCurvilinearCoordinates(gl_def,i+1,j+1,k+1,igx,igy,igz
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
        case ('sha')
! Daniele, 11/08/2017
! Shaped toroidal coordinates system
! x = (R0 + r*sum_i(ai*cos(i*theta))*cos(-phi)
! y = (R0 + r*sum_i(ai*cos(i*theta))*sin(-phi)
! z = r*sum_i(bi*sin(i*theta)
          if (x1 == 0d0) then
            vec = 0d0
          else
          major_r = gp(1)  !Major radius

          ai(1) = gp(2)  !Horizontal minor radius
          bi(1) = gp(3)  !Vertical   minor radius

          do ii=2,10
             ai(ii) = gp(7+ii) ! Higher order cosine coefficients are stored in gp(9) to gp(17)
             bi(ii) = gp(16+ii) ! Higher order sine coefficients are stored in gp(18) to gp(26)
          enddo

          aa = 0d0
          bb = 0d0
          ap = 0d0
          bp = 0d0
          do ii=1,10
             aa = aa + ai(ii)*cos(ii*x2) ! A(theta)=sum_i(ai*cos(i*theta))
             bb = bb + bi(ii)*sin(ii*x2) ! B(theta)=sum_i(bi*sin(i*theta))
             bp = bp + bi(ii)*ii*cos(ii*x2) ! B'(theta)
             ap = ap - ai(ii)*ii*sin(ii*x2) ! A'(theta)
          enddo

            vec(1) = x1*(major_r + aa*x1)*(ap**2 + bp**2)/
     .               (aa*bp-bb*ap)

            vec(2) = (major_r + aa*x1)*(-aa*ap - bb*bp)/
     .               (aa*bp-bb*ap)

            vec(3) = 0d0

            vec(4) = vec(2)

            vec(5) = (major_r + aa*x1)*(aa**2 + bb**2)/
     .               (x1*(aa*bp-bb*ap))

            vec(6) = 0d0

            vec(7) = vec(3)
            vec(8) = vec(6)
            vec(9) = x1*(aa*bp-bb*ap)/(major_r + aa*x1)

          endif
        case ('sin')

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

c     gsup_anal
c     #################################################################
      function gsup_anal(gl_def,curv) result(tensor)

c     -----------------------------------------------------------------
c     Calculates covariant vectors of curvilinear coordinate system
c     in Cartesian coordinates
c     -----------------------------------------------------------------

        implicit none

c     Input variables

        real(8) :: curv(3),tensor(3,3)

        type(grid_mg_def),pointer :: gl_def

c     Local variables

        real(8) :: aa,bb,major_r,mm,kk,eps,phi,jac,x1,x2,x3
        real(8) :: gp(ngparams),vec(9)

        real(8) :: xmax,ymax,zmax,xmin,ymin,zmin

c     Begin program

        gp = gl_def%params

        x1 = curv(1)
        x2 = curv(2)
        x3 = curv(3)

        xmin = gl_def%gxmin
        xmax = gl_def%gxmax

        ymin = gl_def%gymin
        ymax = gl_def%gymax

        zmin = gl_def%gzmin
        zmax = gl_def%gzmax

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

      end function gsup_anal

      end module grid_anal_map
