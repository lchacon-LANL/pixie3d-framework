
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
     .      ,'sha','sh2','sh3','stl','tsq')

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
        real(8) :: xmax,ymax,zmax,xmin,ymin,zmin,shftR,shftZ
     .            ,kappa,delta,zeta,cs,sn,Z0,kappa_sp
! Daniele, 11/07/2017. Local variables needed for shaped toroidal coordinates
        integer :: i
        real(8) :: ai(10),bi(10)
! Daniele, 11/08/2017. Local variables needed for stellarator (helias) coordinates
        real(8) :: nn,uu,vv,d10,d0,d20,d1m1,d2m1,d2m2,rc,zc

c$$$        logical,save :: dump=.true.
        
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
          rr = gp(3)

          if (abs(abs(kk*rr)-nint(abs(kk*rr))*1d0) > 1d-13) then
            write (*,*) "Inconsistent helical input: k*R=",abs(kk*rr)
     $                 ,"not an integer"
            write (*,*) "Aborting..."
            stop
          endif

          phi = (xlog(2)-kk*xlog(3))/mm
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

        case ('sh3')
! Shaped toroidal coordinates system from X. Z. Tang
! x = (R0 + r*a*(cos(theta+asin(delta*r**2*sin(theta)))-dR))*cos(-phi)
! y = (R0 + r*a*(cos(theta+asin(delta*r**2*sin(theta)))-dR))*sin(-phi)
! z =  Z0 + r*kappa*(sin(theta + zeta*r**2*sin(2*theta))-dZ)

          major_r = gp(1)  !R-coord of magnetic axis (major radius)
          Z0      = gp(2)  !Z-coord of magnetic axis
          aa      = gp(3)  !Horizontal minor radius
          kappa   = gp(4)  !Elongation
          delta   = gp(5)  !Triangularity
          zeta    = gp(6)  !Squareness
          kappa_sp= gp(7)  !Elongation at magnetic axis
          shftR   = gp(8)  !Origin shift in R
          shftZ   = gp(9)  !Origin shift in Z

          if (kappa_sp > 0d0) then
            bb = aa*((1d0-xlog(1))*kappa_sp + xlog(1)*kappa)
          else
            bb = aa*kappa
          endif
          
          !Scale shift to account for ellipticity and to decay with r->0
          shftR = shftR*xlog(1)/aa
          shftZ = shftZ*xlog(1)/bb
          
          cs = cos(xlog(2)+asin(delta*xlog(1)**2*sin(xlog(2))))-shftR
          sn = sin(xlog(2)+zeta*xlog(1)**2*sin(2*xlog(2)))     -shftZ
          
          car(1)=(major_r + aa*xlog(1)*cs)*cos(-xlog(3))
          car(2)=(major_r + aa*xlog(1)*cs)*sin(-xlog(3))
          car(3)= Z0      + bb*xlog(1)*sn

c$$$          if (dump) then
c$$$            dump = .false.
c$$$            write (*,*) gp
c$$$            write (*,*) "EFITdiag"
c$$$            do i=1,100
c$$$              xlog(2) = 2*pi/100*i
c$$$              bb = aa*kappa
c$$$              cs = cos(xlog(2)+asin(delta*sin(xlog(2))))-gp(8)/aa
c$$$              sn = sin(xlog(2)+zeta*sin(2*xlog(2)))     -gp(9)/bb
c$$$              write (*,*) "EFITdiag",(major_r + aa*cs),Z0+bb*sn
c$$$            enddo
c$$$          endif
c$$$          stop
        case ('tsq')            !Square torus
          rc = gp(1)  !Domain center R-coord
          zc = gp(2)  !Domain center Z-coord
          aa = gp(3)  !Horizontal thickness
          bb = gp(4)  !Vertical   thickness

          car(1) = (rc + aa*(xlog(1)-0.5))*cos(-xlog(3))
          car(2) = (rc + aa*(xlog(1)-0.5))*sin(-xlog(3))
          car(3) =  zc + bb*(xlog(2)-0.5)
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

c     xi_x_sp
c     #################################################################
      subroutine xi_x_sp(g_def,xlog,xcar,ghost_cell,polar)

c     -----------------------------------------------------------------
c     Find inverse map xi(x) around SP.
c     -----------------------------------------------------------------

      implicit none

c     Input variables

      real(8) :: xcar(3),xlog(3)
      logical :: ghost_cell,polar

      type(grid_mg_def),pointer :: g_def

c     Local variables

      integer :: ig,jg,kg
      real(8) :: x,y,z,gp(ngparams),mm,kk,aa,bb,R0,R1,R2,dR2,rr,th,zz
     .          ,Jac,rhs,da,db,Z2,dZ2,iaa,ibb,da_a,db_b,rold,Rs,gc_cf
     .          ,phi

      integer :: m,n
      real(8) :: ai(10),bi(10)

c     Begin program

      gc_cf = 1d0
      if (ghost_cell) gc_cf = -1d0

      x = xcar(1) ; y = xcar(2) ; z = xcar(3)
      
c     Find radius and angle

      select case (coords)
      case('hel')

        if (polar) then
          rr = gc_cf*sqrt(x*x+y*y)
          th = atan2(y,x)
          if (rr < 0d0) th = th + pi
          zz = z
          th = mod(th+2*pi,2*pi) !Ensure angle is in [0,2pi]
        else
          rr = gc_cf*sqrt(x*x+y*y)
          phi = atan2(y,x)      !Physical angle
          if (rr < 0d0) phi = phi + pi
          mm = g_def%params(1)  !Poloidal mode number
          kk = g_def%params(2)  !Toroidal wavenumber
          th = (mm*phi + kk*z)  !Logical angle
          zz = z
          th = mod(th+2*pi,2*pi) !Ensure angle is in [0,2pi]
        endif
        
      case('cyl')

        rr = gc_cf*sqrt(x*x+y*y)
        th = atan2(y,x)
        if (rr < 0d0) th = th + pi
        zz = z
        th = mod(th+2*pi,2*pi)   !Ensure angle is in [0,2pi]

      case('tor')

        R0 = g_def%params(1) !Major radius
        aa = g_def%params(2) !Horizontal radius
        bb = g_def%params(3) !Vertical   radius

        R1 = (sqrt(x*x+y*y)-R0)/aa          !In-plane radius
        
        rr = gc_cf*sqrt(R1*R1+z*z/(bb*bb))  !Poloidal radius
        zz =-atan2(y,x)                     !Toroidal angle

        th = atan2(z/bb,R1)                 !Poloidal angle
        
        if (rr < 0d0) th = th + pi

        th = mod(th+2*pi,2*pi)   !Ensure angle is in [0,2pi]
        zz = mod(zz+2*pi,2*pi)   !Ensure angle is in [0,2pi]

      case('sha')
! Shaped toroidal coordinates system
! x = (R0 + r*sum_i(ai*cos(i*theta))*cos(-phi)
! y = (R0 + r*sum_i(ai*cos(i*theta))*sin(-phi)
! z = r*sum_i(bi*sin(i*theta)

        R0    = g_def%params(1) !Major radius
        ai(1) = g_def%params(2) !Horizontal radius
        bi(1) = g_def%params(3) !Vertical   radius

        do m=2,10
          ai(m) = g_def%params(7+m)  ! Higher order cosine coefficients are stored in gp(9) to gp(17)
          bi(m) = g_def%params(16+m) ! Higher order sine coefficients are stored in gp(18) to gp(26)
        enddo

        !Poloidal angle
        th = xlog(2)      !Don't expect packing in poloidal angle, so we can trust logical coord
        
        !Toroidal angle
        zz =-atan2(y,x)   !Or x3, for same reason

        !Poloidal radius
        aa = 0d0 ; bb = 0d0
        do m=1,10
          aa = aa + ai(m)*cos(m*th)
          bb = bb + bi(m)*sin(m*th)
        enddo
        aa = aa/cos(th)
        bb = bb/sin(th)
        
        R1 = ((sqrt(x*x+y*y)-R0)/aa)          !In-plane radius
        
        rr = gc_cf*sqrt(R1*R1+z*z/(bb*bb))
        
        if (rr < 0d0) th = th + pi

      case('sh2')
! Shaped toroidal coordinates system
! x = (R0 + r*a1*cos(theta)+r^2*sum_i(ai*cos(i*theta))*cos(-phi)
! y = (R0 + r*a1*cos(theta)+r^2*sum_i(ai*cos(i*theta))*sin(-phi)
!! x = (R0 + r*sum_i(r^(i-1)*ai*cos(i*theta))*cos(-phi)
!! y = (R0 + r*sum_i(r^(i-1)*ai*cos(i*theta))*sin(-phi)
! z = r*sum_i(bi*sin(i*theta) <<< INCORRECT???

        R0    = g_def%params(1) !Major radius
        ai(1) = g_def%params(2) !Horizontal radius
        bi(1) = g_def%params(3) !Vertical   radius

        do m=2,10
          ai(m) = g_def%params(7+m)  ! Higher order cosine coefficients are stored in gp(9) to gp(17)
          bi(m) = g_def%params(16+m) ! Higher order sine coefficients are stored in gp(18) to gp(26)
        enddo

        !Poloidal angle
        th = xlog(2)           !Don't expect packing in poloidal angle, so we can trust logical coord
        
        !Toroidal angle
        zz =-atan2(y,x)   !Or x3, for same reason

        !Poloidal radius (requires iteration!)
        Rs = (sqrt(x*x+y*y)-R0)  !Save for efficiency

        aa = ai(1)
        bb = bi(1)
        R2 = Rs/aa ; R2 = R2*R2  !In-plane radius
        Z2 = z/bb  ; Z2 = Z2*Z2
        rr = sqrt(R2+Z2)

c$$$        do n=1,30   !Picard iteration
c$$$          rold = rr
c$$$          aa = ai(1)*cos(th)
c$$$          bb = bi(1)*sin(th)
c$$$          do m=2,10
c$$$            aa = aa + rr*ai(m)*cos(m*th)
c$$$            bb = bb + rr*bi(m)*sin(m*th)  !THIS IS NOT WHAT IS STATED IN COMMENTS ABOVE
c$$$          enddo
c$$$          iaa = cos(th)/aa
c$$$          ibb = sin(th)/bb
c$$$
c$$$          !Compute Jacobian
c$$$          R2 = Rs*iaa ; R2 = R2*R2 !In-plane radius
c$$$          Z2 = z*z*ibb*ibb
c$$$          rr = sqrt(R2+Z2)
c$$$          write (*,*) "iter",n,"conv",abs(rr-rold)
c$$$          if (abs(rr-rold) < 1d-5) exit
c$$$        enddo
c$$$        write (*,*)

        do n=1,30  !Newton iteration
          aa = ai(1)*cos(th)
          bb = bi(1)*sin(th)
          da = 0d0
          db = 0d0
          do m=2,10
            aa = aa + rr*ai(m)*cos(m*th)
            bb = bb + rr*bi(m)*sin(m*th)  !THIS IS NOT WHAT IS STATED IN COMMENTS ABOVE
            da = da +    ai(m)*cos(m*th)
            db = db +    bi(m)*sin(m*th)
          enddo
          iaa = 1d0/aa
          ibb = 1d0/bb
          da_a = da*iaa
          db_b = db*ibb
          iaa = cos(th)*iaa
          ibb = sin(th)*ibb

          !Compute Jacobian
          R2 = Rs*iaa ; R2 = R2*R2 !In-plane radius
          Z2 = z*ibb  ; Z2 = Z2*Z2 !Z coord
          dR2 = -2*R2*da_a
          dZ2 = -2*Z2*db_b
          Jac = (2*rr - dR2 - dZ2) !(1d0-0.5/sqrt(R2+Z2)*(dR2+dZ2))

          !Compute residual
          rhs = rr*rr - R2 - Z2  !rr-sqrt(R2+Z2)

          !Update solution
          rr = rr - rhs/Jac
cc          write (*,*) n,"rr",rr,"res",rhs,"Jac",Jac,da/aa,db/bb
          if (abs(rhs) < 1d-5) exit
        enddo
cc        write (*,*)

        if (n > 30) then
          write (*,*) "Coordinate iteration did not converge"
          write (*,*) "Aborting in xi_x_sp..."
          stop
        endif

        rr = gc_cf*rr

        if (rr < 0d0) th = th + pi

      case default
        write (*,*) "Coordinates "//coords//" not implemented"
        write (*,*) "Aborting in xi_x_sp..."
        stop
      end select

c     Find logical coords from (rr,th,zz) by inverting packing map (if appropriate)

      if ((.not.polar).and.(coords /= 'lpk')) then
        xlog(1) = scale_coord_inv(g_def,1,rr)
        xlog(2) = scale_coord_inv(g_def,2,th)
        xlog(3) = scale_coord_inv(g_def,3,zz)
      else
        xlog = (/ rr, th, zz /)
      endif

      end subroutine xi_x_sp

c     find_polar_coords
c     ########################################################################
      subroutine find_polar_coords(g_def,igrid,i,j,k,rr,th,zz,coords)

c     ------------------------------------------------------------------------
c     In helical, cylindrical, and toroidal geometries, it finds polar coords
c     (r,theta,z) or (r,theta,phi) from Cartesian components. This routine
c     is needed because, when packing is used in radius, the logical coordinate
c     xi no longer corresponds to the real radial coordinate, r(xi).
c     ------------------------------------------------------------------------
      
      implicit none

c     Call variables

      type(grid_mg_def),pointer :: g_def

      integer,INTENT(IN)  :: i,j,k,igrid
      real(8),INTENT(OUT) :: rr,th,zz

      character(3) :: coords

c     Local variables

      integer :: ig,jg,kg
      real(8) :: xlog(3),xcar(3)

c     Begin program

      call getCurvilinearCoordinates(g_def,i,j,k
     .                              ,igrid,igrid,igrid
     .                              ,ig,jg,kg,xlog(1),xlog(2),xlog(3))

      call getCartesianCoordinates(g_def,i,j,k
     .                            ,igrid,igrid,igrid
     .                            ,ig,jg,kg,xcar(1),xcar(2),xcar(3))

      call xi_x_sp(g_def,xlog,xcar,(xlog(1)<0d0),.true.)

      rr = xlog(1)
      th = xlog(2)
      zz = xlog(3)
 
c     End program

      end subroutine find_polar_coords

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
          call getCurvilinearCoordinates(gl_def,i+1,j,k,igx,igy,igz
     .                                  ,ig,jg,kg,x1p,x2p,x3p)
          x1 = x1+x1p
          x2 = x2+x2p
          x3 = x3+x3p
          
          call getCurvilinearCoordinates(gl_def,i,j+1,k,igx,igy,igz
     .                                  ,ig,jg,kg,x1p,x2p,x3p)
          x1 = x1+x1p
          x2 = x2+x2p
          x3 = x3+x3p

          call getCurvilinearCoordinates(gl_def,i,j,k+1,igx,igy,igz
     .                                  ,ig,jg,kg,x1p,x2p,x3p)
          x1 = x1+x1p
          x2 = x2+x2p
          x3 = x3+x3p

          call getCurvilinearCoordinates(gl_def,i+1,j+1,k,igx,igy,igz
     .                                  ,ig,jg,kg,x1p,x2p,x3p)
          x1 = x1+x1p
          x2 = x2+x2p
          x3 = x3+x3p

          call getCurvilinearCoordinates(gl_def,i+1,j,k+1,igx,igy,igz
     .                                  ,ig,jg,kg,x1p,x2p,x3p)
          x1 = x1+x1p
          x2 = x2+x2p
          x3 = x3+x3p

          call getCurvilinearCoordinates(gl_def,i,j+1,k+1,igx,igy,igz
     .                                  ,ig,jg,kg,x1p,x2p,x3p)
          x1 = x1+x1p
          x2 = x2+x2p
          x3 = x3+x3p

          call getCurvilinearCoordinates(gl_def,i+1,j+1,k+1,igx,igy,igz
     .                                  ,ig,jg,kg,x1p,x2p,x3p)
          x1 = x1+x1p
          x2 = x2+x2p
          x3 = x3+x3p
          
          x1 = 0.125*x1
          x2 = 0.125*x2
          x3 = 0.125*x3
c$$$          call getCurvilinearCoordinates(gl_def,i+1,j+1,k+1,igx,igy,igz
c$$$     .                                  ,ig,jg,kg,x1p,x2p,x3p)
c$$$          x1 = 0.5*(x1+x1p)
c$$$          x2 = 0.5*(x2+x2p)
c$$$          x3 = 0.5*(x3+x3p)
        endif

        if (coords /= 'lpk') then
          curv(1) = scale_coord(gl_def,1,x1)
          curv(2) = scale_coord(gl_def,2,x2)
          curv(3) = scale_coord(gl_def,3,x3)
        else
          curv = (/ x1, x2, x3 /)
        endif

        select case (coords)
        case ('car','lpk')
          jac = 1d0
        case ('cyl')
          jac = curv(1)
        case ('hel')
          mm = gp(1)
          jac = curv(1)/mm
        case ('hl2')
          jac = 1d0
        case ('tor')
          major_r = gp(1)
          aa      = gp(2)
          bb      = gp(3)

          jac = aa*bb*curv(1)*(major_r + aa*curv(1)*cos(curv(2)))
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
             aa = aa + ai(ii)*cos(ii*curv(2)) ! A(theta)=sum_i(ai*cos(i*theta))
             bb = bb + bi(ii)*sin(ii*curv(2)) ! B(theta)=sum_i(bi*sin(i*theta))
             bp = bp + bi(ii)*ii*cos(ii*curv(2)) ! B'(theta)
             ap = ap - ai(ii)*ii*sin(ii*curv(2)) ! A'(theta)
          enddo

          jac = (aa*bp-bb*ap)*curv(1)*(major_r + aa*curv(1))

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

        if (coords /= 'lpk') then
          jac = jac*scale_coord(gl_def,1,x1,deriv=1)
     .             *scale_coord(gl_def,2,x2,deriv=1)
     .             *scale_coord(gl_def,3,x3,deriv=1)
        endif

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
        real(8) :: car(3),curv(3),vec(3),gp(ngparams),xlog(3)
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
          call getCurvilinearCoordinates(gl_def,i+1,j,k,igx,igy,igz
     .                                  ,ig,jg,kg,x1p,x2p,x3p)
          x1 = x1+x1p
          x2 = x2+x2p
          x3 = x3+x3p
          
          call getCurvilinearCoordinates(gl_def,i,j+1,k,igx,igy,igz
     .                                  ,ig,jg,kg,x1p,x2p,x3p)
          x1 = x1+x1p
          x2 = x2+x2p
          x3 = x3+x3p

          call getCurvilinearCoordinates(gl_def,i,j,k+1,igx,igy,igz
     .                                  ,ig,jg,kg,x1p,x2p,x3p)
          x1 = x1+x1p
          x2 = x2+x2p
          x3 = x3+x3p

          call getCurvilinearCoordinates(gl_def,i+1,j+1,k,igx,igy,igz
     .                                  ,ig,jg,kg,x1p,x2p,x3p)
          x1 = x1+x1p
          x2 = x2+x2p
          x3 = x3+x3p

          call getCurvilinearCoordinates(gl_def,i+1,j,k+1,igx,igy,igz
     .                                  ,ig,jg,kg,x1p,x2p,x3p)
          x1 = x1+x1p
          x2 = x2+x2p
          x3 = x3+x3p

          call getCurvilinearCoordinates(gl_def,i,j+1,k+1,igx,igy,igz
     .                                  ,ig,jg,kg,x1p,x2p,x3p)
          x1 = x1+x1p
          x2 = x2+x2p
          x3 = x3+x3p

          call getCurvilinearCoordinates(gl_def,i+1,j+1,k+1,igx,igy,igz
     .                                  ,ig,jg,kg,x1p,x2p,x3p)
          x1 = x1+x1p
          x2 = x2+x2p
          x3 = x3+x3p
          
          x1 = 0.125*x1
          x2 = 0.125*x2
          x3 = 0.125*x3
        endif

        if (coords /= 'lpk') then
          curv(1) = scale_coord(gl_def,1,x1)
          curv(2) = scale_coord(gl_def,2,x2)
          curv(3) = scale_coord(gl_def,3,x3)
        else
          curv = (/ x1, x2, x3 /)
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
          if (curv(1) == 0d0) curv(1) = SP_flsv

          select case (comp)
            case (1)
              vec = (/ cos(curv(2)),sin(curv(2)),0d0 /)
            case (2)
              vec = (/-sin(curv(2)),cos(curv(2)),0d0 /)/curv(1)
            case (3)
              vec = (/ 0d0,0d0,1d0 /)
          end select
        case ('hel')
          mm = gp(1)
          kk = gp(2)
          phi = (curv(2)-kk*curv(3))/mm

          if (curv(1) == 0d0) curv(1) = mm*SP_flsv

          select case (comp)
            case (1)
              vec = (/ cos(phi),sin(phi),0d0 /)
            case (2)
              vec = (/-mm*sin(phi)/curv(1),mm*cos(phi)/curv(1),kk /)
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
          if (curv(1) == 0d0) curv(1) = SP_flsv

          major_r = gp(1)
          aa      = gp(2)
          bb      = gp(3)

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

          if (curv(1) == 0d0) curv(1) = SP_flsv

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
        
        if (coords /= 'lpk') then
          xlog = (/x1,x2,x3/)
          vec = vec/scale_coord(gl_def,comp,xlog(comp),deriv=1)
        endif

      end function covariantVector

c$$$c     cov_vector_anal
c$$$c     #################################################################
c$$$      function cov_vector_anal(gl_def,curv) result(vec)
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
c$$$        real(8) :: curv(3),vec(3,3)
c$$$
c$$$        type(grid_mg_def),pointer :: gl_def
c$$$
c$$$c     Local variables
c$$$
c$$$        real(8) :: aa,bb,major_r,mm,kk,eps,phi,jac,r
c$$$     .            ,x1,x2,x3
c$$$        real(8) :: gp(ngparams)
c$$$
c$$$        real(8) :: xmax,ymax,zmax,xmin,ymin,zmin
c$$$! Daniele, 11/08/2017. Local variables needed for shaped toroidal coordinates
c$$$        integer :: ii
c$$$        real(8) :: ai(10),bi(10)
c$$$        real(8) :: ap,bp
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
c$$$          vec(1,:) = (/ 1d0,0d0,0d0 /)
c$$$          vec(2,:) = (/ 0d0,1d0,0d0 /)
c$$$          vec(3,:) = (/ 0d0,0d0,1d0 /)
c$$$        case ('cyl')
c$$$          if (curv(1) == 0d0) then
c$$$            vec = 0d0
c$$$            return
c$$$          endif
c$$$
c$$$          vec(1,:) = (/ cos(curv(2)),sin(curv(2)),0d0 /)
c$$$          vec(2,:) = (/-sin(curv(2)),cos(curv(2)),0d0 /)/curv(1)
c$$$          vec(3,:) = (/ 0d0,0d0,1d0 /)
c$$$        case ('hel')
c$$$          if (curv(1) == 0d0) then
c$$$            vec = 0d0
c$$$            return
c$$$          endif
c$$$
c$$$          mm = gp(1)
c$$$          kk = gp(2)
c$$$          phi = (curv(2)-aa*curv(3))/mm
c$$$
c$$$          vec(1,:) = (/ cos(phi),sin(phi),0d0 /)
c$$$          vec(2,:) = (/-mm*sin(phi)/curv(1),mm*cos(phi)/curv(1),kk /)
c$$$          vec(3,:) = (/ 0d0,0d0,1d0 /)
c$$$cc        case ('hl2')
c$$$cc          curv = (/ sign(sqrt(2*abs(x1)),x1),x2,x3 /)
c$$$cc
c$$$cc          mm = gp(1)
c$$$cc          kk = gp(2)
c$$$cc          aa = kk/mm
c$$$cc          phi = (curv(2)-aa*curv(3))
c$$$cc
c$$$cc          vec(1,:) = (/ cos(phi),sin(phi),0d0 /)*curv(1)
c$$$cc          vec(2,:) = (/-sin(phi)/curv(1),cos(phi)/curv(1),aa /)
c$$$cc          vec(3,:) = (/ 0d0,0d0,1d0 /)
c$$$        case ('tor')
c$$$          if (curv(1) == 0d0) then
c$$$            vec = 0d0
c$$$            return
c$$$          endif
c$$$
c$$$          major_r = gp(1)
c$$$          aa      = gp(2)
c$$$          bb      = gp(3)
c$$$
c$$$          vec(1,:) = (/  cos(curv(2))*cos(curv(3))/aa
c$$$     .                 ,-cos(curv(2))*sin(curv(3))/aa
c$$$     .                 , sin(curv(2))/bb             /)
c$$$          vec(2,:) = (/-sin(curv(2))*cos(curv(3))/aa
c$$$     .                 ,sin(curv(2))*sin(curv(3))/aa
c$$$     .                 ,cos(curv(2))/bb             /)/curv(1)
c$$$          vec(3,:) = (/ -sin(curv(3)),-cos(curv(3)),0d0 /)
c$$$     .                /(major_r + aa*curv(1)*cos(curv(2)))
c$$$        case ('sha')
c$$$! Daniele, 11/08/2017
c$$$! Shaped toroidal coordinates system
c$$$! x = (R0 + r*sum_i(ai*cos(i*theta))*cos(-phi)
c$$$! y = (R0 + r*sum_i(ai*cos(i*theta))*sin(-phi)
c$$$! z = r*sum_i(bi*sin(i*theta)
c$$$
c$$$          if (curv(1) == 0d0) then
c$$$            vec = 0d0
c$$$            return
c$$$          endif
c$$$
c$$$          major_r = gp(1)  !Major radius
c$$$
c$$$          ai(1) = gp(2)  !Horizontal minor radius
c$$$          bi(1) = gp(3)  !Vertical   minor radius
c$$$
c$$$          do ii=2,10
c$$$             ai(ii) = gp(7+ii) ! Higher order cosine coefficients are stored in gp(9) to gp(17)
c$$$             bi(ii) = gp(16+ii) ! Higher order sine coefficients are stored in gp(18) to gp(26)
c$$$          enddo
c$$$
c$$$          aa = 0d0
c$$$          bb = 0d0
c$$$          ap = 0d0
c$$$          bp = 0d0
c$$$          do ii=1,10
c$$$             aa = aa + ai(ii)*cos(ii*curv(2)) ! A(theta)=sum_i(ai*cos(i*theta))
c$$$             bb = bb + bi(ii)*sin(ii*curv(2)) ! B(theta)=sum_i(bi*sin(i*theta))
c$$$             bp = bp + bi(ii)*ii*cos(ii*curv(2)) ! B'(theta)
c$$$             ap = ap - ai(ii)*ii*sin(ii*curv(2)) ! A'(theta)
c$$$          enddo
c$$$
c$$$          vec(1,:) = (/  bp*cos(curv(3))
c$$$     .                 ,-bp*sin(curv(3))
c$$$     .                 ,-ap             /)/(aa*bp-bb*ap)
c$$$          vec(2,:) = (/ -bb*cos(curv(3))
c$$$     .                 , bb*sin(curv(3))
c$$$     .                 , aa             /)/curv(1)/(aa*bp-bb*ap)
c$$$          vec(3,:) = (/ -sin(curv(3)),-cos(curv(3)),0d0 /)
c$$$     .                /(major_r + aa*curv(1))
c$$$
c$$$        case ('sin')
c$$$
c$$$          eps = gp(1)
c$$$
c$$$          x1 = curv(1)
c$$$          x2 = curv(2)
c$$$          x3 = curv(3)
c$$$          vec(1,:) = (/(xmax*(ymax + 2*eps*Pi*Cos((2*Pi*x2)/ymax)*
c$$$     -                                        Sin((2*Pi*x1)/xmax)))/
c$$$     -                 (xmax*ymax + eps*Pi*(xmax - ymax)
c$$$     -                         *Sin(Pi*((2*x1)/xmax - (2*x2)/ymax))
c$$$     -                            + eps*Pi*(xmax + ymax)
c$$$     .                         *Sin(2*Pi*(x1/xmax + x2/ymax)))
c$$$     .              , (-2*eps*Pi*xmax*Cos((2*Pi*x2)/ymax)
c$$$     .                               *Sin((2*Pi*x1)/xmax))/
c$$$     -                (xmax*ymax + eps*Pi*(xmax - ymax)
c$$$     -                         *Sin(Pi*((2*x1)/xmax - (2*x2)/ymax))  
c$$$     -                           + eps*Pi*(xmax + ymax)
c$$$     .                         *Sin(2*Pi*(x1/xmax + x2/ymax)))
c$$$     .              , 0d0 /)
c$$$          vec(2,:) = (/(-2*eps*Pi*ymax*Cos((2*Pi*x1)/xmax)
c$$$     .                                *Sin((2*Pi*x2)/ymax))/
c$$$     -                 (xmax*ymax + eps*Pi*(xmax - ymax)
c$$$     -                         *Sin(Pi*((2*x1)/xmax - (2*x2)/ymax))
c$$$     -                            + eps*Pi*(xmax + ymax)
c$$$     .                         *Sin(2*Pi*(x1/xmax + x2/ymax)))
c$$$     .              , (ymax*(xmax + 2*eps*Pi*Cos((2*Pi*x1)/xmax)
c$$$     .                                      *Sin((2*Pi*x2)/ymax)))/
c$$$     -                (xmax*ymax + eps*Pi*(xmax - ymax)
c$$$     -                         *Sin(Pi*((2*x1)/xmax - (2*x2)/ymax))
c$$$     -                           + eps*Pi*(xmax + ymax)
c$$$     .                         *Sin(2*Pi*(x1/xmax + x2/ymax)))
c$$$     .              ,0d0 /)
c$$$          vec(3,:) = (/ 0d0,0d0,1d0 /)
c$$$        case default
c$$$          write (*,*) 'Grid not implemented in covariantVector'
c$$$          write (*,*) 'Aborting...'
c$$$          stop
c$$$        end select
c$$$
c$$$        if (coords /= 'lpk') then
c$$$          vec(1,:) = vec(1,:)*scale_coord(gl_def,2,curv(2),deriv=.true.)
c$$$     .                       *scale_coord(gl_def,3,curv(3),deriv=.true.)
c$$$          vec(2,:) = vec(2,:)*scale_coord(gl_def,1,curv(1),deriv=.true.)
c$$$     .                       *scale_coord(gl_def,3,curv(3),deriv=.true.)
c$$$          vec(3,:) = vec(3,:)*scale_coord(gl_def,1,curv(1),deriv=.true.)
c$$$     .                       *scale_coord(gl_def,2,curv(2),deriv=.true.)
c$$$        endif
c$$$
c$$$      end function cov_vector_anal

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
        real(8) :: car(3),curv(3),jac,xlog(3)
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
          call getCurvilinearCoordinates(gl_def,i+1,j,k,igx,igy,igz
     .                                  ,ig,jg,kg,x1p,x2p,x3p)
          x1 = x1+x1p
          x2 = x2+x2p
          x3 = x3+x3p
          
          call getCurvilinearCoordinates(gl_def,i,j+1,k,igx,igy,igz
     .                                  ,ig,jg,kg,x1p,x2p,x3p)
          x1 = x1+x1p
          x2 = x2+x2p
          x3 = x3+x3p

          call getCurvilinearCoordinates(gl_def,i,j,k+1,igx,igy,igz
     .                                  ,ig,jg,kg,x1p,x2p,x3p)
          x1 = x1+x1p
          x2 = x2+x2p
          x3 = x3+x3p

          call getCurvilinearCoordinates(gl_def,i+1,j+1,k,igx,igy,igz
     .                                  ,ig,jg,kg,x1p,x2p,x3p)
          x1 = x1+x1p
          x2 = x2+x2p
          x3 = x3+x3p

          call getCurvilinearCoordinates(gl_def,i+1,j,k+1,igx,igy,igz
     .                                  ,ig,jg,kg,x1p,x2p,x3p)
          x1 = x1+x1p
          x2 = x2+x2p
          x3 = x3+x3p

          call getCurvilinearCoordinates(gl_def,i,j+1,k+1,igx,igy,igz
     .                                  ,ig,jg,kg,x1p,x2p,x3p)
          x1 = x1+x1p
          x2 = x2+x2p
          x3 = x3+x3p

          call getCurvilinearCoordinates(gl_def,i+1,j+1,k+1,igx,igy,igz
     .                                  ,ig,jg,kg,x1p,x2p,x3p)
          x1 = x1+x1p
          x2 = x2+x2p
          x3 = x3+x3p
          
          x1 = 0.125*x1
          x2 = 0.125*x2
          x3 = 0.125*x3
          
c$$$          call getCurvilinearCoordinates(gl_def,i+1,j+1,k+1,igx,igy,igz
c$$$     .                                  ,ig,jg,kg,x1p,x2p,x3p)
c$$$          x1 = 0.5*(x1+x1p)
c$$$          x2 = 0.5*(x2+x2p)
c$$$          x3 = 0.5*(x3+x3p)
        endif

        if (coords /= 'lpk') then
          curv(1) = scale_coord(gl_def,1,x1)
          curv(2) = scale_coord(gl_def,2,x2)
          curv(3) = scale_coord(gl_def,3,x3)
        else
          curv = (/ x1, x2, x3 /)
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
          if (curv(1) == 0d0) curv(1) = SP_flsv

          select case (comp)
            case (1)
              vec = (/ cos(curv(2)),sin(curv(2)),0d0 /)/curv(1)
            case (2)
              vec = (/-sin(curv(2)),cos(curv(2)),0d0 /)
            case (3)
              vec = (/ 0d0,0d0,1d0 /)/curv(1)
          end select
        case ('hel')
          mm = gp(1)
          kk = gp(2)
          phi = (curv(2)-kk*curv(3))/mm

          select case (comp)
            case (1)
              if (curv(1) == 0d0) curv(1) = mm*SP_flsv
              vec = (/ cos(phi),sin(phi),0d0 /)*mm/curv(1)
            case (2)
              vec = (/-sin(phi),cos(phi),0d0 /)
            case (3)
              vec=(/kk*sin(phi)*curv(1)/mm,-kk*cos(phi)*curv(1)/mm,1d0/)
              if (curv(1) == 0d0) curv(1) = mm*SP_flsv
              vec = vec*mm/curv(1)
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
          if (curv(1) == 0d0) curv(1) = SP_flsv

          major_r = gp(1)
          aa      = gp(2)
          bb      = gp(3)

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

          if (curv(1) == 0d0) curv(1) = SP_flsv

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

        if (coords /= 'lpk') then
          xlog = (/x1,x2,x3/)

          vec = vec*scale_coord(gl_def,comp,xlog(comp),deriv=1)
     .            /(scale_coord(gl_def,1,xlog(1),deriv=1)
     .             *scale_coord(gl_def,2,xlog(2),deriv=1)
     .             *scale_coord(gl_def,3,xlog(3),deriv=1))
        endif

      end function contravariantVector

c$$$c     cnv_vector_anal
c$$$c     #################################################################
c$$$      function cnv_vector_anal(gl_def,curv) result(vec)
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
c$$$        real(8) :: curv(3),vec(3,3)
c$$$
c$$$        type(grid_mg_def),pointer :: gl_def
c$$$
c$$$c     Local variables
c$$$
c$$$        real(8) :: aa,bb,major_r,mm,kk,eps,phi,jac,r
c$$$     .            ,x1,x2,x3
c$$$        real(8) :: gp(ngparams)
c$$$
c$$$        real(8) :: xmax,ymax,zmax,xmin,ymin,zmin
c$$$! Daniele, 11/08/2017. Local variables needed for shaped toroidal coordinates
c$$$        integer :: ii
c$$$        real(8) :: ai(10),bi(10)
c$$$        real(8) :: ap,bp
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
c$$$          vec(1,:) = (/ 1d0,0d0,0d0 /)
c$$$          vec(2,:) = (/ 0d0,1d0,0d0 /)
c$$$          vec(3,:) = (/ 0d0,0d0,1d0 /)
c$$$        case ('cyl')
c$$$          if (curv(1) == 0d0) then
c$$$            vec = 0d0
c$$$            return
c$$$          endif
c$$$
c$$$          vec(1,:) = (/ cos(curv(2)),sin(curv(2)),0d0 /)/curv(1)
c$$$          vec(2,:) = (/-sin(curv(2)),cos(curv(2)),0d0 /)
c$$$          vec(3,:) = (/ 0d0,0d0,1d0 /)/curv(1)
c$$$        case ('hel')
c$$$          if (curv(1) == 0d0) then
c$$$            vec = 0d0
c$$$            return
c$$$          endif
c$$$
c$$$          mm = gp(1)
c$$$          kk = gp(2)
c$$$          phi = (curv(2)-aa*curv(3))/mm
c$$$
c$$$          vec(1,:) = (/ cos(phi),sin(phi),0d0 /)*mm/curv(1)
c$$$          vec(2,:) = (/-sin(phi),cos(phi),0d0 /)
c$$$          vec(3,:) = (/ kk*sin(phi),-kk*cos(phi),mm/curv(1) /)
c$$$cc        case ('hl2')
c$$$cc          curv = (/ sign(sqrt(2*abs(x1)),x1),x2,x3 /)
c$$$cc
c$$$cc          mm = gp(1)
c$$$cc          kk = gp(2)
c$$$cc          aa = kk/mm
c$$$cc          phi = (curv(2)-aa*curv(3))
c$$$cc
c$$$cc          vec(1,:) = (/ cos(phi),sin(phi),0d0 /)/curv(1)
c$$$cc          vec(2,:) = (/-sin(phi),cos(phi),0d0 /)*curv(1)
c$$$cc          vec(3,:) = (/ aa*sin(phi)*curv(1),-aa*cos(phi)*curv(1),1d0 /)
c$$$        case ('tor')
c$$$          if (curv(1) == 0d0) then
c$$$            vec = 0d0
c$$$            return
c$$$          endif
c$$$
c$$$          major_r = gp(1)
c$$$          aa      = gp(2)
c$$$          bb      = gp(3)
c$$$
c$$$          vec(1,:) = (/  cos(curv(2))*cos(curv(3))/bb
c$$$     .                ,-cos(curv(2))*sin(curv(3))/bb
c$$$     .                , sin(curv(2))/aa             /)
c$$$     .              /curv(1)/(major_r + aa*curv(1)*cos(curv(2)))
c$$$          vec(2,:) = (/-sin(curv(2))*cos(curv(3))/bb
c$$$     .                ,sin(curv(2))*sin(curv(3))/bb
c$$$     .                ,cos(curv(2))/aa              /)
c$$$     .              /(major_r + aa*curv(1)*cos(curv(2)))
c$$$          vec(3,:) = (/ -sin(curv(3)),-cos(curv(3)),0d0 /)/curv(1)/aa/bb
c$$$
c$$$        case ('sha')
c$$$! Daniele, 11/08/2017
c$$$! Shaped toroidal coordinates system
c$$$! x = (R0 + r*sum_i(ai*cos(i*theta))*cos(-phi)
c$$$! y = (R0 + r*sum_i(ai*cos(i*theta))*sin(-phi)
c$$$! z = r*sum_i(bi*sin(i*theta)
c$$$
c$$$          if (curv(1) == 0d0) then
c$$$            vec = 0d0
c$$$            return
c$$$          endif
c$$$
c$$$          major_r = gp(1)  !Major radius
c$$$
c$$$          ai(1) = gp(2)  !Horizontal minor radius
c$$$          bi(1) = gp(3)  !Vertical   minor radius
c$$$
c$$$          do ii=2,10
c$$$             ai(ii) = gp(7+ii) ! Higher order cosine coefficients are stored in gp(9) to gp(17)
c$$$             bi(ii) = gp(16+ii) ! Higher order sine coefficients are stored in gp(18) to gp(26)
c$$$          enddo
c$$$
c$$$          aa = 0d0
c$$$          bb = 0d0
c$$$          ap = 0d0
c$$$          bp = 0d0
c$$$          do ii=1,10
c$$$             aa = aa + ai(ii)*cos(ii*curv(2)) ! A(theta)=sum_i(ai*cos(i*theta))
c$$$             bb = bb + bi(ii)*sin(ii*curv(2)) ! B(theta)=sum_i(bi*sin(i*theta))
c$$$             bp = bp + bi(ii)*ii*cos(ii*curv(2)) ! B'(theta)
c$$$             ap = ap - ai(ii)*ii*sin(ii*curv(2)) ! A'(theta)
c$$$          enddo
c$$$
c$$$          vec(1,:) = (/ aa*cos(curv(3))
c$$$     .                ,-aa*sin(curv(3))
c$$$     .                , bb             /)
c$$$     .              /curv(1)/(major_r + aa*curv(1))/(aa*bp-bb*ap)
c$$$          vec(2,:) = (/ ap*cos(curv(3))
c$$$     .                ,-ap*sin(curv(3))
c$$$     .                , bp             /)
c$$$     .              /(major_r + aa*curv(1))/(aa*bp-bb*ap)
c$$$          vec(3,:) = (/ -sin(curv(3)),-cos(curv(3)),0d0 /)
c$$$     .              /curv(1)/(aa*bp-bb*ap)
c$$$
c$$$        case ('sin')
c$$$
c$$$          eps = gp(1)
c$$$
c$$$          x1 = curv(1)
c$$$          x2 = curv(2)
c$$$          x3 = curv(3)
c$$$
c$$$          vec(1,:) =(/ (ymax*(xmax + 2*eps*Pi*Cos((2*Pi*x1)/xmax)*
c$$$     -                                        Sin((2*Pi*x2)/ymax)))/
c$$$     -                 (xmax*ymax + eps*Pi*(xmax - ymax)
c$$$     -                         *Sin(Pi*((2*x1)/xmax - (2*x2)/ymax))
c$$$     -                            + eps*Pi*(xmax + ymax)
c$$$     .                         *Sin(2*Pi*(x1/xmax + x2/ymax)))
c$$$     .              , (2*eps*Pi*ymax*Cos((2*Pi*x1)/xmax)
c$$$     .                              *Sin((2*Pi*x2)/ymax))/
c$$$     -                (xmax*ymax + eps*Pi*(xmax - ymax)
c$$$     -                         *Sin(Pi*((2*x1)/xmax - (2*x2)/ymax))
c$$$     -                           + eps*Pi*(xmax + ymax)
c$$$     .                         *Sin(2*Pi*(x1/xmax + x2/ymax)))
c$$$     .              ,0d0 /)
c$$$          vec(2,:) =(/ (2*eps*Pi*xmax*Cos((2*Pi*x2)/ymax)
c$$$     .                               *Sin((2*Pi*x1)/xmax))/
c$$$     -                 (xmax*ymax + eps*Pi*(xmax - ymax)
c$$$     -                         *Sin(Pi*((2*x1)/xmax - (2*x2)/ymax))
c$$$     -                            + eps*Pi*(xmax + ymax)
c$$$     .                         *Sin(2*Pi*(x1/xmax + x2/ymax)))
c$$$     .              ,(xmax*(ymax + 2*eps*Pi*Cos((2*Pi*x2)/ymax)
c$$$     .                                     *Sin((2*Pi*x1)/xmax)))/
c$$$     -               (xmax*ymax + eps*Pi*(xmax - ymax)
c$$$     -                         *Sin(Pi*((2*x1)/xmax - (2*x2)/ymax))
c$$$     -                          + eps*Pi*(xmax + ymax)
c$$$     .                         *Sin(2*Pi*(x1/xmax + x2/ymax)))
c$$$     .              ,0d0 /)
c$$$          vec(3,:) =(/ 0d0
c$$$     .                ,0d0
c$$$     .                ,(xmax*ymax)/(xmax*ymax
c$$$     .                  + eps*Pi*(xmax - ymax)
c$$$     -                       *Sin(Pi*((2*x1)/xmax - (2*x2)/ymax))
c$$$     -                  + eps*Pi*(xmax + ymax)
c$$$     .                       *Sin(2*Pi*(x1/xmax + x2/ymax))) /)
c$$$        case default
c$$$          write (*,*) 'Grid not implemented in covariantVector'
c$$$          write (*,*) 'Aborting...'
c$$$          stop
c$$$        end select
c$$$
c$$$        if (coords /= 'lpk') then
c$$$          vec(1,:) = vec(1,:)*scale_coord(gl_def,1,curv(1),deriv=1)
c$$$          vec(2,:) = vec(2,:)*scale_coord(gl_def,2,curv(2),deriv=1)
c$$$          vec(3,:) = vec(3,:)*scale_coord(gl_def,3,curv(3),deriv=1)
c$$$        endif
c$$$
c$$$      end function cnv_vector_anal

c     hessian22
c     #################################################################
      function hessian22(gl_def,l,i,j,k,igx,igy,igz) result (tensor)

c     -----------------------------------------------------------------
c     Calculates hessian elements of curvilinear coordinate system in
c     the covariant basis, i.e., 
c            hessian[l](n,m) = J^2 <cnv(n)|grad(cov[l])|cnv(m)>
c     -----------------------------------------------------------------

        implicit none

c     Input variables

        integer :: l,i,j,k,igx,igy,igz

        type(grid_mg_def),pointer :: gl_def

        real(8) :: tensor(3,3)
        
c     Local variables

        integer    :: ig,jg,kg
        real(8)    :: car(3),curv(3),vec(9)
        real(8)    :: x1,x2,x3,aa,bb,major_r,mm,kk,eps,phi,gp(ngparams)

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

        call getCurvilinearCoordinates(gl_def,i,j,k,igx,igy,igz
     .                                ,ig,jg,kg,x1,x2,x3)

        if (coords /= 'lpk') then
          curv(1) = scale_coord(gl_def,1,x1)
          curv(2) = scale_coord(gl_def,2,x2)
          curv(3) = scale_coord(gl_def,3,x3)
        else
          curv = (/ x1, x2, x3 /)
        endif

        select case (coords)
        case ('car','lpk')
          vec = (/ 0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0 /)
        case ('cyl')
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
          mm = gp(1)
          kk = gp(2)
          phi = (curv(2)-kk*curv(3))/mm
          select case (l)
            case (1)
              vec(1) = 0d0
              vec(2) = 0d0
              vec(3) = 0d0
              vec(4) = vec(2)
              vec(5) = curv(1)/(mm*mm)
              vec(6) = -kk*curv(1)/(mm*mm)
              vec(7) = vec(3)
              vec(8) = vec(6)
              vec(9) = kk*kk*curv(1)/(mm*mm)
            case (2)
              vec(1) = 0d0
              vec(2) = -1./curv(1)
              vec(3) = kk/curv(1)
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

        integer :: n,m,ig,jg,kg
        real(8) :: curv(3),xlog(3)

c     Begin program

        tensor(1,:,:) = hessian22(gl_def,1,i,j,k,igx,igy,igz)
        tensor(2,:,:) = hessian22(gl_def,2,i,j,k,igx,igy,igz)
        tensor(3,:,:) = hessian22(gl_def,3,i,j,k,igx,igy,igz)

        if (coords /= 'lpk') then
          call getCurvilinearCoordinates(gl_def,i,j,k,igx,igy,igz
     .                                ,ig,jg,kg,xlog(1),xlog(2),xlog(3))

          curv(1) = scale_coord(gl_def,1,xlog(1))
          curv(2) = scale_coord(gl_def,2,xlog(2))
          curv(3) = scale_coord(gl_def,3,xlog(3))

          do l=1,3
            do n=1,3
              do m=1,3
                tensor(l,m,n)=tensor(l,m,n)
     .                 *scale_coord(gl_def,m,xlog(m),deriv=1)
     .                 *scale_coord(gl_def,n,xlog(n),deriv=1)
     .                 /scale_coord(gl_def,l,xlog(l),deriv=1)
                if (m==n.and.m==l) then
                  tensor(l,m,n)=tensor(l,m,n) 
     .                +scale_coord(gl_def,m,xlog(m),deriv=2)
     .                /scale_coord(gl_def,l,xlog(l),deriv=1)
                endif
              enddo
            enddo
          enddo
        endif

c     End program

      end function christ_2knd

c$$$c     hessian_cnv
c$$$c     #################################################################
c$$$      function hessian_cnv(gl_def,k,x1,x2,x3) result (tensor)
c$$$
c$$$c     -----------------------------------------------------------------
c$$$c     Calculates elements of tensor grad(cnv) in a mixed coordinate 
c$$$c     system:
c$$$c              hessian_cnv[k](i,j) = J^2 <cnv(i)|grad(cnv[k])|cov(j)>
c$$$c     -----------------------------------------------------------------
c$$$
c$$$        implicit none
c$$$
c$$$c     Input variables
c$$$
c$$$        integer    :: k
c$$$        real(8)    :: x1,x2,x3,tensor(3,3)
c$$$
c$$$        type(grid_mg_def),pointer :: gl_def
c$$$
c$$$c     Local variables
c$$$
c$$$        real(8)    :: vec(9),car(3),aa,bb,major_r,mm,kk,eps,gp(ngparams)
c$$$
c$$$c     Begin program
c$$$
c$$$        gp = gl_def%params
c$$$
c$$$        select case (coords)
c$$$        case ('car','lpk')
c$$$          vec = (/ 0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0 /)
c$$$        case ('cyl')
c$$$          select case (k)
c$$$            case (1)
c$$$              vec(1) = -1./x1
c$$$              vec(2) = 0d0
c$$$              vec(3) = 0d0
c$$$              vec(4) = 0d0
c$$$              vec(5) = 1./x1
c$$$              vec(6) = 0d0
c$$$              vec(7) = 0d0
c$$$              vec(8) = 0d0
c$$$              vec(9) = 0d0
c$$$            case (2)
c$$$              vec(1) = 0d0
c$$$              vec(2) = 0d0
c$$$              vec(3) = 0d0
c$$$              vec(4) = -x1
c$$$              vec(5) = 0d0
c$$$              vec(6) = 0d0
c$$$              vec(7) = 0d0
c$$$              vec(8) = 0d0
c$$$              vec(9) = 0d0
c$$$            case (3)
c$$$              vec(1) = 0d0
c$$$              vec(2) = 0d0
c$$$              vec(3) = -1./x1
c$$$              vec(4) = 0d0
c$$$              vec(5) = 0d0
c$$$              vec(6) = 0d0
c$$$              vec(7) = 0d0
c$$$              vec(8) = 0d0
c$$$              vec(9) = 0d0
c$$$          end select
c$$$        case ('hel')
c$$$          mm = gp(1)
c$$$          kk = gp(2)
c$$$          aa = kk/mm 
c$$$          select case (k)
c$$$          case (1)
c$$$            vec(1) = -1./x1
c$$$            vec(2) = 0d0
c$$$            vec(3) = 0d0
c$$$            vec(4) = 0d0
c$$$            vec(5) = 1./x1
c$$$            vec(6) = 0d0
c$$$            vec(7) = 0d0
c$$$            vec(8) = -kk/x1
c$$$            vec(9) = 0d0
c$$$          case (2)
c$$$            vec(1) = 0d0
c$$$            vec(2) = 0d0
c$$$            vec(3) = 0d0
c$$$            vec(4) = -x1/(mm*mm)
c$$$            vec(5) = 0d0
c$$$            vec(6) = 0d0
c$$$            vec(7) = kk*x1/(mm*mm)
c$$$            vec(8) = 0d0
c$$$            vec(9) = 0d0
c$$$          case (3)
c$$$            vec(1) = 0d0
c$$$            vec(2) = -kk/x1
c$$$            vec(3) = -1./x1
c$$$            vec(4) = kk*x1/(mm*mm)
c$$$            vec(5) = 0d0
c$$$            vec(6) = 0d0
c$$$            vec(7) = -kk*kk*x1/(mm*mm)
c$$$            vec(8) = 0d0
c$$$            vec(9) = 0d0
c$$$          end select
c$$$        case ('tor')
c$$$          major_r = gp(1)
c$$$          select case (k)
c$$$            case (1)
c$$$              vec(1) = -(major_r + 2*x1*sin(x2))
c$$$     .                  /x1/(major_r + x1*sin(x2))
c$$$              vec(2) = 0d0
c$$$              vec(3) = 0d0
c$$$              vec(4) = -x1*cos(x2)/(major_r + x1*sin(x2))
c$$$              vec(5) = 1./x1
c$$$              vec(6) = 0d0
c$$$              vec(7) = 0d0
c$$$              vec(8) = 0d0
c$$$              vec(9) = sin(x2)/(major_r + x1*sin(x2))
c$$$            case (2)
c$$$              vec(1) = 0d0
c$$$              vec(2) = -sin(x2)/(major_r + x1*Sin(x2))
c$$$              vec(3) = 0d0
c$$$              vec(4) = -x1
c$$$              vec(5) = -x1*cos(x2)/(major_r + x1*Sin(x2))
c$$$              vec(6) = 0d0
c$$$              vec(7) = 0d0
c$$$              vec(8) = 0d0
c$$$              vec(9) =  x1*cos(x2)/(major_r + x1*Sin(x2))
c$$$            case (3)
c$$$              vec(1) =  0d0
c$$$              vec(2) =  0d0
c$$$              vec(3) =  -1./x1
c$$$              vec(4) = 0d0
c$$$              vec(5) = 0d0
c$$$              vec(6) = 0d0
c$$$              vec(7) = -sin(x2)*(major_r + x1*Sin(x2))
c$$$              vec(8) = -cos(x2)*(major_r + x1*Sin(x2))/x1
c$$$              vec(9) = 0d0
c$$$          end select
c$$$
c$$$        case default
c$$$          write (*,*) 'Grid not implemented in hessian_cnv'
c$$$          write (*,*) 'Aborting...'
c$$$          stop
c$$$        end select
c$$$
c$$$        tensor = transpose(reshape(vec, (/3,3/)))
c$$$
c$$$      end function hessian_cnv

c$$$c     g_sub
c$$$c     #################################################################
c$$$      function g_sub(gl_def,i,j,k,igx,igy,igz,vertex) result (tensor)
c$$$
c$$$c     -----------------------------------------------------------------
c$$$c     Calculates contravariant metric tensor of curvilinear coordinate 
c$$$c     system
c$$$c     -----------------------------------------------------------------
c$$$
c$$$        implicit none
c$$$
c$$$c     Input variables
c$$$
c$$$        integer :: i,j,k,igx,igy,igz
c$$$        real(8) :: tensor(3,3)
c$$$        logical,optional :: vertex
c$$$
c$$$        type(grid_mg_def),pointer :: gl_def
c$$$
c$$$c     Local variables
c$$$
c$$$        integer :: ig,jg,kg
c$$$        real(8) :: vec(9),jac,car(3),curv(3)
c$$$        real(8) :: x1,x2,x3,aa,bb,major_r,mm,kk,eps,gp(ngparams)
c$$$     .            ,x1p,x2p,x3p,i_x1
c$$$        logical :: vrtx
c$$$
c$$$        real(8) :: xmax,ymax,zmax,xmin,ymin,zmin
c$$$! Daniele, 11/08/2017. Local variables needed for shaped toroidal coordinates
c$$$        integer :: ii
c$$$        real(8) :: ai(10),bi(10)
c$$$        real(8) :: ap,bp
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
c$$$        if (PRESENT(vertex)) then
c$$$          vrtx = vertex
c$$$        else
c$$$          vrtx = .false.
c$$$        endif
c$$$
c$$$        gp = gl_def%params
c$$$
c$$$        call getCurvilinearCoordinates(gl_def,i,j,k,igx,igy,igz,ig,jg,kg
c$$$     .                                ,x1,x2,x3)
c$$$
c$$$        if (vrtx) then
c$$$          call getCurvilinearCoordinates(gl_def,i+1,j+1,k+1,igx,igy,igz
c$$$     .                                  ,ig,jg,kg,x1p,x2p,x3p)
c$$$          x1 = 0.5*(x1+x1p)
c$$$          x2 = 0.5*(x2+x2p)
c$$$          x3 = 0.5*(x3+x3p)
c$$$        endif
c$$$
c$$$        select case (coords)
c$$$        case ('car','lpk')
c$$$          vec = (/ 1d0, 0d0, 0d0
c$$$     .            ,0d0, 1d0, 0d0
c$$$     .            ,0d0, 0d0, 1d0 /)
c$$$        case ('cyl')
c$$$          if (x1 == 0d0) then
c$$$            vec = 0d0
c$$$          else
c$$$            vec = (/ 1d0/x1, 0d0  , 0d0
c$$$     .              ,0d0   , x1   , 0d0
c$$$     .              ,0d0   , 0d0  , 1d0/x1 /)
c$$$          endif
c$$$        case ('hel')
c$$$          mm = gp(1)
c$$$          kk = gp(2)
c$$$          if (x1 == 0d0) then
c$$$            vec = 0d0
c$$$          else
c$$$            vec = (/ mm/x1, 0d0      ,    0d0
c$$$     .              ,0d0  , x1/mm    ,  -kk*x1/mm
c$$$     .              ,0d0  ,-kk*x1/mm , mm/x1 + kk*kk*x1/mm /)
c$$$          endif
c$$$        case ('hl2')
c$$$          curv = (/ abs(2*x1),x2,x3 /)
c$$$
c$$$          mm = gp(1)
c$$$          kk = gp(2)
c$$$          aa = kk/mm
c$$$          vec = (/ 1./curv(1),0d0, 0d0
c$$$     .            ,0d0, curv(1) ,-aa*curv(1)
c$$$     .            ,0d0,-aa*curv(1), 1 + aa**2*curv(1) /)
c$$$        case ('tor')
c$$$          if (x1 == 0d0) then
c$$$            vec = 0d0
c$$$          else
c$$$            major_r = gp(1)
c$$$            aa      = gp(2)
c$$$            bb      = gp(3)
c$$$
c$$$            vec(1) = (aa**2 + bb**2 + (aa**2 - bb**2)*Cos(2*x2))/
c$$$     -               (2.*aa*bb*x1*(major_r + aa*x1*Cos(x2)))
c$$$
c$$$            vec(2) = ((-aa**2 + bb**2)*Cos(x2)*Sin(x2))
c$$$     .               /(aa*bb*(major_r + aa*x1*Cos(x2)))
c$$$
c$$$            vec(3) = 0d0
c$$$
c$$$            vec(4) = vec(2)
c$$$
c$$$            vec(5) = (x1*(aa**2 + bb**2 + (-aa**2 + bb**2)*Cos(2*x2)))/
c$$$     -               (2.*aa*bb*(major_r + aa*x1*Cos(x2)))
c$$$
c$$$            vec(6) = 0d0
c$$$
c$$$            vec(7) = vec(3)
c$$$            vec(8) = vec(6)
c$$$            vec(9) = (major_r + aa*x1*Cos(x2))/(aa*bb*x1)
c$$$          endif
c$$$        case ('sha')
c$$$! Daniele, 11/08/2017
c$$$! Shaped toroidal coordinates system
c$$$! x = (R0 + r*sum_i(ai*cos(i*theta))*cos(-phi)
c$$$! y = (R0 + r*sum_i(ai*cos(i*theta))*sin(-phi)
c$$$! z = r*sum_i(bi*sin(i*theta)
c$$$          if (x1 == 0d0) then
c$$$            vec = 0d0
c$$$          else
c$$$          major_r = gp(1)  !Major radius
c$$$
c$$$          ai(1) = gp(2)  !Horizontal minor radius
c$$$          bi(1) = gp(3)  !Vertical   minor radius
c$$$
c$$$          do ii=2,10
c$$$             ai(ii) = gp(7+ii) ! Higher order cosine coefficients are stored in gp(9) to gp(17)
c$$$             bi(ii) = gp(16+ii) ! Higher order sine coefficients are stored in gp(18) to gp(26)
c$$$          enddo
c$$$
c$$$          aa = 0d0
c$$$          bb = 0d0
c$$$          ap = 0d0
c$$$          bp = 0d0
c$$$          do ii=1,10
c$$$             aa = aa + ai(ii)*cos(ii*x2) ! A(theta)=sum_i(ai*cos(i*theta))
c$$$             bb = bb + bi(ii)*sin(ii*x2) ! B(theta)=sum_i(bi*sin(i*theta))
c$$$             bp = bp + bi(ii)*ii*cos(ii*x2) ! B'(theta)
c$$$             ap = ap - ai(ii)*ii*sin(ii*x2) ! A'(theta)
c$$$          enddo
c$$$
c$$$            vec(1) = (aa**2 + bb**2)/
c$$$     .               (x1*(major_r + aa*x1)*(aa*bp-bb*ap))
c$$$
c$$$            vec(2) = (aa*ap + bb*bp)/
c$$$     .               ((major_r + aa*x1)*(aa*bp-bb*ap))
c$$$
c$$$            vec(3) = 0d0
c$$$
c$$$            vec(4) = vec(2)
c$$$
c$$$            vec(5) = x1*(ap**2 + bp**2)/
c$$$     .               ((major_r + aa*x1)*(aa*bp-bb*ap))
c$$$
c$$$            vec(6) = 0d0
c$$$
c$$$            vec(7) = vec(3)
c$$$            vec(8) = vec(6)
c$$$            vec(9) = (major_r + aa*x1)/x1/(aa*bp-bb*ap)
c$$$
c$$$          endif
c$$$        case ('sin')
c$$$
c$$$          eps = gp(1)
c$$$
c$$$          vec(1) =
c$$$     .         (ymax*(2*eps**2*Pi**2 + xmax**2 + 
c$$$     -        2*eps**2*Pi**2*Cos((4*Pi*x1)/xmax) - 
c$$$     -        eps**2*Pi**2*Cos(Pi*((4*x1)/xmax - (4*x2)/ymax)) - 
c$$$     -        eps**2*Pi**2*Cos(4*Pi*(x1/xmax + x2/ymax)) - 
c$$$     -        2*eps**2*Pi**2*Cos((4*Pi*x2)/ymax) - 
c$$$     -        2*eps*Pi*xmax*Sin(Pi*((2*x1)/xmax - (2*x2)/ymax)) + 
c$$$     -        2*eps*Pi*xmax*Sin(2*Pi*(x1/xmax + x2/ymax))))/
c$$$     -    (xmax*(xmax*ymax + eps*Pi*(xmax - ymax)*
c$$$     -         Sin(Pi*((2*x1)/xmax - (2*x2)/ymax)) + 
c$$$     -        eps*Pi*(xmax + ymax)*Sin(2*Pi*(x1/xmax + x2/ymax))))
c$$$
c$$$          vec(2) =
c$$$     .         (eps*Pi*(eps*Pi*Cos(Pi*((4*x1)/xmax - (4*x2)/ymax)) - 
c$$$     -        eps*Pi*Cos(4*Pi*(x1/xmax + x2/ymax)) + 
c$$$     -        xmax*Sin(Pi*((2*x1)/xmax - (2*x2)/ymax)) - 
c$$$     -        ymax*Sin(Pi*((2*x1)/xmax - (2*x2)/ymax)) + 
c$$$     -        xmax*Sin(2*Pi*(x1/xmax + x2/ymax)) + 
c$$$     -        ymax*Sin(2*Pi*(x1/xmax + x2/ymax))))/
c$$$     -    (xmax*ymax + eps*Pi*(xmax - ymax)*
c$$$     -       Sin(Pi*((2*x1)/xmax - (2*x2)/ymax)) + 
c$$$     -      eps*Pi*(xmax + ymax)*Sin(2*Pi*(x1/xmax + x2/ymax)))
c$$$
c$$$          vec(3) = 0d0
c$$$          vec(4) = vec(2)
c$$$
c$$$          vec(5) =
c$$$     .         (xmax*(2*eps**2*Pi**2 + ymax**2 - 
c$$$     -        2*eps**2*Pi**2*Cos((4*Pi*x1)/xmax) - 
c$$$     -        eps**2*Pi**2*Cos(Pi*((4*x1)/xmax - (4*x2)/ymax)) - 
c$$$     -        eps**2*Pi**2*Cos(4*Pi*(x1/xmax + x2/ymax)) + 
c$$$     -        2*eps**2*Pi**2*Cos((4*Pi*x2)/ymax) + 
c$$$     -        2*eps*Pi*ymax*Sin(Pi*((2*x1)/xmax - (2*x2)/ymax)) + 
c$$$     -        2*eps*Pi*ymax*Sin(2*Pi*(x1/xmax + x2/ymax))))/
c$$$     -    (ymax*(xmax*ymax + eps*Pi*(xmax - ymax)*
c$$$     -         Sin(Pi*((2*x1)/xmax - (2*x2)/ymax)) + 
c$$$     -        eps*Pi*(xmax + ymax)*Sin(2*Pi*(x1/xmax + x2/ymax))))!
c$$$
c$$$          vec(6) = 0d0
c$$$          vec(7) = vec(3)
c$$$          vec(8) = vec(6)
c$$$
c$$$          vec(9) = (xmax*ymax)/
c$$$     -    (xmax*ymax + eps*Pi*(xmax - ymax)*
c$$$     -       Sin(Pi*((2*x1)/xmax - (2*x2)/ymax)) + 
c$$$     -      eps*Pi*(xmax + ymax)*Sin(2*Pi*(x1/xmax + x2/ymax)))
c$$$
c$$$        case default
c$$$          write (*,*) 'Grid not implemented in g_sub'
c$$$          write (*,*) 'Aborting...'
c$$$          stop
c$$$        end select
c$$$
c$$$        tensor = reshape(vec,(/3,3/))
c$$$
c$$$      end function g_sub

c     g_sub
c     #################################################################
      function g_sub(gl_def,i,j,k,igx,igy,igz,vertex) result (gsub)

c     -----------------------------------------------------------------
c     Calculates contravariant metric tensor of curvilinear coordinate 
c     system
c     -----------------------------------------------------------------

        implicit none

c     Input variables

        integer :: i,j,k,igx,igy,igz
        real(8) :: gsub(3,3)
        logical,optional :: vertex

        type(grid_mg_def),pointer :: gl_def

c     Local variables

        integer :: m,n
        real(8) :: cnv(3,3),jac

c     Begin program

        do m = 1,3
          cnv(m,:) = contravariantVector(gl_def,m,i,j,k,igx,igy,igz
     .                                  ,vertex=vertex)
        enddo

        jac = jacobian(gl_def,i,j,k,igx,igy,igz,vertex=vertex)

        do m=1,3
          do n=m,3
            gsub(n,m) = jac*dot_product(cnv(n,:),cnv(m,:))
            gsub(m,n) = gsub(n,m) !Symmetry
          enddo
        enddo

      end function g_sub

c$$$c     g_sup
c$$$c     #################################################################
c$$$      function g_sup(gl_def,i,j,k,igx,igy,igz,vertex) result (tensor)
c$$$
c$$$c     -----------------------------------------------------------------
c$$$c     Calculates covariant metric tensor of curvilinear coordinate 
c$$$c     system
c$$$c     -----------------------------------------------------------------
c$$$
c$$$        implicit none
c$$$
c$$$c     Input variables
c$$$
c$$$        integer :: i,j,k,igx,igy,igz
c$$$        real(8) :: tensor(3,3)
c$$$        logical,optional :: vertex
c$$$
c$$$        type(grid_mg_def),pointer :: gl_def
c$$$
c$$$c     Local variables
c$$$
c$$$        integer :: ig,jg,kg
c$$$        real(8) :: vec(9),jac,car(3),curv(3)
c$$$        real(8) :: x1,x2,x3,aa,bb,major_r,mm,kk,eps,gp(ngparams)
c$$$     .            ,x1p,x2p,x3p
c$$$        logical :: vrtx
c$$$
c$$$        real(8) :: xmax,ymax,zmax,xmin,ymin,zmin
c$$$! Daniele, 11/08/2017. Local variables needed for shaped toroidal coordinates
c$$$        integer :: ii
c$$$        real(8) :: ai(10),bi(10)
c$$$        real(8) :: ap,bp
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
c$$$        if (PRESENT(vertex)) then
c$$$          vrtx = vertex
c$$$        else
c$$$          vrtx = .false.
c$$$        endif
c$$$
c$$$        gp = gl_def%params
c$$$
c$$$        call getCurvilinearCoordinates(gl_def,i,j,k,igx,igy,igz,ig,jg,kg
c$$$     .                                ,x1,x2,x3)
c$$$
c$$$        if (vrtx) then
c$$$          call getCurvilinearCoordinates(gl_def,i+1,j+1,k+1,igx,igy,igz
c$$$     .                                  ,ig,jg,kg,x1p,x2p,x3p)
c$$$          x1 = 0.5*(x1+x1p)
c$$$          x2 = 0.5*(x2+x2p)
c$$$          x3 = 0.5*(x3+x3p)
c$$$        endif
c$$$
c$$$        select case (coords)
c$$$        case ('car','lpk')
c$$$          vec = (/ 1d0, 0d0, 0d0
c$$$     .            ,0d0, 1d0, 0d0
c$$$     .            ,0d0, 0d0, 1d0 /)
c$$$        case ('cyl')
c$$$          if (x1 == 0d0) then
c$$$            vec = 0d0
c$$$          else
c$$$            vec = (/ x1 , 0d0   , 0d0
c$$$     .              ,0d0, 1d0/x1, 0d0
c$$$     .              ,0d0, 0d0   , x1 /)
c$$$          endif
c$$$        case ('hel')
c$$$          mm = gp(1)
c$$$          kk = gp(2)
c$$$          if (x1 == 0d0) then
c$$$            vec = 0d0
c$$$          else
c$$$            vec = (/  x1/mm , 0d0                , 0d0
c$$$     .              , 0d0   , mm/x1 + kk*kk*x1/mm, kk*x1/mm
c$$$     .              , 0d0   , kk*x1/mm           , x1/mm /)
c$$$          endif
c$$$        case ('hl2')
c$$$          mm = gp(1)
c$$$          kk = gp(2)
c$$$          aa = kk/mm
c$$$
c$$$          x1 = abs(2*x1)
c$$$
c$$$          if (x1 == 0d0) then
c$$$            vec = 0d0
c$$$          else
c$$$            vec = (/ x1 , 0d0          , 0d0
c$$$     .              ,0d0, 1./x1 + aa**2, aa
c$$$     .              ,0d0, aa           , 1d0 /)
c$$$          endif
c$$$        case ('tor')
c$$$          if (x1 == 0d0) then
c$$$            vec = 0d0
c$$$          else
c$$$            major_r = gp(1)
c$$$            aa      = gp(2)
c$$$            bb      = gp(3)
c$$$
c$$$            vec(1) = -(x1*(major_r + aa*x1*Cos(x2))*
c$$$     .               (-aa**2 - bb**2 + (aa**2 - bb**2)*Cos(2*x2)))
c$$$     .               /(2.*aa*bb)
c$$$
c$$$            vec(2) = ((aa**2 - bb**2)*Cos(x2)
c$$$     .               *(major_r + aa*x1*Cos(x2))*Sin(x2))
c$$$     .               /(aa*bb)
c$$$
c$$$            vec(3) = 0d0
c$$$
c$$$            vec(4) = vec(2)
c$$$
c$$$            vec(5) = ((major_r + aa*x1*Cos(x2))
c$$$     .              *(aa**2 + bb**2 + (aa**2 - bb**2)*Cos(2*x2)))
c$$$     .              /(2.*aa*bb*x1)
c$$$
c$$$            vec(6) = 0d0
c$$$
c$$$            vec(7) = vec(3)
c$$$            vec(8) = vec(6)
c$$$            vec(9) = (aa*bb*x1)/(major_r + aa*x1*Cos(x2))
c$$$          endif
c$$$        case ('sha')
c$$$! Daniele, 11/08/2017
c$$$! Shaped toroidal coordinates system
c$$$! x = (R0 + r*sum_i(ai*cos(i*theta))*cos(-phi)
c$$$! y = (R0 + r*sum_i(ai*cos(i*theta))*sin(-phi)
c$$$! z = r*sum_i(bi*sin(i*theta)
c$$$          if (x1 == 0d0) then
c$$$            vec = 0d0
c$$$          else
c$$$          major_r = gp(1)  !Major radius
c$$$
c$$$          ai(1) = gp(2)  !Horizontal minor radius
c$$$          bi(1) = gp(3)  !Vertical   minor radius
c$$$
c$$$          do ii=2,10
c$$$             ai(ii) = gp(7+ii) ! Higher order cosine coefficients are stored in gp(9) to gp(17)
c$$$             bi(ii) = gp(16+ii) ! Higher order sine coefficients are stored in gp(18) to gp(26)
c$$$          enddo
c$$$
c$$$          aa = 0d0
c$$$          bb = 0d0
c$$$          ap = 0d0
c$$$          bp = 0d0
c$$$          do ii=1,10
c$$$             aa = aa + ai(ii)*cos(ii*x2) ! A(theta)=sum_i(ai*cos(i*theta))
c$$$             bb = bb + bi(ii)*sin(ii*x2) ! B(theta)=sum_i(bi*sin(i*theta))
c$$$             bp = bp + bi(ii)*ii*cos(ii*x2) ! B'(theta)
c$$$             ap = ap - ai(ii)*ii*sin(ii*x2) ! A'(theta)
c$$$          enddo
c$$$
c$$$            vec(1) = x1*(major_r + aa*x1)*(ap**2 + bp**2)/
c$$$     .               (aa*bp-bb*ap)
c$$$
c$$$            vec(2) = (major_r + aa*x1)*(-aa*ap - bb*bp)/
c$$$     .               (aa*bp-bb*ap)
c$$$
c$$$            vec(3) = 0d0
c$$$
c$$$            vec(4) = vec(2)
c$$$
c$$$            vec(5) = (major_r + aa*x1)*(aa**2 + bb**2)/
c$$$     .               (x1*(aa*bp-bb*ap))
c$$$
c$$$            vec(6) = 0d0
c$$$
c$$$            vec(7) = vec(3)
c$$$            vec(8) = vec(6)
c$$$            vec(9) = x1*(aa*bp-bb*ap)/(major_r + aa*x1)
c$$$
c$$$          endif
c$$$        case ('sin')
c$$$
c$$$          eps = gp(1)
c$$$
c$$$          vec(1) =
c$$$     .         (xmax*(2*eps**2*Pi**2 + ymax**2 - 
c$$$     -        2*eps**2*Pi**2*Cos((4*Pi*x1)/xmax) - 
c$$$     -        eps**2*Pi**2*Cos(Pi*((4*x1)/xmax - (4*x2)/ymax)) - 
c$$$     -        eps**2*Pi**2*Cos(4*Pi*(x1/xmax + x2/ymax)) + 
c$$$     -        2*eps**2*Pi**2*Cos((4*Pi*x2)/ymax) + 
c$$$     -        2*eps*Pi*ymax*Sin(Pi*((2*x1)/xmax - (2*x2)/ymax)) + 
c$$$     -        2*eps*Pi*ymax*Sin(2*Pi*(x1/xmax + x2/ymax))))/
c$$$     -    (ymax*(xmax*ymax + eps*Pi*(xmax - ymax)*
c$$$     -         Sin(Pi*((2*x1)/xmax - (2*x2)/ymax)) + 
c$$$     -        eps*Pi*(xmax + ymax)*Sin(2*Pi*(x1/xmax + x2/ymax))))
c$$$
c$$$          vec(2) =
c$$$     .         -((eps*Pi*(eps*Pi*Cos(Pi*((4*x1)/xmax - (4*x2)/ymax)) - 
c$$$     -          eps*Pi*Cos(4*Pi*(x1/xmax + x2/ymax)) + 
c$$$     -          xmax*Sin(Pi*((2*x1)/xmax - (2*x2)/ymax)) - 
c$$$     -          ymax*Sin(Pi*((2*x1)/xmax - (2*x2)/ymax)) + 
c$$$     -          xmax*Sin(2*Pi*(x1/xmax + x2/ymax)) + 
c$$$     -          ymax*Sin(2*Pi*(x1/xmax + x2/ymax))))/
c$$$     -      (xmax*ymax + eps*Pi*(xmax - ymax)*
c$$$     -         Sin(Pi*((2*x1)/xmax - (2*x2)/ymax)) + 
c$$$     -        eps*Pi*(xmax + ymax)*Sin(2*Pi*(x1/xmax + x2/ymax))))
c$$$
c$$$          vec(3) = 0d0
c$$$          vec(4) = vec(2)
c$$$
c$$$          vec(5) =
c$$$     .         (ymax*(2*eps**2*Pi**2 + xmax**2 + 
c$$$     -        2*eps**2*Pi**2*Cos((4*Pi*x1)/xmax) - 
c$$$     -        eps**2*Pi**2*Cos(Pi*((4*x1)/xmax - (4*x2)/ymax)) - 
c$$$     -        eps**2*Pi**2*Cos(4*Pi*(x1/xmax + x2/ymax)) - 
c$$$     -        2*eps**2*Pi**2*Cos((4*Pi*x2)/ymax) - 
c$$$     -        2*eps*Pi*xmax*Sin(Pi*((2*x1)/xmax - (2*x2)/ymax)) + 
c$$$     -        2*eps*Pi*xmax*Sin(2*Pi*(x1/xmax + x2/ymax))))/
c$$$     -    (xmax*(xmax*ymax + eps*Pi*(xmax - ymax)*
c$$$     -         Sin(Pi*((2*x1)/xmax - (2*x2)/ymax)) + 
c$$$     -        eps*Pi*(xmax + ymax)*Sin(2*Pi*(x1/xmax + x2/ymax))))
c$$$
c$$$          vec(6) = 0d0
c$$$          vec(7) = vec(3)
c$$$          vec(8) = vec(6)
c$$$
c$$$          vec(9) = (xmax*ymax + eps*Pi*(xmax - ymax)*
c$$$     -       Sin(Pi*((2*x1)/xmax - (2*x2)/ymax)) + 
c$$$     -      eps*Pi*(xmax + ymax)*Sin(2*Pi*(x1/xmax + x2/ymax)))
c$$$     .         /(xmax*ymax)
c$$$
c$$$        case default
c$$$          write (*,*) 'Grid not implemented in g_sup'
c$$$          write (*,*) 'Aborting...'
c$$$          stop
c$$$        end select
c$$$
c$$$        tensor = reshape(vec,(/3,3/))
c$$$
c$$$      end function g_sup

c$$$c     gsup_anal
c$$$c     #################################################################
c$$$      function gsup_anal(gl_def,curv) result(tensor)
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
c$$$        real(8) :: curv(3),tensor(3,3)
c$$$
c$$$        type(grid_mg_def),pointer :: gl_def
c$$$
c$$$c     Local variables
c$$$
c$$$        real(8) :: aa,bb,major_r,mm,kk,eps,phi,jac,x1,x2,x3
c$$$        real(8) :: gp(ngparams),vec(9)
c$$$
c$$$        real(8) :: xmax,ymax,zmax,xmin,ymin,zmin
c$$$
c$$$c     Begin program
c$$$
c$$$        gp = gl_def%params
c$$$
c$$$        x1 = curv(1)
c$$$        x2 = curv(2)
c$$$        x3 = curv(3)
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
c$$$        select case (coords)
c$$$        case ('car','lpk')
c$$$          vec = (/ 1d0, 0d0, 0d0
c$$$     .            ,0d0, 1d0, 0d0
c$$$     .            ,0d0, 0d0, 1d0 /)
c$$$        case ('cyl')
c$$$          if (x1 == 0d0) then
c$$$            vec = 0d0
c$$$          else
c$$$            vec = (/ x1 , 0d0   , 0d0
c$$$     .              ,0d0, 1d0/x1, 0d0
c$$$     .              ,0d0, 0d0   , x1 /)
c$$$          endif
c$$$        case ('hel')
c$$$          mm = gp(1)
c$$$          kk = gp(2)
c$$$          if (x1 == 0d0) then
c$$$            vec = 0d0
c$$$          else
c$$$            vec = (/  x1/mm , 0d0                , 0d0
c$$$     .              , 0d0   , mm/x1 + kk*kk*x1/mm, kk*x1/mm
c$$$     .              , 0d0   , kk*x1/mm           , x1/mm /)
c$$$          endif
c$$$        case ('hl2')
c$$$          mm = gp(1)
c$$$          kk = gp(2)
c$$$          aa = kk/mm
c$$$
c$$$          x1 = abs(2*x1)
c$$$
c$$$          if (x1 == 0d0) then
c$$$            vec = 0d0
c$$$          else
c$$$            vec = (/ x1 , 0d0          , 0d0
c$$$     .              ,0d0, 1./x1 + aa**2, aa
c$$$     .              ,0d0, aa           , 1d0 /)
c$$$          endif
c$$$        case ('tor')
c$$$          if (x1 == 0d0) then
c$$$            vec = 0d0
c$$$          else
c$$$            major_r = gp(1)
c$$$            aa      = gp(2)
c$$$            bb      = gp(3)
c$$$
c$$$            vec(1) = -(x1*(major_r + aa*x1*Cos(x2))*
c$$$     .               (-aa**2 - bb**2 + (aa**2 - bb**2)*Cos(2*x2)))
c$$$     .               /(2.*aa*bb)
c$$$
c$$$            vec(2) = ((aa**2 - bb**2)*Cos(x2)
c$$$     .               *(major_r + aa*x1*Cos(x2))*Sin(x2))
c$$$     .               /(aa*bb)
c$$$
c$$$            vec(3) = 0d0
c$$$
c$$$            vec(4) = vec(2)
c$$$
c$$$            vec(5) = ((major_r + aa*x1*Cos(x2))
c$$$     .              *(aa**2 + bb**2 + (aa**2 - bb**2)*Cos(2*x2)))
c$$$     .              /(2.*aa*bb*x1)
c$$$
c$$$            vec(6) = 0d0
c$$$
c$$$            vec(7) = vec(3)
c$$$            vec(8) = vec(6)
c$$$            vec(9) = (aa*bb*x1)/(major_r + aa*x1*Cos(x2))
c$$$          endif
c$$$        case ('sin')
c$$$
c$$$          eps = gp(1)
c$$$
c$$$          vec(1) =
c$$$     .         (xmax*(2*eps**2*Pi**2 + ymax**2 - 
c$$$     -        2*eps**2*Pi**2*Cos((4*Pi*x1)/xmax) - 
c$$$     -        eps**2*Pi**2*Cos(Pi*((4*x1)/xmax - (4*x2)/ymax)) - 
c$$$     -        eps**2*Pi**2*Cos(4*Pi*(x1/xmax + x2/ymax)) + 
c$$$     -        2*eps**2*Pi**2*Cos((4*Pi*x2)/ymax) + 
c$$$     -        2*eps*Pi*ymax*Sin(Pi*((2*x1)/xmax - (2*x2)/ymax)) + 
c$$$     -        2*eps*Pi*ymax*Sin(2*Pi*(x1/xmax + x2/ymax))))/
c$$$     -    (ymax*(xmax*ymax + eps*Pi*(xmax - ymax)*
c$$$     -         Sin(Pi*((2*x1)/xmax - (2*x2)/ymax)) + 
c$$$     -        eps*Pi*(xmax + ymax)*Sin(2*Pi*(x1/xmax + x2/ymax))))
c$$$
c$$$          vec(2) =
c$$$     .         -((eps*Pi*(eps*Pi*Cos(Pi*((4*x1)/xmax - (4*x2)/ymax)) - 
c$$$     -          eps*Pi*Cos(4*Pi*(x1/xmax + x2/ymax)) + 
c$$$     -          xmax*Sin(Pi*((2*x1)/xmax - (2*x2)/ymax)) - 
c$$$     -          ymax*Sin(Pi*((2*x1)/xmax - (2*x2)/ymax)) + 
c$$$     -          xmax*Sin(2*Pi*(x1/xmax + x2/ymax)) + 
c$$$     -          ymax*Sin(2*Pi*(x1/xmax + x2/ymax))))/
c$$$     -      (xmax*ymax + eps*Pi*(xmax - ymax)*
c$$$     -         Sin(Pi*((2*x1)/xmax - (2*x2)/ymax)) + 
c$$$     -        eps*Pi*(xmax + ymax)*Sin(2*Pi*(x1/xmax + x2/ymax))))
c$$$
c$$$          vec(3) = 0d0
c$$$          vec(4) = vec(2)
c$$$
c$$$          vec(5) =
c$$$     .         (ymax*(2*eps**2*Pi**2 + xmax**2 + 
c$$$     -        2*eps**2*Pi**2*Cos((4*Pi*x1)/xmax) - 
c$$$     -        eps**2*Pi**2*Cos(Pi*((4*x1)/xmax - (4*x2)/ymax)) - 
c$$$     -        eps**2*Pi**2*Cos(4*Pi*(x1/xmax + x2/ymax)) - 
c$$$     -        2*eps**2*Pi**2*Cos((4*Pi*x2)/ymax) - 
c$$$     -        2*eps*Pi*xmax*Sin(Pi*((2*x1)/xmax - (2*x2)/ymax)) + 
c$$$     -        2*eps*Pi*xmax*Sin(2*Pi*(x1/xmax + x2/ymax))))/
c$$$     -    (xmax*(xmax*ymax + eps*Pi*(xmax - ymax)*
c$$$     -         Sin(Pi*((2*x1)/xmax - (2*x2)/ymax)) + 
c$$$     -        eps*Pi*(xmax + ymax)*Sin(2*Pi*(x1/xmax + x2/ymax))))
c$$$
c$$$          vec(6) = 0d0
c$$$          vec(7) = vec(3)
c$$$          vec(8) = vec(6)
c$$$
c$$$          vec(9) = (xmax*ymax + eps*Pi*(xmax - ymax)*
c$$$     -       Sin(Pi*((2*x1)/xmax - (2*x2)/ymax)) + 
c$$$     -      eps*Pi*(xmax + ymax)*Sin(2*Pi*(x1/xmax + x2/ymax)))
c$$$     .         /(xmax*ymax)
c$$$
c$$$        case default
c$$$          write (*,*) 'Grid not implemented in g_sup'
c$$$          write (*,*) 'Aborting...'
c$$$          stop
c$$$        end select
c$$$
c$$$        tensor = reshape(vec,(/3,3/))
c$$$
c$$$      end function gsup_anal

c     g_sup
c     #################################################################
      function g_sup(gl_def,i,j,k,igx,igy,igz,vertex) result (gsup)

c     -----------------------------------------------------------------
c     Calculates contravariant metric tensor of curvilinear coordinate 
c     system
c     -----------------------------------------------------------------

        implicit none

c     Input variables

        integer :: i,j,k,igx,igy,igz
        real(8) :: gsup(3,3)
        logical,optional :: vertex

        type(grid_mg_def),pointer :: gl_def

c     Local variables

        integer :: m,n
        real(8) :: cov(3,3),jac

c     Begin program

        do m = 1,3
          cov(m,:) = covariantVector(gl_def,m,i,j,k,igx,igy,igz
     .                              ,vertex=vertex)
        enddo

        jac = jacobian(gl_def,i,j,k,igx,igy,igz,vertex=vertex)

        do m=1,3
          do n=m,3
            gsup(n,m) = jac*dot_product(cov(n,:),cov(m,:))
            gsup(m,n) = gsup(n,m) !Symmetry
          enddo
        enddo

        end function g_sup
      
      end module grid_anal_map
