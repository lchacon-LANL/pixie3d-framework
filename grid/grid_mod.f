c module grid_core
c #####################################################################
      module grid_core

        implicit none

        real(8)    :: gparams(5)

        character*(3) :: coords

        real(8)    :: xmax,ymax,zmax,xmin,ymin,zmin  !3D grid dimension

        real(8),private :: pi,lambda,cc,ypp,eps,mm,kk,aa,phi,major_r

      contains

c     xi_x
c     #################################################################
      function xi_x(xx,yy,zz) result(curv)

c     -----------------------------------------------------------------
c     Gives curvilinear coordinates from Cartesian coordinates
c     -----------------------------------------------------------------

        implicit none

c     Input variables

        real(8)    :: xx,yy,zz,curv(3)

c     Local variables

        real(8)    :: x1,x2,x3

c     Begin program

        select case (coords)
        case ('car')
          x1 = xx
          x2 = yy
          x3 = zz
        case ('scl')
          lambda = gparams(1)
          cc = 0.5/lambda
          cc = 1./tanh(cc)

          x1 = xx
          x2 = 0.5*ymax*(1.+cc*tanh((yy-0.5)/lambda))
          x3 = zz
        case ('cyl')
          x1 = sqrt(xx**2 + yy**2)
          x2 = atan(yy/xx)
          x3 = zz
        case ('hel')
          mm = gparams(1)
          kk = gparams(2)
          aa = kk/mm

          x1 = sqrt(xx**2 + yy**2)
          x2 = atan(yy/xx) + aa*zz
          x3 = zz
        case ('tor')
          major_r = gparams(1)

          x1 = sqrt( zz**2 + (sqrt(xx**2 + yy**2) - major_r)**2 )
          x2 = atan( (sqrt(xx**2 + yy**2) - major_r)/zz)
          x3 = atan( yy/xx )
        case ('sin')

          pi = acos(-1d0)

          eps = gparams(1)
          x1 = xx + eps*sin(2*pi*xx/xmax)*sin(2*pi*yy/ymax)
          x2 = yy + eps*sin(2*pi*xx/xmax)*sin(2*pi*yy/ymax)
          x3 = zz
        case default
          write (*,*) 'Grid not implemented'
          write (*,*) 'Aborting...'
          stop
        end select

        curv = (/ x1,x2,x3 /)

      end function xi_x

c     x_xi
c     #################################################################
      function x_xi(x1,x2,x3) result(cartesian)

c     -----------------------------------------------------------------
c     Gives Cartesian coordinates from curvilinear coordinates
c     -----------------------------------------------------------------

        implicit none

c     Input variables

        real(8)    :: x1,x2,x3,cartesian(3)

c     Local variables

        integer(4) :: inewt,ic
        real(8)    :: xx,yy,zz,jac_mat(3,3),rhs(3),dx(3),rr,rr0

c     Begin program

        select case (coords)
        case ('car')
          xx = x1
          yy = x2
          zz = x3
        case ('scl')
          lambda = gparams(1)
          cc = 0.5/lambda
          cc = 1./tanh(cc)
          ypp = (2*x2/ymax-1.)

          xx = x1
          yy = 0.5+lambda*atanh(ypp/cc)
          zz = x3
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
        case ('tor')
          major_r = gparams(1)

          xx = (major_r + x1*sin(x2))*cos(x3)
          yy = (major_r + x1*sin(x2))*sin(x3)
          zz = x1*cos(x2)
        case ('sin')

          !Initial guess
          xx = x1
          yy = x2
          zz = x3

          !Initial residual
          rhs = -( (/x1,x2,x3/) - xi_x(xx,yy,zz) )

          rr0 = sqrt(sum(rhs*rhs))
          rr  = rr0

cc          write (*,*)

          do inewt = 1,10

cc            write (*,*) xx,yy,zz,rr

            if (rr < (1d-10 + 1d-4*rr0)) exit

          !Jacobian
            do ic =1,3
              jac_mat(ic,:) = -covariantVector(ic,xx,yy,zz,.true.)
            enddo

          !Solve Jacobian system

            call solve(3,jac_mat,rhs,dx)

          !Update Newton solution

            xx = xx + dx(1)
            yy = yy + dx(2)
            zz = zz + dx(3)

          !Check residual

            rhs = -( (/x1,x2,x3/) - xi_x(xx,yy,zz) )

            rr = sqrt(sum(rhs*rhs))

          enddo

          if (inewt >= 10) then
            write (*,*) 'Newton did not converge in',inewt,' iterations'
            write (*,*) 'Problem in x_xi'
            write (*,*) 'Aborting'
            stop
          endif

        case default
          write (*,*) 'Grid not implemented'
          write (*,*) 'Aborting...'
          stop
         end select

        cartesian = (/ xx,yy,zz /)

      contains

c     atanh
c     #################################################################
      real(8) function atanh(x)

        real(8) :: x

        atanh = 0.5*(log( (1+x)/(1-x) ) )

      end function atanh

c     solve
c     #################################################################
      subroutine solve(size,mat,rhs,x)

        implicit none

        integer(4) :: size

        real(8)    :: mat(size,size),rhs(size),x(size)

        integer(4) :: ipiv(size),info

        external dgesv

        x = rhs

        call dgesv(size,1,mat,size,ipiv,x,size,info)  !LAPACK routine

        if (info /= 0) then
          if (info < 0) then
            write (*,*) 'Problem in factorization in argument',-info
            write (*,*) 'Error in x_xi'
            write (*,*) 'Aborting'
            stop
          else
            write (*,*) 'Matrix is singular'
            write (*,*) 'Error in x_xi'
            write (*,*) 'Aborting'
            stop
          endif
        endif

      end subroutine solve

      end function x_xi

c     jacobian
c     #################################################################
      function jacobian(x1,x2,x3,cartesian) result(jac)

c     -----------------------------------------------------------------
c     Calculates Jacobian of curvilinear coordinate system
c     -----------------------------------------------------------------

        implicit none

c     Input variables

        real(8) :: x1,x2,x3,jac
        logical :: cartesian

c     Local variables

        real(8) :: car(3),curv(3)

c     Begin program

        select case (coords)
        case ('car')
          jac = 1d0
        case ('scl')
          if (cartesian) then
            curv = xi_x(x1,x2,x3)
          else
            curv = (/ x1,x2,x3 /)
          endif
          lambda = gparams(1)
          cc = 0.5/lambda
          cc = 1./tanh(cc)
          ypp = (2*curv(2)/ymax-1.)
          jac = cc*lambda/(cc**2-ypp**2)
        case ('cyl')
          if (cartesian) then
            curv = xi_x(x1,x2,x3)
          else
            curv = (/ x1,x2,x3 /)
          endif
          jac = curv(1)
        case ('hel')
          if (cartesian) then
            curv = xi_x(x1,x2,x3)
          else
            curv = (/ x1,x2,x3 /)
          endif
          jac = curv(1)
        case ('tor')
          if (cartesian) then
            curv = xi_x(x1,x2,x3)
          else
            curv = (/ x1,x2,x3 /)
          endif
          major_r = gparams(1)
          jac = curv(1)*(major_r + curv(1)*sin(curv(2)))
        case ('sin')

          pi = acos(-1d0)

          if (cartesian) then
            car = (/ x1,x2,x3 /)
          else
            car = x_xi(x1,x2,x3)
          endif
          eps = gparams(1)
cc          jac = 1./(1 + 2*Pi*eps*Sin(2*Pi*(car(1) + car(2))))
          jac = (xmax*ymax)/
     -          (xmax*ymax +
     -        eps*Pi*xmax*Sin((2*car(1)*Pi)/xmax - (2*car(2)*Pi)/ymax) - 
     -        eps*Pi*ymax*Sin((2*car(1)*Pi)/xmax - (2*car(2)*Pi)/ymax) + 
     -        eps*Pi*xmax*Sin((2*car(1)*Pi)/xmax + (2*car(2)*Pi)/ymax) + 
     -        eps*Pi*ymax*Sin((2*car(1)*Pi)/xmax + (2*car(2)*Pi)/ymax))
        case default
          write (*,*) 'Grid not implemented'
          write (*,*) 'Aborting...'
          stop
        end select

      end function jacobian

c     covariantVector
c     #################################################################
      function covariantVector(comp,x1,x2,x3,cartesian) result (vec)

c     -----------------------------------------------------------------
c     Calculates covariant vectors of curvilinear coordinate system
c     in Cartesian coordinates
c     -----------------------------------------------------------------

        implicit none

c     Input variables

        integer(4) :: comp
        real(8)    :: x1,x2,x3,vec(3)
        logical    :: cartesian

c     Local variables

        real(8)    :: car(3),curv(3),jac

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
          if (cartesian) then
            curv = xi_x(x1,x2,x3)
          else
            curv = (/ x1,x2,x3 /)
          endif
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
          if (cartesian) then
            curv = xi_x(x1,x2,x3)
          else
            curv = (/ x1,x2,x3 /)
          endif
          select case (comp)
            case (1)
              vec = (/ cos(curv(2)),sin(curv(2)),0d0 /)
            case (2)
              vec = (/-sin(curv(2)),cos(curv(2)),0d0 /)/curv(1)
            case (3)
              vec = (/ 0d0,0d0,1d0 /)
          end select
        case ('hel')
          if (cartesian) then
            curv = xi_x(x1,x2,x3)
          else
            curv = (/ x1,x2,x3 /)
          endif
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
        case ('tor')
          major_r = gparams(1)
          if (cartesian) then
            curv = xi_x(x1,x2,x3)
          else
            curv = (/ x1,x2,x3 /)
          endif
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

          if (cartesian) then
            car = (/ x1,x2,x3 /)
          else
            car = x_xi(x1,x2,x3)
          endif
          eps = gparams(1)
          select case (comp)
            case (1)
              vec = (/ 1 + (2*eps*Pi*Cos((2*car(1)*Pi)/xmax)
     .                              *Sin((2*car(2)*Pi)/ymax))/xmax
     .                ,    (2*eps*Pi*Cos((2*car(2)*Pi)/ymax)
     .                              *Sin((2*car(1)*Pi)/xmax))/ymax
     .                ,0d0 /)
            case (2)
              vec = (/     (2*eps*Pi*Cos((2*car(1)*Pi)/xmax)
     .                              *Sin((2*car(2)*Pi)/ymax))/xmax
     .                ,1 + (2*eps*Pi*Cos((2*car(2)*Pi)/ymax)
     .                              *Sin((2*car(1)*Pi)/xmax))/ymax
     .                ,0d0 /)
            case (3)
              vec = (/ 0d0,0d0,1d0 /)
          end select
        case default
          write (*,*) 'Grid not implemented'
          write (*,*) 'Aborting...'
          stop
        end select

      end function covariantVector

c     contravariantVector
c     #################################################################
      function contravariantVector(comp,x1,x2,x3,cartesian) result (vec)

c     -----------------------------------------------------------------
c     Calculates contravariant vectors of curvilinear coordinate system
c     in Cartesian coordinates
c     -----------------------------------------------------------------

        implicit none

c     Input variables

        integer(4) :: comp
        real(8)    :: x1,x2,x3,vec(3)
        logical    :: cartesian

c     Local variables

        real(8)    :: car(3),curv(3),jac

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
          if (cartesian) then
            curv = xi_x(x1,x2,x3)
          else
            curv = (/ x1,x2,x3 /)
          endif
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
          if (cartesian) then
            curv = xi_x(x1,x2,x3)
          else
            curv = (/ x1,x2,x3 /)
          endif
          select case (comp)
            case (1)
              vec = (/ cos(curv(2)),sin(curv(2)),0d0 /)/curv(1)
            case (2)
              vec = (/-sin(curv(2)),cos(curv(2)),0d0 /)
            case (3)
              vec = (/ 0d0,0d0,1d0 /)/curv(1)
          end select
        case ('hel')
          if (cartesian) then
            curv = xi_x(x1,x2,x3)
          else
            curv = (/ x1,x2,x3 /)
          endif
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
        case ('tor')
          major_r = gparams(1)
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

          if (cartesian) then
            car = (/ x1,x2,x3 /)
          else
            car = x_xi(x1,x2,x3)
          endif
          eps = gparams(1)
          select case (comp)
            case (1)
              vec = (/ 1 + (2*eps*Pi*Cos((2*car(2)*Pi)/ymax)
     .                              *Sin((2*car(1)*Pi)/xmax))/ymax
     .                ,   -(2*eps*Pi*Cos((2*car(1)*Pi)/xmax)
     .                              *Sin((2*car(2)*Pi)/ymax))/xmax
     .                ,0d0 /)
            case (2)
              vec = (/    -(2*eps*Pi*Cos((2*car(2)*Pi)/ymax)
     .                              *Sin((2*car(1)*Pi)/xmax))/ymax
     .                ,1 + (2*eps*Pi*Cos((2*car(1)*Pi)/xmax)
     .                              *Sin((2*car(2)*Pi)/ymax))/xmax
     .                ,0d0 /)
            case (3)
              vec = (/ 0d0
     .                ,0d0
     .                ,1 + (2*eps*Pi*Cos((2*car(2)*Pi)/ymax)
     .                              *Sin((2*car(1)*Pi)/xmax))/ymax 
     -                   + (2*eps*Pi*Cos((2*car(1)*Pi)/xmax)
     .                              *Sin((2*car(2)*Pi)/ymax))/xmax/)
          end select
        case default
          write (*,*) 'Grid not implemented'
          write (*,*) 'Aborting...'
          stop
        end select

      end function contravariantVector

c     hessian
c     #################################################################
      function hessian(k,x1,x2,x3,cartesian) result (tensor)

c     -----------------------------------------------------------------
c     Calculates hessian elements of curvilinear coordinate system in
c     the covariant basis, i.e., 
c            hessian[k](i,j) = J^2 <cnv(i)|grad(cov[k])|cnv(j)>
c     -----------------------------------------------------------------

        implicit none

c     Input variables

        integer(4) :: k
        real(8)    :: x1,x2,x3,tensor(3,3)
        logical    :: cartesian

c     Local variables

        real(8)    :: car(3),curv(3),vec(9)

c     Begin program

        select case (coords)
        case ('car')
          vec = (/ 0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0 /)
        case ('scl')
          if (cartesian) then
            curv = xi_x(x1,x2,x3)
          else
            curv = (/ x1,x2,x3 /)
          endif
          lambda = gparams(1)
          cc = 0.5/lambda
          cc = 1./tanh(cc)
          ypp = (2*curv(2)/ymax-1.)

          select case (k)
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
        case ('cyl')
          if (cartesian) then
            curv = xi_x(x1,x2,x3)
          else
            curv = (/ x1,x2,x3 /)
          endif
          select case (k)
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
        case ('hel')
          if (cartesian) then
            curv = xi_x(x1,x2,x3)
          else
            curv = (/ x1,x2,x3 /)
          endif
          mm = gparams(1)
          kk = gparams(2)
          aa = kk/mm
          phi = (curv(2)-aa*curv(3))
          select case (k)
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
        case ('tor')
          if (cartesian) then
            curv = xi_x(x1,x2,x3)
          else
            curv = (/ x1,x2,x3 /)
          endif
          major_r = gparams(1)
          select case (k)
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
        case ('sin')

          pi = acos(-1d0)

          if (cartesian) then
            car = (/ x1,x2,x3 /)
          else
            car = x_xi(x1,x2,x3)
          endif
          eps = gparams(1)
          select case (k)
            case (1)
              vec(1) =
     .    (4*(eps**3*Pi**4*Cos((2*car(1)*Pi)/xmax - (6*car(2)*Pi)/ymax)+ 
     -      2*eps**3*Pi**4*Cos((2*car(1)*Pi)/xmax - (2*car(2)*Pi)/ymax)+ 
     -      eps*Pi**2*ymax**2*
     -       Cos((2*car(1)*Pi)/xmax - (2*car(2)*Pi)/ymax) + 
     -      eps**3*Pi**4*Cos((6*car(1)*Pi)/xmax - (2*car(2)*Pi)/ymax) - 
     -      2*eps**3*Pi**4*Cos((2*car(1)*Pi)/xmax + (2*car(2)*Pi)/ymax)- 
     -      eps*Pi**2*ymax**2*
     -       Cos((2*car(1)*Pi)/xmax + (2*car(2)*Pi)/ymax) - 
     -      eps**3*Pi**4*Cos((6*car(1)*Pi)/xmax + (2*car(2)*Pi)/ymax) - 
     -      eps**3*Pi**4*Cos((2*car(1)*Pi)/xmax + (6*car(2)*Pi)/ymax) + 
     -      4*eps**2*Pi**3*ymax*Sin((4*car(2)*Pi)/ymax)))/
     -  (-2*eps**2*Pi**2*xmax**2 - 2*eps**2*Pi**2*ymax**2 - 
     -    2*xmax**2*ymax**2 + 
     -    2*eps**2*Pi**2*xmax**2*Cos((4*car(1)*Pi)/xmax) - 
     -    2*eps**2*Pi**2*ymax**2*Cos((4*car(1)*Pi)/xmax) + 
     -    eps**2*Pi**2*xmax**2*
     -     Cos((4*car(1)*Pi)/xmax - (4*car(2)*Pi)/ymax) - 
     -    2*eps**2*Pi**2*xmax*ymax*
     -     Cos((4*car(1)*Pi)/xmax - (4*car(2)*Pi)/ymax) + 
     -    eps**2*Pi**2*ymax**2*
     -     Cos((4*car(1)*Pi)/xmax - (4*car(2)*Pi)/ymax) + 
     -    eps**2*Pi**2*xmax**2*
     -     Cos((4*car(1)*Pi)/xmax + (4*car(2)*Pi)/ymax) + 
     -    2*eps**2*Pi**2*xmax*ymax*
     -     Cos((4*car(1)*Pi)/xmax + (4*car(2)*Pi)/ymax) + 
     -    eps**2*Pi**2*ymax**2*
     -     Cos((4*car(1)*Pi)/xmax + (4*car(2)*Pi)/ymax) - 
     -    2*eps**2*Pi**2*xmax**2*Cos((4*car(2)*Pi)/ymax) + 
     -    2*eps**2*Pi**2*ymax**2*Cos((4*car(2)*Pi)/ymax) - 
     -    4*eps*Pi*xmax**2*ymax*
     -     Sin((2*car(1)*Pi)/xmax - (2*car(2)*Pi)/ymax) + 
     -    4*eps*Pi*xmax*ymax**2*
     -     Sin((2*car(1)*Pi)/xmax - (2*car(2)*Pi)/ymax) - 
     -    4*eps*Pi*xmax**2*ymax*
     -     Sin((2*car(1)*Pi)/xmax + (2*car(2)*Pi)/ymax) - 
     -    4*eps*Pi*xmax*ymax**2*
     -     Sin((2*car(1)*Pi)/xmax + (2*car(2)*Pi)/ymax))

              vec(2) =
     .   (-4*(eps**3*Pi**4*Cos((2*car(1)*Pi)/xmax - (6*car(2)*Pi)/ymax)+ 
     -      2*eps**3*Pi**4*Cos((2*car(1)*Pi)/xmax - (2*car(2)*Pi)/ymax)+ 
     -      eps*Pi**2*xmax*ymax*
     -       Cos((2*car(1)*Pi)/xmax - (2*car(2)*Pi)/ymax) + 
     -      eps**3*Pi**4*Cos((6*car(1)*Pi)/xmax - (2*car(2)*Pi)/ymax) - 
     -      2*eps**3*Pi**4*Cos((2*car(1)*Pi)/xmax + (2*car(2)*Pi)/ymax)+ 
     -      eps*Pi**2*xmax*ymax*
     -       Cos((2*car(1)*Pi)/xmax + (2*car(2)*Pi)/ymax) - 
     -      eps**3*Pi**4*Cos((6*car(1)*Pi)/xmax + (2*car(2)*Pi)/ymax) - 
     -      eps**3*Pi**4*Cos((2*car(1)*Pi)/xmax + (6*car(2)*Pi)/ymax) + 
     -      2*eps**2*Pi**3*xmax*Sin((4*car(1)*Pi)/xmax) + 
     -      2*eps**2*Pi**3*ymax*Sin((4*car(2)*Pi)/ymax)))/
     -  (-2*eps**2*Pi**2*xmax**2 - 2*eps**2*Pi**2*ymax**2 - 
     -    2*xmax**2*ymax**2 + 
     -    2*eps**2*Pi**2*xmax**2*Cos((4*car(1)*Pi)/xmax) - 
     -    2*eps**2*Pi**2*ymax**2*Cos((4*car(1)*Pi)/xmax) + 
     -    eps**2*Pi**2*xmax**2*
     -     Cos((4*car(1)*Pi)/xmax - (4*car(2)*Pi)/ymax) - 
     -    2*eps**2*Pi**2*xmax*ymax*
     -     Cos((4*car(1)*Pi)/xmax - (4*car(2)*Pi)/ymax) + 
     -    eps**2*Pi**2*ymax**2*
     -     Cos((4*car(1)*Pi)/xmax - (4*car(2)*Pi)/ymax) + 
     -    eps**2*Pi**2*xmax**2*
     -     Cos((4*car(1)*Pi)/xmax + (4*car(2)*Pi)/ymax) + 
     -    2*eps**2*Pi**2*xmax*ymax*
     -     Cos((4*car(1)*Pi)/xmax + (4*car(2)*Pi)/ymax) + 
     -    eps**2*Pi**2*ymax**2*
     -     Cos((4*car(1)*Pi)/xmax + (4*car(2)*Pi)/ymax) - 
     -    2*eps**2*Pi**2*xmax**2*Cos((4*car(2)*Pi)/ymax) + 
     -    2*eps**2*Pi**2*ymax**2*Cos((4*car(2)*Pi)/ymax) - 
     -    4*eps*Pi*xmax**2*ymax*
     -     Sin((2*car(1)*Pi)/xmax - (2*car(2)*Pi)/ymax) + 
     -    4*eps*Pi*xmax*ymax**2*
     -     Sin((2*car(1)*Pi)/xmax - (2*car(2)*Pi)/ymax) - 
     -    4*eps*Pi*xmax**2*ymax*
     -     Sin((2*car(1)*Pi)/xmax + (2*car(2)*Pi)/ymax) - 
     -    4*eps*Pi*xmax*ymax**2*
     -     Sin((2*car(1)*Pi)/xmax + (2*car(2)*Pi)/ymax))

              vec(3) = 0d0
              vec(4) = vec(2)

              vec(5) =
     .   (4*(eps**3*Pi**4*Cos((2*car(1)*Pi)/xmax - (6*car(2)*Pi)/ymax) + 
     -      2*eps**3*Pi**4*Cos((2*car(1)*Pi)/xmax - (2*car(2)*Pi)/ymax)+ 
     -      eps*Pi**2*xmax**2*
     -       Cos((2*car(1)*Pi)/xmax - (2*car(2)*Pi)/ymax) + 
     -      eps**3*Pi**4*Cos((6*car(1)*Pi)/xmax - (2*car(2)*Pi)/ymax) - 
     -      2*eps**3*Pi**4*Cos((2*car(1)*Pi)/xmax + (2*car(2)*Pi)/ymax)- 
     -      eps*Pi**2*xmax**2*
     -       Cos((2*car(1)*Pi)/xmax + (2*car(2)*Pi)/ymax) - 
     -      eps**3*Pi**4*Cos((6*car(1)*Pi)/xmax + (2*car(2)*Pi)/ymax) - 
     -      eps**3*Pi**4*Cos((2*car(1)*Pi)/xmax + (6*car(2)*Pi)/ymax) + 
     -      4*eps**2*Pi**3*xmax*Sin((4*car(1)*Pi)/xmax)))/
     -  (-2*eps**2*Pi**2*xmax**2 - 2*eps**2*Pi**2*ymax**2 - 
     -    2*xmax**2*ymax**2 + 
     -    2*eps**2*Pi**2*xmax**2*Cos((4*car(1)*Pi)/xmax) - 
     -    2*eps**2*Pi**2*ymax**2*Cos((4*car(1)*Pi)/xmax) + 
     -    eps**2*Pi**2*xmax**2*
     -     Cos((4*car(1)*Pi)/xmax - (4*car(2)*Pi)/ymax) - 
     -    2*eps**2*Pi**2*xmax*ymax*
     -     Cos((4*car(1)*Pi)/xmax - (4*car(2)*Pi)/ymax) + 
     -    eps**2*Pi**2*ymax**2*
     -     Cos((4*car(1)*Pi)/xmax - (4*car(2)*Pi)/ymax) + 
     -    eps**2*Pi**2*xmax**2*
     -     Cos((4*car(1)*Pi)/xmax + (4*car(2)*Pi)/ymax) + 
     -    2*eps**2*Pi**2*xmax*ymax*
     -     Cos((4*car(1)*Pi)/xmax + (4*car(2)*Pi)/ymax) + 
     -    eps**2*Pi**2*ymax**2*
     -     Cos((4*car(1)*Pi)/xmax + (4*car(2)*Pi)/ymax) - 
     -    2*eps**2*Pi**2*xmax**2*Cos((4*car(2)*Pi)/ymax) + 
     -    2*eps**2*Pi**2*ymax**2*Cos((4*car(2)*Pi)/ymax) - 
     -    4*eps*Pi*xmax**2*ymax*
     -     Sin((2*car(1)*Pi)/xmax - (2*car(2)*Pi)/ymax) + 
     -    4*eps*Pi*xmax*ymax**2*
     -     Sin((2*car(1)*Pi)/xmax - (2*car(2)*Pi)/ymax) - 
     -    4*eps*Pi*xmax**2*ymax*
     -     Sin((2*car(1)*Pi)/xmax + (2*car(2)*Pi)/ymax) - 
     -    4*eps*Pi*xmax*ymax**2*
     -     Sin((2*car(1)*Pi)/xmax + (2*car(2)*Pi)/ymax))

              vec(6) = 0d0
              vec(7) = vec(3)
              vec(8) = vec(6)
              vec(9) = 0d0
            case (2)
              vec(1) =
     .    (4*(eps**3*Pi**4*Cos((2*car(1)*Pi)/xmax - (6*car(2)*Pi)/ymax)+ 
     -      2*eps**3*Pi**4*Cos((2*car(1)*Pi)/xmax - (2*car(2)*Pi)/ymax)+ 
     -      eps*Pi**2*ymax**2*
     -       Cos((2*car(1)*Pi)/xmax - (2*car(2)*Pi)/ymax) + 
     -      eps**3*Pi**4*Cos((6*car(1)*Pi)/xmax - (2*car(2)*Pi)/ymax) - 
     -      2*eps**3*Pi**4*Cos((2*car(1)*Pi)/xmax + (2*car(2)*Pi)/ymax)- 
     -      eps*Pi**2*ymax**2*
     -       Cos((2*car(1)*Pi)/xmax + (2*car(2)*Pi)/ymax) - 
     -      eps**3*Pi**4*Cos((6*car(1)*Pi)/xmax + (2*car(2)*Pi)/ymax) - 
     -      eps**3*Pi**4*Cos((2*car(1)*Pi)/xmax + (6*car(2)*Pi)/ymax) + 
     -      4*eps**2*Pi**3*ymax*Sin((4*car(2)*Pi)/ymax)))/
     -  (-2*eps**2*Pi**2*xmax**2 - 2*eps**2*Pi**2*ymax**2 - 
     -    2*xmax**2*ymax**2 + 
     -    2*eps**2*Pi**2*xmax**2*Cos((4*car(1)*Pi)/xmax) - 
     -    2*eps**2*Pi**2*ymax**2*Cos((4*car(1)*Pi)/xmax) + 
     -    eps**2*Pi**2*xmax**2*
     -     Cos((4*car(1)*Pi)/xmax - (4*car(2)*Pi)/ymax) - 
     -    2*eps**2*Pi**2*xmax*ymax*
     -     Cos((4*car(1)*Pi)/xmax - (4*car(2)*Pi)/ymax) + 
     -    eps**2*Pi**2*ymax**2*
     -     Cos((4*car(1)*Pi)/xmax - (4*car(2)*Pi)/ymax) + 
     -    eps**2*Pi**2*xmax**2*
     -     Cos((4*car(1)*Pi)/xmax + (4*car(2)*Pi)/ymax) + 
     -    2*eps**2*Pi**2*xmax*ymax*
     -     Cos((4*car(1)*Pi)/xmax + (4*car(2)*Pi)/ymax) + 
     -    eps**2*Pi**2*ymax**2*
     -     Cos((4*car(1)*Pi)/xmax + (4*car(2)*Pi)/ymax) - 
     -    2*eps**2*Pi**2*xmax**2*Cos((4*car(2)*Pi)/ymax) + 
     -    2*eps**2*Pi**2*ymax**2*Cos((4*car(2)*Pi)/ymax) - 
     -    4*eps*Pi*xmax**2*ymax*
     -     Sin((2*car(1)*Pi)/xmax - (2*car(2)*Pi)/ymax) + 
     -    4*eps*Pi*xmax*ymax**2*
     -     Sin((2*car(1)*Pi)/xmax - (2*car(2)*Pi)/ymax) - 
     -    4*eps*Pi*xmax**2*ymax*
     -     Sin((2*car(1)*Pi)/xmax + (2*car(2)*Pi)/ymax) - 
     -    4*eps*Pi*xmax*ymax**2*
     -     Sin((2*car(1)*Pi)/xmax + (2*car(2)*Pi)/ymax))

              vec(2) =
     .   (-4*(eps**3*Pi**4*Cos((2*car(1)*Pi)/xmax - (6*car(2)*Pi)/ymax)+ 
     -      2*eps**3*Pi**4*Cos((2*car(1)*Pi)/xmax - (2*car(2)*Pi)/ymax)+ 
     -      eps*Pi**2*xmax*ymax*
     -       Cos((2*car(1)*Pi)/xmax - (2*car(2)*Pi)/ymax) + 
     -      eps**3*Pi**4*Cos((6*car(1)*Pi)/xmax - (2*car(2)*Pi)/ymax) - 
     -      2*eps**3*Pi**4*Cos((2*car(1)*Pi)/xmax + (2*car(2)*Pi)/ymax)+ 
     -      eps*Pi**2*xmax*ymax*
     -       Cos((2*car(1)*Pi)/xmax + (2*car(2)*Pi)/ymax) - 
     -      eps**3*Pi**4*Cos((6*car(1)*Pi)/xmax + (2*car(2)*Pi)/ymax) - 
     -      eps**3*Pi**4*Cos((2*car(1)*Pi)/xmax + (6*car(2)*Pi)/ymax) + 
     -      2*eps**2*Pi**3*xmax*Sin((4*car(1)*Pi)/xmax) + 
     -      2*eps**2*Pi**3*ymax*Sin((4*car(2)*Pi)/ymax)))/
     -  (-2*eps**2*Pi**2*xmax**2 - 2*eps**2*Pi**2*ymax**2 - 
     -    2*xmax**2*ymax**2 + 
     -    2*eps**2*Pi**2*xmax**2*Cos((4*car(1)*Pi)/xmax) - 
     -    2*eps**2*Pi**2*ymax**2*Cos((4*car(1)*Pi)/xmax) + 
     -    eps**2*Pi**2*xmax**2*
     -     Cos((4*car(1)*Pi)/xmax - (4*car(2)*Pi)/ymax) - 
     -    2*eps**2*Pi**2*xmax*ymax*
     -     Cos((4*car(1)*Pi)/xmax - (4*car(2)*Pi)/ymax) + 
     -    eps**2*Pi**2*ymax**2*
     -     Cos((4*car(1)*Pi)/xmax - (4*car(2)*Pi)/ymax) + 
     -    eps**2*Pi**2*xmax**2*
     -     Cos((4*car(1)*Pi)/xmax + (4*car(2)*Pi)/ymax) + 
     -    2*eps**2*Pi**2*xmax*ymax*
     -     Cos((4*car(1)*Pi)/xmax + (4*car(2)*Pi)/ymax) + 
     -    eps**2*Pi**2*ymax**2*
     -     Cos((4*car(1)*Pi)/xmax + (4*car(2)*Pi)/ymax) - 
     -    2*eps**2*Pi**2*xmax**2*Cos((4*car(2)*Pi)/ymax) + 
     -    2*eps**2*Pi**2*ymax**2*Cos((4*car(2)*Pi)/ymax) - 
     -    4*eps*Pi*xmax**2*ymax*
     -     Sin((2*car(1)*Pi)/xmax - (2*car(2)*Pi)/ymax) + 
     -    4*eps*Pi*xmax*ymax**2*
     -     Sin((2*car(1)*Pi)/xmax - (2*car(2)*Pi)/ymax) - 
     -    4*eps*Pi*xmax**2*ymax*
     -     Sin((2*car(1)*Pi)/xmax + (2*car(2)*Pi)/ymax) - 
     -    4*eps*Pi*xmax*ymax**2*
     -     Sin((2*car(1)*Pi)/xmax + (2*car(2)*Pi)/ymax))

              vec(3) = 0d0
              vec(4) = vec(2)

              vec(5) =
     .    (4*(eps**3*Pi**4*Cos((2*car(1)*Pi)/xmax - (6*car(2)*Pi)/ymax)+ 
     -      2*eps**3*Pi**4*Cos((2*car(1)*Pi)/xmax - (2*car(2)*Pi)/ymax)+ 
     -      eps*Pi**2*xmax**2*
     -       Cos((2*car(1)*Pi)/xmax - (2*car(2)*Pi)/ymax) + 
     -      eps**3*Pi**4*Cos((6*car(1)*Pi)/xmax - (2*car(2)*Pi)/ymax) - 
     -      2*eps**3*Pi**4*Cos((2*car(1)*Pi)/xmax + (2*car(2)*Pi)/ymax)- 
     -      eps*Pi**2*xmax**2*
     -       Cos((2*car(1)*Pi)/xmax + (2*car(2)*Pi)/ymax) - 
     -      eps**3*Pi**4*Cos((6*car(1)*Pi)/xmax + (2*car(2)*Pi)/ymax) - 
     -      eps**3*Pi**4*Cos((2*car(1)*Pi)/xmax + (6*car(2)*Pi)/ymax) + 
     -      4*eps**2*Pi**3*xmax*Sin((4*car(1)*Pi)/xmax)))/
     -  (-2*eps**2*Pi**2*xmax**2 - 2*eps**2*Pi**2*ymax**2 - 
     -    2*xmax**2*ymax**2 + 
     -    2*eps**2*Pi**2*xmax**2*Cos((4*car(1)*Pi)/xmax) - 
     -    2*eps**2*Pi**2*ymax**2*Cos((4*car(1)*Pi)/xmax) + 
     -    eps**2*Pi**2*xmax**2*
     -     Cos((4*car(1)*Pi)/xmax - (4*car(2)*Pi)/ymax) - 
     -    2*eps**2*Pi**2*xmax*ymax*
     -     Cos((4*car(1)*Pi)/xmax - (4*car(2)*Pi)/ymax) + 
     -    eps**2*Pi**2*ymax**2*
     -     Cos((4*car(1)*Pi)/xmax - (4*car(2)*Pi)/ymax) + 
     -    eps**2*Pi**2*xmax**2*
     -     Cos((4*car(1)*Pi)/xmax + (4*car(2)*Pi)/ymax) + 
     -    2*eps**2*Pi**2*xmax*ymax*
     -     Cos((4*car(1)*Pi)/xmax + (4*car(2)*Pi)/ymax) + 
     -    eps**2*Pi**2*ymax**2*
     -     Cos((4*car(1)*Pi)/xmax + (4*car(2)*Pi)/ymax) - 
     -    2*eps**2*Pi**2*xmax**2*Cos((4*car(2)*Pi)/ymax) + 
     -    2*eps**2*Pi**2*ymax**2*Cos((4*car(2)*Pi)/ymax) - 
     -    4*eps*Pi*xmax**2*ymax*
     -     Sin((2*car(1)*Pi)/xmax - (2*car(2)*Pi)/ymax) + 
     -    4*eps*Pi*xmax*ymax**2*
     -     Sin((2*car(1)*Pi)/xmax - (2*car(2)*Pi)/ymax) - 
     -    4*eps*Pi*xmax**2*ymax*
     -     Sin((2*car(1)*Pi)/xmax + (2*car(2)*Pi)/ymax) - 
     -    4*eps*Pi*xmax*ymax**2*
     -     Sin((2*car(1)*Pi)/xmax + (2*car(2)*Pi)/ymax))

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
          write (*,*) 'Grid not implemented'
          write (*,*) 'Aborting...'
          stop
        end select

        tensor = transpose(reshape(vec, (/3,3/)))

      end function hessian

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
        case ('sin')

          pi = acos(-1d0)

          car = x_xi(x1,x2,x3)
          eps = gparams(1)
          select case (k)
            case (1)
              vec(1) =
     .   (-4*(eps**3*Pi**4*Cos((2*car(1)*Pi)/xmax - (6*car(2)*Pi)/ymax)+ 
     -      2*eps**3*Pi**4*Cos((2*car(1)*Pi)/xmax - (2*car(2)*Pi)/ymax)+ 
     -      eps*Pi**2*xmax*ymax*
     -       Cos((2*car(1)*Pi)/xmax - (2*car(2)*Pi)/ymax) + 
     -      eps**3*Pi**4*Cos((6*car(1)*Pi)/xmax - (2*car(2)*Pi)/ymax) - 
     -      2*eps**3*Pi**4*Cos((2*car(1)*Pi)/xmax + (2*car(2)*Pi)/ymax)+ 
     -      eps*Pi**2*xmax*ymax*
     -       Cos((2*car(1)*Pi)/xmax + (2*car(2)*Pi)/ymax) - 
     -      eps**3*Pi**4*Cos((6*car(1)*Pi)/xmax + (2*car(2)*Pi)/ymax) - 
     -      eps**3*Pi**4*Cos((2*car(1)*Pi)/xmax + (6*car(2)*Pi)/ymax) + 
     -      2*eps**2*Pi**3*xmax*Sin((4*car(1)*Pi)/xmax) + 
     -      2*eps**2*Pi**3*ymax*Sin((4*car(2)*Pi)/ymax)))/
     -  (-2*eps**2*Pi**2*xmax**2 - 2*eps**2*Pi**2*ymax**2 - 
     -    2*xmax**2*ymax**2 + 
     -    2*eps**2*Pi**2*xmax**2*Cos((4*car(1)*Pi)/xmax) - 
     -    2*eps**2*Pi**2*ymax**2*Cos((4*car(1)*Pi)/xmax) + 
     -    eps**2*Pi**2*xmax**2*
     -     Cos((4*car(1)*Pi)/xmax - (4*car(2)*Pi)/ymax) - 
     -    2*eps**2*Pi**2*xmax*ymax*
     -     Cos((4*car(1)*Pi)/xmax - (4*car(2)*Pi)/ymax) + 
     -    eps**2*Pi**2*ymax**2*
     -     Cos((4*car(1)*Pi)/xmax - (4*car(2)*Pi)/ymax) + 
     -    eps**2*Pi**2*xmax**2*
     -     Cos((4*car(1)*Pi)/xmax + (4*car(2)*Pi)/ymax) + 
     -    2*eps**2*Pi**2*xmax*ymax*
     -     Cos((4*car(1)*Pi)/xmax + (4*car(2)*Pi)/ymax) + 
     -    eps**2*Pi**2*ymax**2*
     -     Cos((4*car(1)*Pi)/xmax + (4*car(2)*Pi)/ymax) - 
     -    2*eps**2*Pi**2*xmax**2*Cos((4*car(2)*Pi)/ymax) + 
     -    2*eps**2*Pi**2*ymax**2*Cos((4*car(2)*Pi)/ymax) - 
     -    4*eps*Pi*xmax**2*ymax*
     -     Sin((2*car(1)*Pi)/xmax - (2*car(2)*Pi)/ymax) + 
     -    4*eps*Pi*xmax*ymax**2*
     -     Sin((2*car(1)*Pi)/xmax - (2*car(2)*Pi)/ymax) - 
     -    4*eps*Pi*xmax**2*ymax*
     -     Sin((2*car(1)*Pi)/xmax + (2*car(2)*Pi)/ymax) - 
     -    4*eps*Pi*xmax*ymax**2*
     -     Sin((2*car(1)*Pi)/xmax + (2*car(2)*Pi)/ymax))
              vec(2) =
     .   (-4*(eps**3*Pi**4*Cos((2*car(1)*Pi)/xmax - (6*car(2)*Pi)/ymax)+ 
     -      2*eps**3*Pi**4*Cos((2*car(1)*Pi)/xmax - (2*car(2)*Pi)/ymax)+ 
     -      eps*Pi**2*ymax**2*
     -       Cos((2*car(1)*Pi)/xmax - (2*car(2)*Pi)/ymax) + 
     -      eps**3*Pi**4*Cos((6*car(1)*Pi)/xmax - (2*car(2)*Pi)/ymax) - 
     -      2*eps**3*Pi**4*Cos((2*car(1)*Pi)/xmax + (2*car(2)*Pi)/ymax)- 
     -      eps*Pi**2*ymax**2*
     -       Cos((2*car(1)*Pi)/xmax + (2*car(2)*Pi)/ymax) - 
     -      eps**3*Pi**4*Cos((6*car(1)*Pi)/xmax + (2*car(2)*Pi)/ymax) - 
     -      eps**3*Pi**4*Cos((2*car(1)*Pi)/xmax + (6*car(2)*Pi)/ymax) + 
     -      4*eps**2*Pi**3*ymax*Sin((4*car(2)*Pi)/ymax)))/
     -  (-2*eps**2*Pi**2*xmax**2 - 2*eps**2*Pi**2*ymax**2 - 
     -    2*xmax**2*ymax**2 + 
     -    2*eps**2*Pi**2*xmax**2*Cos((4*car(1)*Pi)/xmax) - 
     -    2*eps**2*Pi**2*ymax**2*Cos((4*car(1)*Pi)/xmax) + 
     -    eps**2*Pi**2*xmax**2*
     -     Cos((4*car(1)*Pi)/xmax - (4*car(2)*Pi)/ymax) - 
     -    2*eps**2*Pi**2*xmax*ymax*
     -     Cos((4*car(1)*Pi)/xmax - (4*car(2)*Pi)/ymax) + 
     -    eps**2*Pi**2*ymax**2*
     -     Cos((4*car(1)*Pi)/xmax - (4*car(2)*Pi)/ymax) + 
     -    eps**2*Pi**2*xmax**2*
     -     Cos((4*car(1)*Pi)/xmax + (4*car(2)*Pi)/ymax) + 
     -    2*eps**2*Pi**2*xmax*ymax*
     -     Cos((4*car(1)*Pi)/xmax + (4*car(2)*Pi)/ymax) + 
     -    eps**2*Pi**2*ymax**2*
     -     Cos((4*car(1)*Pi)/xmax + (4*car(2)*Pi)/ymax) - 
     -    2*eps**2*Pi**2*xmax**2*Cos((4*car(2)*Pi)/ymax) + 
     -    2*eps**2*Pi**2*ymax**2*Cos((4*car(2)*Pi)/ymax) - 
     -    4*eps*Pi*xmax**2*ymax*
     -     Sin((2*car(1)*Pi)/xmax - (2*car(2)*Pi)/ymax) + 
     -    4*eps*Pi*xmax*ymax**2*
     -     Sin((2*car(1)*Pi)/xmax - (2*car(2)*Pi)/ymax) - 
     -    4*eps*Pi*xmax**2*ymax*
     -     Sin((2*car(1)*Pi)/xmax + (2*car(2)*Pi)/ymax) - 
     -    4*eps*Pi*xmax*ymax**2*
     -     Sin((2*car(1)*Pi)/xmax + (2*car(2)*Pi)/ymax))
              vec(3) = 0d0
              vec(4) =
     .    (4*(eps**3*Pi**4*Cos((2*car(1)*Pi)/xmax - (6*car(2)*Pi)/ymax)+ 
     -      2*eps**3*Pi**4*Cos((2*car(1)*Pi)/xmax - (2*car(2)*Pi)/ymax)+ 
     -      eps*Pi**2*xmax**2*
     -       Cos((2*car(1)*Pi)/xmax - (2*car(2)*Pi)/ymax) + 
     -      eps**3*Pi**4*Cos((6*car(1)*Pi)/xmax - (2*car(2)*Pi)/ymax) - 
     -      2*eps**3*Pi**4*Cos((2*car(1)*Pi)/xmax + (2*car(2)*Pi)/ymax)- 
     -      eps*Pi**2*xmax**2*
     -       Cos((2*car(1)*Pi)/xmax + (2*car(2)*Pi)/ymax) - 
     -      eps**3*Pi**4*Cos((6*car(1)*Pi)/xmax + (2*car(2)*Pi)/ymax) - 
     -      eps**3*Pi**4*Cos((2*car(1)*Pi)/xmax + (6*car(2)*Pi)/ymax) + 
     -      4*eps**2*Pi**3*xmax*Sin((4*car(1)*Pi)/xmax)))/
     -  (-2*eps**2*Pi**2*xmax**2 - 2*eps**2*Pi**2*ymax**2 - 
     -    2*xmax**2*ymax**2 + 
     -    2*eps**2*Pi**2*xmax**2*Cos((4*car(1)*Pi)/xmax) - 
     -    2*eps**2*Pi**2*ymax**2*Cos((4*car(1)*Pi)/xmax) + 
     -    eps**2*Pi**2*xmax**2*
     -     Cos((4*car(1)*Pi)/xmax - (4*car(2)*Pi)/ymax) - 
     -    2*eps**2*Pi**2*xmax*ymax*
     -     Cos((4*car(1)*Pi)/xmax - (4*car(2)*Pi)/ymax) + 
     -    eps**2*Pi**2*ymax**2*
     -     Cos((4*car(1)*Pi)/xmax - (4*car(2)*Pi)/ymax) + 
     -    eps**2*Pi**2*xmax**2*
     -     Cos((4*car(1)*Pi)/xmax + (4*car(2)*Pi)/ymax) + 
     -    2*eps**2*Pi**2*xmax*ymax*
     -     Cos((4*car(1)*Pi)/xmax + (4*car(2)*Pi)/ymax) + 
     -    eps**2*Pi**2*ymax**2*
     -     Cos((4*car(1)*Pi)/xmax + (4*car(2)*Pi)/ymax) - 
     -    2*eps**2*Pi**2*xmax**2*Cos((4*car(2)*Pi)/ymax) + 
     -    2*eps**2*Pi**2*ymax**2*Cos((4*car(2)*Pi)/ymax) - 
     -    4*eps*Pi*xmax**2*ymax*
     -     Sin((2*car(1)*Pi)/xmax - (2*car(2)*Pi)/ymax) + 
     -    4*eps*Pi*xmax*ymax**2*
     -     Sin((2*car(1)*Pi)/xmax - (2*car(2)*Pi)/ymax) - 
     -    4*eps*Pi*xmax**2*ymax*
     -     Sin((2*car(1)*Pi)/xmax + (2*car(2)*Pi)/ymax) - 
     -    4*eps*Pi*xmax*ymax**2*
     -     Sin((2*car(1)*Pi)/xmax + (2*car(2)*Pi)/ymax))
              vec(5) =
     .    (4*(eps**3*Pi**4*Cos((2*car(1)*Pi)/xmax - (6*car(2)*Pi)/ymax)+ 
     -      2*eps**3*Pi**4*Cos((2*car(1)*Pi)/xmax - (2*car(2)*Pi)/ymax)+ 
     -      eps*Pi**2*xmax*ymax*
     -       Cos((2*car(1)*Pi)/xmax - (2*car(2)*Pi)/ymax) + 
     -      eps**3*Pi**4*Cos((6*car(1)*Pi)/xmax - (2*car(2)*Pi)/ymax) - 
     -      2*eps**3*Pi**4*Cos((2*car(1)*Pi)/xmax + (2*car(2)*Pi)/ymax)+ 
     -      eps*Pi**2*xmax*ymax*
     -       Cos((2*car(1)*Pi)/xmax + (2*car(2)*Pi)/ymax) - 
     -      eps**3*Pi**4*Cos((6*car(1)*Pi)/xmax + (2*car(2)*Pi)/ymax) - 
     -      eps**3*Pi**4*Cos((2*car(1)*Pi)/xmax + (6*car(2)*Pi)/ymax) + 
     -      2*eps**2*Pi**3*xmax*Sin((4*car(1)*Pi)/xmax) + 
     -      2*eps**2*Pi**3*ymax*Sin((4*car(2)*Pi)/ymax)))/
     -  (-2*eps**2*Pi**2*xmax**2 - 2*eps**2*Pi**2*ymax**2 - 
     -    2*xmax**2*ymax**2 + 
     -    2*eps**2*Pi**2*xmax**2*Cos((4*car(1)*Pi)/xmax) - 
     -    2*eps**2*Pi**2*ymax**2*Cos((4*car(1)*Pi)/xmax) + 
     -    eps**2*Pi**2*xmax**2*
     -     Cos((4*car(1)*Pi)/xmax - (4*car(2)*Pi)/ymax) - 
     -    2*eps**2*Pi**2*xmax*ymax*
     -     Cos((4*car(1)*Pi)/xmax - (4*car(2)*Pi)/ymax) + 
     -    eps**2*Pi**2*ymax**2*
     -     Cos((4*car(1)*Pi)/xmax - (4*car(2)*Pi)/ymax) + 
     -    eps**2*Pi**2*xmax**2*
     -     Cos((4*car(1)*Pi)/xmax + (4*car(2)*Pi)/ymax) + 
     -    2*eps**2*Pi**2*xmax*ymax*
     -     Cos((4*car(1)*Pi)/xmax + (4*car(2)*Pi)/ymax) + 
     -    eps**2*Pi**2*ymax**2*
     -     Cos((4*car(1)*Pi)/xmax + (4*car(2)*Pi)/ymax) - 
     -    2*eps**2*Pi**2*xmax**2*Cos((4*car(2)*Pi)/ymax) + 
     -    2*eps**2*Pi**2*ymax**2*Cos((4*car(2)*Pi)/ymax) - 
     -    4*eps*Pi*xmax**2*ymax*
     -     Sin((2*car(1)*Pi)/xmax - (2*car(2)*Pi)/ymax) + 
     -    4*eps*Pi*xmax*ymax**2*
     -     Sin((2*car(1)*Pi)/xmax - (2*car(2)*Pi)/ymax) - 
     -    4*eps*Pi*xmax**2*ymax*
     -     Sin((2*car(1)*Pi)/xmax + (2*car(2)*Pi)/ymax) - 
     -    4*eps*Pi*xmax*ymax**2*
     -     Sin((2*car(1)*Pi)/xmax + (2*car(2)*Pi)/ymax))
              vec(6) = 0d0
              vec(7) = 0d0
              vec(8) = 0d0
              vec(9) = 0d0
            case (2)
              vec(1) =
     .    (4*(eps**3*Pi**4*Cos((2*car(1)*Pi)/xmax - (6*car(2)*Pi)/ymax)+ 
     -      2*eps**3*Pi**4*Cos((2*car(1)*Pi)/xmax - (2*car(2)*Pi)/ymax)+ 
     -      eps*Pi**2*xmax*ymax*
     -       Cos((2*car(1)*Pi)/xmax - (2*car(2)*Pi)/ymax) + 
     -      eps**3*Pi**4*Cos((6*car(1)*Pi)/xmax - (2*car(2)*Pi)/ymax) - 
     -      2*eps**3*Pi**4*Cos((2*car(1)*Pi)/xmax + (2*car(2)*Pi)/ymax)+ 
     -      eps*Pi**2*xmax*ymax*
     -       Cos((2*car(1)*Pi)/xmax + (2*car(2)*Pi)/ymax) - 
     -      eps**3*Pi**4*Cos((6*car(1)*Pi)/xmax + (2*car(2)*Pi)/ymax) - 
     -      eps**3*Pi**4*Cos((2*car(1)*Pi)/xmax + (6*car(2)*Pi)/ymax) + 
     -      2*eps**2*Pi**3*xmax*Sin((4*car(1)*Pi)/xmax) + 
     -      2*eps**2*Pi**3*ymax*Sin((4*car(2)*Pi)/ymax)))/
     -  (-2*eps**2*Pi**2*xmax**2 - 2*eps**2*Pi**2*ymax**2 - 
     -    2*xmax**2*ymax**2 + 
     -    2*eps**2*Pi**2*xmax**2*Cos((4*car(1)*Pi)/xmax) - 
     -    2*eps**2*Pi**2*ymax**2*Cos((4*car(1)*Pi)/xmax) + 
     -    eps**2*Pi**2*xmax**2*
     -     Cos((4*car(1)*Pi)/xmax - (4*car(2)*Pi)/ymax) - 
     -    2*eps**2*Pi**2*xmax*ymax*
     -     Cos((4*car(1)*Pi)/xmax - (4*car(2)*Pi)/ymax) + 
     -    eps**2*Pi**2*ymax**2*
     -     Cos((4*car(1)*Pi)/xmax - (4*car(2)*Pi)/ymax) + 
     -    eps**2*Pi**2*xmax**2*
     -     Cos((4*car(1)*Pi)/xmax + (4*car(2)*Pi)/ymax) + 
     -    2*eps**2*Pi**2*xmax*ymax*
     -     Cos((4*car(1)*Pi)/xmax + (4*car(2)*Pi)/ymax) + 
     -    eps**2*Pi**2*ymax**2*
     -     Cos((4*car(1)*Pi)/xmax + (4*car(2)*Pi)/ymax) - 
     -    2*eps**2*Pi**2*xmax**2*Cos((4*car(2)*Pi)/ymax) + 
     -    2*eps**2*Pi**2*ymax**2*Cos((4*car(2)*Pi)/ymax) - 
     -    4*eps*Pi*xmax**2*ymax*
     -     Sin((2*car(1)*Pi)/xmax - (2*car(2)*Pi)/ymax) + 
     -    4*eps*Pi*xmax*ymax**2*
     -     Sin((2*car(1)*Pi)/xmax - (2*car(2)*Pi)/ymax) - 
     -    4*eps*Pi*xmax**2*ymax*
     -     Sin((2*car(1)*Pi)/xmax + (2*car(2)*Pi)/ymax) - 
     -    4*eps*Pi*xmax*ymax**2*
     -     Sin((2*car(1)*Pi)/xmax + (2*car(2)*Pi)/ymax))
              vec(2) =
     .   (4*(eps**3*Pi**4*Cos((2*car(1)*Pi)/xmax - (6*car(2)*Pi)/ymax) + 
     -     2*eps**3*Pi**4*Cos((2*car(1)*Pi)/xmax - (2*car(2)*Pi)/ymax) + 
     -      eps*Pi**2*ymax**2*
     -       Cos((2*car(1)*Pi)/xmax - (2*car(2)*Pi)/ymax) + 
     -      eps**3*Pi**4*Cos((6*car(1)*Pi)/xmax - (2*car(2)*Pi)/ymax) - 
     -     2*eps**3*Pi**4*Cos((2*car(1)*Pi)/xmax + (2*car(2)*Pi)/ymax) - 
     -      eps*Pi**2*ymax**2*
     -       Cos((2*car(1)*Pi)/xmax + (2*car(2)*Pi)/ymax) - 
     -      eps**3*Pi**4*Cos((6*car(1)*Pi)/xmax + (2*car(2)*Pi)/ymax) - 
     -      eps**3*Pi**4*Cos((2*car(1)*Pi)/xmax + (6*car(2)*Pi)/ymax) + 
     -      4*eps**2*Pi**3*ymax*Sin((4*car(2)*Pi)/ymax)))/
     -  (-2*eps**2*Pi**2*xmax**2 - 2*eps**2*Pi**2*ymax**2 - 
     -    2*xmax**2*ymax**2 + 
     -    2*eps**2*Pi**2*xmax**2*Cos((4*car(1)*Pi)/xmax) - 
     -    2*eps**2*Pi**2*ymax**2*Cos((4*car(1)*Pi)/xmax) + 
     -    eps**2*Pi**2*xmax**2*
     -     Cos((4*car(1)*Pi)/xmax - (4*car(2)*Pi)/ymax) - 
     -    2*eps**2*Pi**2*xmax*ymax*
     -     Cos((4*car(1)*Pi)/xmax - (4*car(2)*Pi)/ymax) + 
     -    eps**2*Pi**2*ymax**2*
     -     Cos((4*car(1)*Pi)/xmax - (4*car(2)*Pi)/ymax) + 
     -    eps**2*Pi**2*xmax**2*
     -     Cos((4*car(1)*Pi)/xmax + (4*car(2)*Pi)/ymax) + 
     -    2*eps**2*Pi**2*xmax*ymax*
     -     Cos((4*car(1)*Pi)/xmax + (4*car(2)*Pi)/ymax) + 
     -    eps**2*Pi**2*ymax**2*
     -     Cos((4*car(1)*Pi)/xmax + (4*car(2)*Pi)/ymax) - 
     -    2*eps**2*Pi**2*xmax**2*Cos((4*car(2)*Pi)/ymax) + 
     -    2*eps**2*Pi**2*ymax**2*Cos((4*car(2)*Pi)/ymax) - 
     -    4*eps*Pi*xmax**2*ymax*
     -     Sin((2*car(1)*Pi)/xmax - (2*car(2)*Pi)/ymax) + 
     -    4*eps*Pi*xmax*ymax**2*
     -     Sin((2*car(1)*Pi)/xmax - (2*car(2)*Pi)/ymax) - 
     -    4*eps*Pi*xmax**2*ymax*
     -     Sin((2*car(1)*Pi)/xmax + (2*car(2)*Pi)/ymax) - 
     -    4*eps*Pi*xmax*ymax**2*
     -     Sin((2*car(1)*Pi)/xmax + (2*car(2)*Pi)/ymax))
              vec(3) = 0d0
              vec(4) =
     .   (-4*(eps**3*Pi**4*Cos((2*car(1)*Pi)/xmax - (6*car(2)*Pi)/ymax)+ 
     -      2*eps**3*Pi**4*Cos((2*car(1)*Pi)/xmax - (2*car(2)*Pi)/ymax)+ 
     -      eps*Pi**2*xmax**2*
     -       Cos((2*car(1)*Pi)/xmax - (2*car(2)*Pi)/ymax) + 
     -      eps**3*Pi**4*Cos((6*car(1)*Pi)/xmax - (2*car(2)*Pi)/ymax) - 
     -      2*eps**3*Pi**4*Cos((2*car(1)*Pi)/xmax + (2*car(2)*Pi)/ymax)- 
     -      eps*Pi**2*xmax**2*
     -       Cos((2*car(1)*Pi)/xmax + (2*car(2)*Pi)/ymax) - 
     -      eps**3*Pi**4*Cos((6*car(1)*Pi)/xmax + (2*car(2)*Pi)/ymax) - 
     -      eps**3*Pi**4*Cos((2*car(1)*Pi)/xmax + (6*car(2)*Pi)/ymax) + 
     -      4*eps**2*Pi**3*xmax*Sin((4*car(1)*Pi)/xmax)))/
     -  (-2*eps**2*Pi**2*xmax**2 - 2*eps**2*Pi**2*ymax**2 - 
     -    2*xmax**2*ymax**2 + 
     -    2*eps**2*Pi**2*xmax**2*Cos((4*car(1)*Pi)/xmax) - 
     -    2*eps**2*Pi**2*ymax**2*Cos((4*car(1)*Pi)/xmax) + 
     -    eps**2*Pi**2*xmax**2*
     -     Cos((4*car(1)*Pi)/xmax - (4*car(2)*Pi)/ymax) - 
     -    2*eps**2*Pi**2*xmax*ymax*
     -     Cos((4*car(1)*Pi)/xmax - (4*car(2)*Pi)/ymax) + 
     -    eps**2*Pi**2*ymax**2*
     -     Cos((4*car(1)*Pi)/xmax - (4*car(2)*Pi)/ymax) + 
     -    eps**2*Pi**2*xmax**2*
     -     Cos((4*car(1)*Pi)/xmax + (4*car(2)*Pi)/ymax) + 
     -    2*eps**2*Pi**2*xmax*ymax*
     -     Cos((4*car(1)*Pi)/xmax + (4*car(2)*Pi)/ymax) + 
     -    eps**2*Pi**2*ymax**2*
     -     Cos((4*car(1)*Pi)/xmax + (4*car(2)*Pi)/ymax) - 
     -    2*eps**2*Pi**2*xmax**2*Cos((4*car(2)*Pi)/ymax) + 
     -    2*eps**2*Pi**2*ymax**2*Cos((4*car(2)*Pi)/ymax) - 
     -    4*eps*Pi*xmax**2*ymax*
     -     Sin((2*car(1)*Pi)/xmax - (2*car(2)*Pi)/ymax) + 
     -    4*eps*Pi*xmax*ymax**2*
     -     Sin((2*car(1)*Pi)/xmax - (2*car(2)*Pi)/ymax) - 
     -    4*eps*Pi*xmax**2*ymax*
     -     Sin((2*car(1)*Pi)/xmax + (2*car(2)*Pi)/ymax) - 
     -    4*eps*Pi*xmax*ymax**2*
     -     Sin((2*car(1)*Pi)/xmax + (2*car(2)*Pi)/ymax))
              vec(5) =
     .   (-4*(eps**3*Pi**4*Cos((2*car(1)*Pi)/xmax - (6*car(2)*Pi)/ymax)+ 
     -      2*eps**3*Pi**4*Cos((2*car(1)*Pi)/xmax - (2*car(2)*Pi)/ymax)+ 
     -      eps*Pi**2*xmax*ymax*
     -       Cos((2*car(1)*Pi)/xmax - (2*car(2)*Pi)/ymax) + 
     -      eps**3*Pi**4*Cos((6*car(1)*Pi)/xmax - (2*car(2)*Pi)/ymax) - 
     -      2*eps**3*Pi**4*Cos((2*car(1)*Pi)/xmax + (2*car(2)*Pi)/ymax)+ 
     -      eps*Pi**2*xmax*ymax*
     -       Cos((2*car(1)*Pi)/xmax + (2*car(2)*Pi)/ymax) - 
     -      eps**3*Pi**4*Cos((6*car(1)*Pi)/xmax + (2*car(2)*Pi)/ymax) - 
     -      eps**3*Pi**4*Cos((2*car(1)*Pi)/xmax + (6*car(2)*Pi)/ymax) + 
     -      2*eps**2*Pi**3*xmax*Sin((4*car(1)*Pi)/xmax) + 
     -      2*eps**2*Pi**3*ymax*Sin((4*car(2)*Pi)/ymax)))/
     -  (-2*eps**2*Pi**2*xmax**2 - 2*eps**2*Pi**2*ymax**2 - 
     -    2*xmax**2*ymax**2 + 
     -    2*eps**2*Pi**2*xmax**2*Cos((4*car(1)*Pi)/xmax) - 
     -    2*eps**2*Pi**2*ymax**2*Cos((4*car(1)*Pi)/xmax) + 
     -    eps**2*Pi**2*xmax**2*
     -     Cos((4*car(1)*Pi)/xmax - (4*car(2)*Pi)/ymax) - 
     -    2*eps**2*Pi**2*xmax*ymax*
     -     Cos((4*car(1)*Pi)/xmax - (4*car(2)*Pi)/ymax) + 
     -    eps**2*Pi**2*ymax**2*
     -     Cos((4*car(1)*Pi)/xmax - (4*car(2)*Pi)/ymax) + 
     -    eps**2*Pi**2*xmax**2*
     -     Cos((4*car(1)*Pi)/xmax + (4*car(2)*Pi)/ymax) + 
     -    2*eps**2*Pi**2*xmax*ymax*
     -     Cos((4*car(1)*Pi)/xmax + (4*car(2)*Pi)/ymax) + 
     -    eps**2*Pi**2*ymax**2*
     -     Cos((4*car(1)*Pi)/xmax + (4*car(2)*Pi)/ymax) - 
     -    2*eps**2*Pi**2*xmax**2*Cos((4*car(2)*Pi)/ymax) + 
     -    2*eps**2*Pi**2*ymax**2*Cos((4*car(2)*Pi)/ymax) - 
     -    4*eps*Pi*xmax**2*ymax*
     -     Sin((2*car(1)*Pi)/xmax - (2*car(2)*Pi)/ymax) + 
     -    4*eps*Pi*xmax*ymax**2*
     -     Sin((2*car(1)*Pi)/xmax - (2*car(2)*Pi)/ymax) - 
     -    4*eps*Pi*xmax**2*ymax*
     -     Sin((2*car(1)*Pi)/xmax + (2*car(2)*Pi)/ymax) - 
     -    4*eps*Pi*xmax*ymax**2*
     -     Sin((2*car(1)*Pi)/xmax + (2*car(2)*Pi)/ymax))
              vec(6) = 0d0
              vec(7) = 0d0
              vec(8) = 0d0
              vec(9) = 0d0
            case (3)
              vec(1) = 0d0
              vec(2) = 0d0
              vec(3) =
     .   (-4*(eps*Pi**2*xmax*ymax*
     -       Cos((2*car(1)*Pi)/xmax - (2*car(2)*Pi)/ymax) - 
     -      eps*Pi**2*ymax**2*
     -       Cos((2*car(1)*Pi)/xmax - (2*car(2)*Pi)/ymax) + 
     -      eps*Pi**2*xmax*ymax*
     -       Cos((2*car(1)*Pi)/xmax + (2*car(2)*Pi)/ymax) + 
     -      eps*Pi**2*ymax**2*
     -       Cos((2*car(1)*Pi)/xmax + (2*car(2)*Pi)/ymax) + 
     -      2*eps**2*Pi**3*xmax*Sin((4*car(1)*Pi)/xmax) - 
     -      2*eps**2*Pi**3*ymax*Sin((4*car(2)*Pi)/ymax)))/
     -  (-2*eps**2*Pi**2*xmax**2 - 2*eps**2*Pi**2*ymax**2 - 
     -    2*xmax**2*ymax**2 + 
     -    2*eps**2*Pi**2*xmax**2*Cos((4*car(1)*Pi)/xmax) - 
     -    2*eps**2*Pi**2*ymax**2*Cos((4*car(1)*Pi)/xmax) + 
     -    eps**2*Pi**2*xmax**2*
     -     Cos((4*car(1)*Pi)/xmax - (4*car(2)*Pi)/ymax) - 
     -    2*eps**2*Pi**2*xmax*ymax*
     -     Cos((4*car(1)*Pi)/xmax - (4*car(2)*Pi)/ymax) + 
     -    eps**2*Pi**2*ymax**2*
     -     Cos((4*car(1)*Pi)/xmax - (4*car(2)*Pi)/ymax) + 
     -    eps**2*Pi**2*xmax**2*
     -     Cos((4*car(1)*Pi)/xmax + (4*car(2)*Pi)/ymax) + 
     -    2*eps**2*Pi**2*xmax*ymax*
     -     Cos((4*car(1)*Pi)/xmax + (4*car(2)*Pi)/ymax) + 
     -    eps**2*Pi**2*ymax**2*
     -     Cos((4*car(1)*Pi)/xmax + (4*car(2)*Pi)/ymax) - 
     -    2*eps**2*Pi**2*xmax**2*Cos((4*car(2)*Pi)/ymax) + 
     -    2*eps**2*Pi**2*ymax**2*Cos((4*car(2)*Pi)/ymax) - 
     -    4*eps*Pi*xmax**2*ymax*
     -     Sin((2*car(1)*Pi)/xmax - (2*car(2)*Pi)/ymax) + 
     -    4*eps*Pi*xmax*ymax**2*
     -     Sin((2*car(1)*Pi)/xmax - (2*car(2)*Pi)/ymax) - 
     -    4*eps*Pi*xmax**2*ymax*
     -     Sin((2*car(1)*Pi)/xmax + (2*car(2)*Pi)/ymax) - 
     -    4*eps*Pi*xmax*ymax**2*
     -     Sin((2*car(1)*Pi)/xmax + (2*car(2)*Pi)/ymax))
              vec(4) = 0d0
              vec(5) = 0d0
              vec(6) =
     .   (4*(eps*Pi**2*xmax**2*
     -       Cos((2*car(1)*Pi)/xmax - (2*car(2)*Pi)/ymax) - 
     -      eps*Pi**2*xmax*ymax*
     -       Cos((2*car(1)*Pi)/xmax - (2*car(2)*Pi)/ymax) - 
     -      eps*Pi**2*xmax**2*
     -       Cos((2*car(1)*Pi)/xmax + (2*car(2)*Pi)/ymax) - 
     -      eps*Pi**2*xmax*ymax*
     -       Cos((2*car(1)*Pi)/xmax + (2*car(2)*Pi)/ymax) + 
     -      2*eps**2*Pi**3*xmax*Sin((4*car(1)*Pi)/xmax) - 
     -      2*eps**2*Pi**3*ymax*Sin((4*car(2)*Pi)/ymax)))/
     -  (-2*eps**2*Pi**2*xmax**2 - 2*eps**2*Pi**2*ymax**2 - 
     -    2*xmax**2*ymax**2 + 
     -    2*eps**2*Pi**2*xmax**2*Cos((4*car(1)*Pi)/xmax) - 
     -    2*eps**2*Pi**2*ymax**2*Cos((4*car(1)*Pi)/xmax) + 
     -    eps**2*Pi**2*xmax**2*
     -     Cos((4*car(1)*Pi)/xmax - (4*car(2)*Pi)/ymax) - 
     -    2*eps**2*Pi**2*xmax*ymax*
     -     Cos((4*car(1)*Pi)/xmax - (4*car(2)*Pi)/ymax) + 
     -    eps**2*Pi**2*ymax**2*
     -     Cos((4*car(1)*Pi)/xmax - (4*car(2)*Pi)/ymax) + 
     -    eps**2*Pi**2*xmax**2*
     -     Cos((4*car(1)*Pi)/xmax + (4*car(2)*Pi)/ymax) + 
     -    2*eps**2*Pi**2*xmax*ymax*
     -     Cos((4*car(1)*Pi)/xmax + (4*car(2)*Pi)/ymax) + 
     -    eps**2*Pi**2*ymax**2*
     -     Cos((4*car(1)*Pi)/xmax + (4*car(2)*Pi)/ymax) - 
     -    2*eps**2*Pi**2*xmax**2*Cos((4*car(2)*Pi)/ymax) + 
     -    2*eps**2*Pi**2*ymax**2*Cos((4*car(2)*Pi)/ymax) - 
     -    4*eps*Pi*xmax**2*ymax*
     -     Sin((2*car(1)*Pi)/xmax - (2*car(2)*Pi)/ymax) + 
     -    4*eps*Pi*xmax*ymax**2*
     -     Sin((2*car(1)*Pi)/xmax - (2*car(2)*Pi)/ymax) - 
     -    4*eps*Pi*xmax**2*ymax*
     -     Sin((2*car(1)*Pi)/xmax + (2*car(2)*Pi)/ymax) - 
     -    4*eps*Pi*xmax*ymax**2*
     -     Sin((2*car(1)*Pi)/xmax + (2*car(2)*Pi)/ymax))
              vec(7) = 0d0
              vec(8) = 0d0
              vec(9) = 0d0
          end select
        case default
          write (*,*) 'Grid not implemented'
          write (*,*) 'Aborting...'
          stop
        end select

        tensor = transpose(reshape(vec, (/3,3/)))

      end function hessian_cnv

c     G_sub
c     #################################################################
      function G_sub(x1,x2,x3,cartesian) result (tensor)

c     -----------------------------------------------------------------
c     Calculates contravariant metric tensor of curvilinear coordinate 
c     system
c     -----------------------------------------------------------------

        implicit none

c     Input variables

        real(8)    :: x1,x2,x3,tensor(3,3)
        logical    :: cartesian

c     Local variables

        real(8)    :: vec(9),jac,car(3),curv(3)

c     Begin program

        select case (coords)
        case ('car')
          vec = (/ 1d0, 0d0, 0d0
     .            ,0d0, 1d0, 0d0
     .            ,0d0, 0d0, 1d0 /)
        case ('scl')
          if (cartesian) then
            curv = xi_x(x1,x2,x3)
          else
            curv = (/ x1,x2,x3 /)
          endif
          lambda = gparams(1)
          cc = 0.5/lambda
          cc = 1./tanh(cc)
          ypp = (2*curv(2)/ymax-1.)
          jac = cc*lambda/(cc**2-ypp**2)
          vec = (/ 1d0/jac, 0d0, 0d0
     .            ,0d0    , jac, 0d0
     .            ,0d0    , 0d0, 1d0/jac /)
        case ('cyl')
          if (cartesian) then
            curv = xi_x(x1,x2,x3)
          else
            curv = (/ x1,x2,x3 /)
          endif
          vec = (/ 1d0/curv(1), 0d0     , 0d0
     .            ,0d0        , curv(1) , 0d0
     .            ,0d0        , 0d0     , 1d0/curv(1) /)
        case ('hel')
          if (cartesian) then
            curv = xi_x(x1,x2,x3)
          else
            curv = (/ x1,x2,x3 /)
          endif
          mm = gparams(1)
          kk = gparams(2)
          aa = kk/mm
          vec = (/ 1./curv(1),0d0, 0d0
     .            ,0d0, curv(1) ,-aa*curv(1)
     .            ,0d0,-aa*curv(1), 1./curv(1) + aa**2*curv(1) /)
        case ('tor')
          if (cartesian) then
            curv = xi_x(x1,x2,x3)
          else
            curv = (/ x1,x2,x3 /)
          endif
          major_r = gparams(1)
          vec = (/ 1d0/curv(1)/(major_r + curv(1)*sin(curv(2))),0d0,0d0
     .            ,0d0, curv(1)/(major_r + curv(1)*sin(curv(2))), 0d0
     .            ,0d0, 0d0, (major_r/curv(1) + sin(curv(2))) /)
        case ('sin')

          pi = acos(-1d0)

          if (cartesian) then
            car = (/ x1,x2,x3 /)
          else
            car = x_xi(x1,x2,x3)
          endif
          eps = gparams(1)
          vec(1) =
     .  (2*eps**2*Pi**2*xmax**2 + 2*eps**2*Pi**2*ymax**2 + 
     -    2*xmax**2*ymax**2 - 
     -    2*eps**2*Pi**2*xmax**2*Cos((4*car(1)*Pi)/xmax) + 
     -    2*eps**2*Pi**2*ymax**2*Cos((4*car(1)*Pi)/xmax) - 
     -    eps**2*Pi**2*xmax**2*
     -     Cos((4*car(1)*Pi)/xmax - (4*car(2)*Pi)/ymax) - 
     -    eps**2*Pi**2*ymax**2*
     -     Cos((4*car(1)*Pi)/xmax - (4*car(2)*Pi)/ymax) - 
     -    eps**2*Pi**2*xmax**2*
     -     Cos((4*car(1)*Pi)/xmax + (4*car(2)*Pi)/ymax) - 
     -    eps**2*Pi**2*ymax**2*
     -     Cos((4*car(1)*Pi)/xmax + (4*car(2)*Pi)/ymax) + 
     -    2*eps**2*Pi**2*xmax**2*Cos((4*car(2)*Pi)/ymax) - 
     -    2*eps**2*Pi**2*ymax**2*Cos((4*car(2)*Pi)/ymax) + 
     -    4*eps*Pi*xmax**2*ymax*
     -     Sin((2*car(1)*Pi)/xmax - (2*car(2)*Pi)/ymax) + 
     -    4*eps*Pi*xmax**2*ymax*
     -     Sin((2*car(1)*Pi)/xmax + (2*car(2)*Pi)/ymax))/
     -  (2.*xmax*ymax*(xmax*ymax + 
     -      eps*Pi*xmax*Sin((2*car(1)*Pi)/xmax - (2*car(2)*Pi)/ymax) - 
     -      eps*Pi*ymax*Sin((2*car(1)*Pi)/xmax - (2*car(2)*Pi)/ymax) + 
     -      eps*Pi*xmax*Sin((2*car(1)*Pi)/xmax + (2*car(2)*Pi)/ymax) + 
     -      eps*Pi*ymax*Sin((2*car(1)*Pi)/xmax + (2*car(2)*Pi)/ymax)))

          vec(2) =
     .  (-2*eps**2*Pi**2*xmax**2 - 2*eps**2*Pi**2*ymax**2 + 
     -    2*eps**2*Pi**2*xmax**2*Cos((4*car(1)*Pi)/xmax) - 
     -    2*eps**2*Pi**2*ymax**2*Cos((4*car(1)*Pi)/xmax) + 
     -    eps**2*Pi**2*xmax**2*
     -     Cos((4*car(1)*Pi)/xmax - (4*car(2)*Pi)/ymax) + 
     -    eps**2*Pi**2*ymax**2*
     -     Cos((4*car(1)*Pi)/xmax - (4*car(2)*Pi)/ymax) + 
     -    eps**2*Pi**2*xmax**2*
     -     Cos((4*car(1)*Pi)/xmax + (4*car(2)*Pi)/ymax) + 
     -    eps**2*Pi**2*ymax**2*
     -     Cos((4*car(1)*Pi)/xmax + (4*car(2)*Pi)/ymax) - 
     -    2*eps**2*Pi**2*xmax**2*Cos((4*car(2)*Pi)/ymax) + 
     -    2*eps**2*Pi**2*ymax**2*Cos((4*car(2)*Pi)/ymax) - 
     -    2*eps*Pi*xmax**2*ymax*
     -     Sin((2*car(1)*Pi)/xmax - (2*car(2)*Pi)/ymax) + 
     -    2*eps*Pi*xmax*ymax**2*
     -     Sin((2*car(1)*Pi)/xmax - (2*car(2)*Pi)/ymax) - 
     -    2*eps*Pi*xmax**2*ymax*
     -     Sin((2*car(1)*Pi)/xmax + (2*car(2)*Pi)/ymax) - 
     -    2*eps*Pi*xmax*ymax**2*
     -     Sin((2*car(1)*Pi)/xmax + (2*car(2)*Pi)/ymax))/
     -  (2.*xmax*ymax*(xmax*ymax + 
     -      eps*Pi*xmax*Sin((2*car(1)*Pi)/xmax - (2*car(2)*Pi)/ymax) - 
     -      eps*Pi*ymax*Sin((2*car(1)*Pi)/xmax - (2*car(2)*Pi)/ymax) + 
     -      eps*Pi*xmax*Sin((2*car(1)*Pi)/xmax + (2*car(2)*Pi)/ymax) + 
     -      eps*Pi*ymax*Sin((2*car(1)*Pi)/xmax + (2*car(2)*Pi)/ymax)))

          vec(3) = 0d0
          vec(4) = vec(2)

          vec(5) =
     .  (2*eps**2*Pi**2*xmax**2 + 2*eps**2*Pi**2*ymax**2 + 
     -    2*xmax**2*ymax**2 - 
     -    2*eps**2*Pi**2*xmax**2*Cos((4*car(1)*Pi)/xmax) + 
     -    2*eps**2*Pi**2*ymax**2*Cos((4*car(1)*Pi)/xmax) - 
     -    eps**2*Pi**2*xmax**2*
     -     Cos((4*car(1)*Pi)/xmax - (4*car(2)*Pi)/ymax) - 
     -    eps**2*Pi**2*ymax**2*
     -     Cos((4*car(1)*Pi)/xmax - (4*car(2)*Pi)/ymax) - 
     -    eps**2*Pi**2*xmax**2*
     -     Cos((4*car(1)*Pi)/xmax + (4*car(2)*Pi)/ymax) - 
     -    eps**2*Pi**2*ymax**2*
     -     Cos((4*car(1)*Pi)/xmax + (4*car(2)*Pi)/ymax) + 
     -    2*eps**2*Pi**2*xmax**2*Cos((4*car(2)*Pi)/ymax) - 
     -    2*eps**2*Pi**2*ymax**2*Cos((4*car(2)*Pi)/ymax) - 
     -    4*eps*Pi*xmax*ymax**2*
     -     Sin((2*car(1)*Pi)/xmax - (2*car(2)*Pi)/ymax) + 
     -    4*eps*Pi*xmax*ymax**2*
     -     Sin((2*car(1)*Pi)/xmax + (2*car(2)*Pi)/ymax))/
     -  (2.*xmax*ymax*(xmax*ymax + 
     -      eps*Pi*xmax*Sin((2*car(1)*Pi)/xmax - (2*car(2)*Pi)/ymax) - 
     -      eps*Pi*ymax*Sin((2*car(1)*Pi)/xmax - (2*car(2)*Pi)/ymax) + 
     -      eps*Pi*xmax*Sin((2*car(1)*Pi)/xmax + (2*car(2)*Pi)/ymax) + 
     -      eps*Pi*ymax*Sin((2*car(1)*Pi)/xmax + (2*car(2)*Pi)/ymax)))

          vec(6) = 0d0
          vec(7) = vec(3)
          vec(8) = vec(6)

          vec(9) = (xmax*ymax + eps*Pi*xmax*
     -     Sin((2*car(1)*Pi)/xmax - (2*car(2)*Pi)/ymax) - 
     -    eps*Pi*ymax*Sin((2*car(1)*Pi)/xmax - (2*car(2)*Pi)/ymax) + 
     -    eps*Pi*xmax*Sin((2*car(1)*Pi)/xmax + (2*car(2)*Pi)/ymax) + 
     -    eps*Pi*ymax*Sin((2*car(1)*Pi)/xmax + (2*car(2)*Pi)/ymax))/
     -     (xmax*ymax)
        case default
          write (*,*) 'Grid not implemented'
          write (*,*) 'Aborting...'
          stop
        end select

        tensor = reshape(vec, (/3,3/))

      end function G_sub

c     G_sub_elem
c     #################################################################
      function G_sub_elem(i,j,x1,x2,x3,cartesian) result (element)

c     -----------------------------------------------------------------
c     Calculates contravariant metric tensor of curvilinear coordinate 
c     system
c     -----------------------------------------------------------------

        implicit none

c     Input variables

        integer(4) :: i,j
        real(8)    :: x1,x2,x3,element
        logical    :: cartesian

c     Local variables

        real(8)    :: tensor(3,3)

c     Begin program

        tensor  = G_sub(x1,x2,x3,cartesian)
        element = tensor(i,j)

      end function G_sub_elem

c     G_super
c     #################################################################
      function G_super(x1,x2,x3,cartesian) result (tensor)

c     -----------------------------------------------------------------
c     Calculates covariant metric tensor of curvilinear coordinate 
c     system
c     -----------------------------------------------------------------

        implicit none

c     Input variables

        real(8)    :: x1,x2,x3,tensor(3,3)
        logical    :: cartesian

c     Local variables

        real(8)    :: vec(9),jac,car(3),curv(3)

c     Begin program

        select case (coords)
        case ('car')
          vec = (/ 1d0, 0d0, 0d0
     .            ,0d0, 1d0, 0d0
     .            ,0d0, 0d0, 1d0 /)
       case ('scl')
          if (cartesian) then
            curv = xi_x(x1,x2,x3)
          else
            curv = (/ x1,x2,x3 /)
          endif
          lambda = gparams(1)
          cc = 0.5/lambda
          cc = 1./tanh(cc)
          ypp = (2*curv(2)/ymax-1.)
          jac = cc*lambda/(cc**2-ypp**2)
          vec = (/ jac, 0d0    , 0d0
     .            ,0d0, 1d0/jac, 0d0
     .            ,0d0, 0d0    , jac /)
        case ('cyl')
          if (cartesian) then
            curv = xi_x(x1,x2,x3)
          else
            curv = (/ x1,x2,x3 /)
          endif
          vec = (/ curv(1) , 0d0   , 0d0
     .            ,0d0, 1d0/curv(1), 0d0
     .            ,0d0, 0d0   , curv(1) /)
        case ('hel')
          if (cartesian) then
            curv = xi_x(x1,x2,x3)
          else
            curv = (/ x1,x2,x3 /)
          endif
          mm = gparams(1)
          kk = gparams(2)
          aa = kk/mm
          vec = (/ curv(1),0d0,0d0
     .            ,0d0,1./curv(1) + aa**2*curv(1), aa*curv(1)
     .            ,0d0, aa*curv(1), curv(1) /)
        case ('tor')
          if (cartesian) then
            curv = xi_x(x1,x2,x3)
          else
            curv = (/ x1,x2,x3 /)
          endif
          major_r = gparams(1)
          vec = (/ curv(1)*(major_r + curv(1)*sin(curv(2))), 0d0, 0d0
     .            ,0d0, (major_r/curv(1) + sin(curv(2))), 0d0
     .            ,0d0, 0d0, curv(1)/(major_r + curv(1)*sin(curv(2))) /)
        case ('sin')

          pi = acos(-1d0)

          if (cartesian) then
            car = (/ x1,x2,x3 /)
          else
            car = x_xi(x1,x2,x3)
          endif
          eps = gparams(1)
          vec(1) =
     .  (2*eps**2*Pi**2*xmax**2 + 2*eps**2*Pi**2*ymax**2 + 
     -    2*xmax**2*ymax**2 - 
     -    2*eps**2*Pi**2*xmax**2*Cos((4*car(1)*Pi)/xmax) + 
     -    2*eps**2*Pi**2*ymax**2*Cos((4*car(1)*Pi)/xmax) - 
     -    eps**2*Pi**2*xmax**2*
     -     Cos((4*car(1)*Pi)/xmax - (4*car(2)*Pi)/ymax) - 
     -    eps**2*Pi**2*ymax**2*
     -     Cos((4*car(1)*Pi)/xmax - (4*car(2)*Pi)/ymax) - 
     -    eps**2*Pi**2*xmax**2*
     -     Cos((4*car(1)*Pi)/xmax + (4*car(2)*Pi)/ymax) - 
     -    eps**2*Pi**2*ymax**2*
     -     Cos((4*car(1)*Pi)/xmax + (4*car(2)*Pi)/ymax) + 
     -    2*eps**2*Pi**2*xmax**2*Cos((4*car(2)*Pi)/ymax) - 
     -    2*eps**2*Pi**2*ymax**2*Cos((4*car(2)*Pi)/ymax) - 
     -    4*eps*Pi*xmax*ymax**2*
     -     Sin((2*car(1)*Pi)/xmax - (2*car(2)*Pi)/ymax) + 
     -    4*eps*Pi*xmax*ymax**2*
     -     Sin((2*car(1)*Pi)/xmax + (2*car(2)*Pi)/ymax))/
     -  (2.*xmax*ymax*(xmax*ymax + 
     -      eps*Pi*xmax*Sin((2*car(1)*Pi)/xmax - (2*car(2)*Pi)/ymax) - 
     -      eps*Pi*ymax*Sin((2*car(1)*Pi)/xmax - (2*car(2)*Pi)/ymax) + 
     -      eps*Pi*xmax*Sin((2*car(1)*Pi)/xmax + (2*car(2)*Pi)/ymax) + 
     -      eps*Pi*ymax*Sin((2*car(1)*Pi)/xmax + (2*car(2)*Pi)/ymax)))

          vec(2) =
     .  (2*eps**2*Pi**2*xmax**2 + 2*eps**2*Pi**2*ymax**2 - 
     -    2*eps**2*Pi**2*xmax**2*Cos((4*car(1)*Pi)/xmax) + 
     -    2*eps**2*Pi**2*ymax**2*Cos((4*car(1)*Pi)/xmax) - 
     -    eps**2*Pi**2*xmax**2*
     -     Cos((4*car(1)*Pi)/xmax - (4*car(2)*Pi)/ymax) - 
     -    eps**2*Pi**2*ymax**2*
     -     Cos((4*car(1)*Pi)/xmax - (4*car(2)*Pi)/ymax) - 
     -    eps**2*Pi**2*xmax**2*
     -     Cos((4*car(1)*Pi)/xmax + (4*car(2)*Pi)/ymax) - 
     -    eps**2*Pi**2*ymax**2*
     -     Cos((4*car(1)*Pi)/xmax + (4*car(2)*Pi)/ymax) + 
     -    2*eps**2*Pi**2*xmax**2*Cos((4*car(2)*Pi)/ymax) - 
     -    2*eps**2*Pi**2*ymax**2*Cos((4*car(2)*Pi)/ymax) + 
     -    2*eps*Pi*xmax**2*ymax*
     -     Sin((2*car(1)*Pi)/xmax - (2*car(2)*Pi)/ymax) - 
     -    2*eps*Pi*xmax*ymax**2*
     -     Sin((2*car(1)*Pi)/xmax - (2*car(2)*Pi)/ymax) + 
     -    2*eps*Pi*xmax**2*ymax*
     -     Sin((2*car(1)*Pi)/xmax + (2*car(2)*Pi)/ymax) + 
     -    2*eps*Pi*xmax*ymax**2*
     -     Sin((2*car(1)*Pi)/xmax + (2*car(2)*Pi)/ymax))/
     -  (2.*xmax*ymax*(xmax*ymax + 
     -      eps*Pi*xmax*Sin((2*car(1)*Pi)/xmax - (2*car(2)*Pi)/ymax) - 
     -      eps*Pi*ymax*Sin((2*car(1)*Pi)/xmax - (2*car(2)*Pi)/ymax) + 
     -      eps*Pi*xmax*Sin((2*car(1)*Pi)/xmax + (2*car(2)*Pi)/ymax) + 
     -      eps*Pi*ymax*Sin((2*car(1)*Pi)/xmax + (2*car(2)*Pi)/ymax)))

          vec(3) = 0d0
          vec(4) = vec(2)

          vec(5) =
     .  (2*eps**2*Pi**2*xmax**2 + 2*eps**2*Pi**2*ymax**2 + 
     -    2*xmax**2*ymax**2 - 
     -    2*eps**2*Pi**2*xmax**2*Cos((4*car(1)*Pi)/xmax) + 
     -    2*eps**2*Pi**2*ymax**2*Cos((4*car(1)*Pi)/xmax) - 
     -    eps**2*Pi**2*xmax**2*
     -     Cos((4*car(1)*Pi)/xmax - (4*car(2)*Pi)/ymax) - 
     -    eps**2*Pi**2*ymax**2*
     -     Cos((4*car(1)*Pi)/xmax - (4*car(2)*Pi)/ymax) - 
     -    eps**2*Pi**2*xmax**2*
     -     Cos((4*car(1)*Pi)/xmax + (4*car(2)*Pi)/ymax) - 
     -    eps**2*Pi**2*ymax**2*
     -     Cos((4*car(1)*Pi)/xmax + (4*car(2)*Pi)/ymax) + 
     -    2*eps**2*Pi**2*xmax**2*Cos((4*car(2)*Pi)/ymax) - 
     -    2*eps**2*Pi**2*ymax**2*Cos((4*car(2)*Pi)/ymax) + 
     -    4*eps*Pi*xmax**2*ymax*
     -     Sin((2*car(1)*Pi)/xmax - (2*car(2)*Pi)/ymax) + 
     -    4*eps*Pi*xmax**2*ymax*
     -     Sin((2*car(1)*Pi)/xmax + (2*car(2)*Pi)/ymax))/
     -  (2.*xmax*ymax*(xmax*ymax + 
     -      eps*Pi*xmax*Sin((2*car(1)*Pi)/xmax - (2*car(2)*Pi)/ymax) - 
     -      eps*Pi*ymax*Sin((2*car(1)*Pi)/xmax - (2*car(2)*Pi)/ymax) + 
     -      eps*Pi*xmax*Sin((2*car(1)*Pi)/xmax + (2*car(2)*Pi)/ymax) + 
     -      eps*Pi*ymax*Sin((2*car(1)*Pi)/xmax + (2*car(2)*Pi)/ymax)))

          vec(6) = 0d0
          vec(7) = vec(3)
          vec(8) = vec(6)

          vec(9) = (xmax*ymax)/
     -  (xmax*ymax + eps*Pi*xmax*
     -     Sin((2*car(1)*Pi)/xmax - (2*car(2)*Pi)/ymax) - 
     -    eps*Pi*ymax*Sin((2*car(1)*Pi)/xmax - (2*car(2)*Pi)/ymax) + 
     -    eps*Pi*xmax*Sin((2*car(1)*Pi)/xmax + (2*car(2)*Pi)/ymax) + 
     -    eps*Pi*ymax*Sin((2*car(1)*Pi)/xmax + (2*car(2)*Pi)/ymax))
        case default
          write (*,*) 'Grid not implemented'
          write (*,*) 'Aborting...'
          stop
        end select

        tensor = reshape(vec, (/3,3/))

      end function G_super

c     G_super_elem
c     #################################################################
      function G_super_elem(i,j,x1,x2,x3,cartesian) result (element)

c     -----------------------------------------------------------------
c     Calculates covariant metric tensor of curvilinear coordinate 
c     system
c     -----------------------------------------------------------------

        implicit none

c     Input variables

        integer(4) :: i,j
        real(8)    :: x1,x2,x3,element
        logical    :: cartesian

c     Local variables

        real(8)    :: tensor(3,3)

c     Begin program

        tensor  = G_super(x1,x2,x3,cartesian)
        element = tensor(i,j)

      end function G_super_elem

      end module grid_core

c module grid_interface
c #####################################################################
      module grid_interface

        use grid_core

        implicit none

        type :: grid_def
          integer(4) :: ngrdx                             !# meshes in X
          integer(4) :: ngrdy                             !# meshes in Y
          integer(4) :: ngrdz                             !# meshes in Y
          integer(4) :: ngrid                             !# meshes for MG
          integer(4),pointer,dimension(:)  :: iline       !Restrict ops. to i=iline in MG
          integer(4),pointer,dimension(:)  :: jline       !Restrict ops. to j=jline in MG
          integer(4),pointer,dimension(:)  :: kline       !Restrict ops. to k=kline in MG
          real(8)   ,pointer,dimension(:)  :: xx          !Grid node positions in X (all grids)
          real(8)   ,pointer,dimension(:)  :: yy          !Grid node positions in Y (")
          real(8)   ,pointer,dimension(:)  :: zz          !Grid node positions in Z (")
          real(8)   ,pointer,dimension(:)  :: dx          !Grid spacings in X for integer mesh (")
          real(8)   ,pointer,dimension(:)  :: dy          !Grid spacings in Y for integer mesh (")
          real(8)   ,pointer,dimension(:)  :: dz          !Grid spacings in Z for integer mesh (")
          real(8)   ,pointer,dimension(:)  :: dxh         !Grid spacings in X for half mesh (")
          real(8)   ,pointer,dimension(:)  :: dyh         !Grid spacings in Y for half mesh (")
          real(8)   ,pointer,dimension(:)  :: dzh         !Grid spacings in Z for half mesh (")
          integer(4),pointer,dimension(:)  :: nxv         !# of grid nodes in X  (")
          integer(4),pointer,dimension(:)  :: nyv         !# of grid nodes in Y  (")
          integer(4),pointer,dimension(:)  :: nzv         !# of grid nodes in Z  (")
          integer(4),pointer,dimension(:)  :: ntotv       !Total # of grid nodes (")
          integer(4),pointer,dimension(:)  :: istartx     !Pointer for MG vectors in X
          integer(4),pointer,dimension(:)  :: istarty     !Pointer for MG vectors in Y
          integer(4),pointer,dimension(:)  :: istartz     !Pointer for MG vectors in Z
          integer(4),pointer,dimension(:)  :: istartp     !Pointer for global MG vectors
          integer(4),pointer,dimension(:)  :: mg_ratio_x  !MG coarsening ratio in X
          integer(4),pointer,dimension(:)  :: mg_ratio_y  !MG coarsening ratio in Y
          integer(4),pointer,dimension(:)  :: mg_ratio_z  !MG coarsening ratio in Z
          real(8)                          :: params(5)   !Grid configuration parameters
        end type grid_def

        type (grid_def) :: grid_params

        type :: grid_metrics
          real(8),pointer,dimension(:,:,:)       :: jac   !Jacobian factor at grid cells
          real(8),pointer,dimension(:,:,:)       :: vol   !Cell volumes
          real(8),pointer,dimension(:,:,:,:,:)   :: gsub  !Covariant metric tensor at cells
          real(8),pointer,dimension(:,:,:,:,:)   :: gsup  !Contravariant metric tensor at cells
          real(8),pointer,dimension(:,:,:,:,:,:) :: Gamma !Christoffel symbol at cell centers
        end type grid_metrics

        type :: MG_grid_metrics
          type(grid_metrics),pointer,dimension(:) :: grid
        end type MG_grid_metrics

        type(MG_grid_metrics) :: gmetric

      contains

c     getMGmap
c     #################################################################
      subroutine getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)

c     -----------------------------------------------------------------
c     Inverts curvilinear coordinates to give Cartesian coordinates
c     -----------------------------------------------------------------

        implicit none

c     Input variables

        integer(4) :: i,j,k,igx,igy,igz,ig,jg,kg

c     Local variables

c     Begin program

        ig = i + grid_params%istartx(igx)
        jg = j + grid_params%istarty(igy)
        kg = k + grid_params%istartz(igz)

      end subroutine getMGmap

c     getCoordinates
c     #################################################################
      subroutine getCoordinates(i,j,k,igx,igy,igz,ig,jg,kg,x1,y1,z1
     .                         ,cartesian)

c     -----------------------------------------------------------------
c     Finds optimal coordinates (cartesian,curvilinear) for calculation
c     of grid quantities, depending on grid definition. If analytical
c     inverse map is available, gives curvilinear (cartesian=.false.);
c     otherwise, gives cartesian.
c     -----------------------------------------------------------------

        implicit none

c     Input variables

        integer(4) :: i,j,k,igx,igy,igz,ig,jg,kg
        real(8)    :: x1,y1,z1
        logical    :: cartesian

c     Local variables

        real(8)    :: car(3)

c     Begin program

        select case (coords)
        case ('sin')
          call getCartesianCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
     .                                ,x1,y1,z1)
          cartesian = .true.
        case default
          call getCurvilinearCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
     .                                  ,x1,y1,z1)
          cartesian = .false.
        end select

      end subroutine getCoordinates

c     getCurvilinearCoordinates
c     #################################################################
      subroutine getCurvilinearCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
     .                                    ,x1,y1,z1)

c     -----------------------------------------------------------------
c     Finds curvilinear coordinates for position (i,j,k)
c     -----------------------------------------------------------------

        implicit none

c     Input variables

        integer(4) :: i,j,k,igx,igy,igz,ig,jg,kg
        real(8)    :: x1,y1,z1

c     Local variables

        integer(4) :: ii,jj,ny

c     Begin program

        call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)

        x1 = grid_params%xx(ig)
        y1 = grid_params%yy(jg)
        z1 = grid_params%zz(kg)

      end subroutine getCurvilinearCoordinates

c     getCartesianCoordinates
c     #################################################################
      subroutine getCartesianCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
     .                                  ,x1,y1,z1)

c     -----------------------------------------------------------------
c     Inverts curvilinear coordinates to give Cartesian coordinates
c     -----------------------------------------------------------------

        implicit none

c     Input variables

        integer(4) :: i,j,k,igx,igy,igz,ig,jg,kg
        real(8)    :: x1,y1,z1

c     Local variables

        real(8)    :: car(3)

c     Begin program

        call getCurvilinearCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
     .                                ,x1,y1,z1)

        car = x_xi(x1,y1,z1)

        x1 = car(1)
        y1 = car(2)
        z1 = car(3)

      end subroutine getCartesianCoordinates

c     transformVectorToCartesian
c     #################################################################
      subroutine transformVectorToCartesian(i,j,k,igx,igy,igz
     .                                     ,c1,c2,c3,covariant
     .                                     ,cx,cy,cz)

c     -----------------------------------------------------------------
c     Interface for transformVecToCar that reads logical indeces (i,j,k)
c     instead of spatial coordinates.
c     -----------------------------------------------------------------

        implicit none

c     Input variables

        integer(4) :: i,j,k,igx,igy,igz
        real(8)    :: c1,c2,c3,cx,cy,cz
        logical    :: covariant
        

c     Local variables

        integer(4) :: ig,jg,kg
        real(8)    :: x,y,z
        logical    :: cartesian

c     Begin program

        call getCoordinates(i,j,k,igx,igy,igz,ig,jg,kg,x,y,z,cartesian)

        call transformVecToCar(x,y,z,cartesian,c1,c2,c3,covariant
     .                        ,cx,cy,cz)

      end subroutine transformVectorToCartesian

c     transformVecToCar
c     #################################################################
      subroutine transformVecToCar(x,y,z,cartesian
     .                            ,c1,c2,c3,covariant
     .                            ,cx,cy,cz)

c     -----------------------------------------------------------------
c     Transforms a curvilinear vector (c1,c2,c3) to Cartesian (cx,cy,cz)
c     at spatial coordinates (x,y,z) (cartesian if cartesian=true, logical
c     otherwise).
c     -----------------------------------------------------------------

        implicit none

c     Input variables

        real(8)    :: x,y,z,c1,c2,c3,cx,cy,cz
        logical    :: covariant,cartesian
        

c     Local variables

        integer(4) :: ic
        real(8)    :: T_to_car(3,3),vec(3)

c     Begin program

        if (covariant) then

          do ic =1,3
            T_to_car(:,ic) = covariantVector(ic,x,y,z,cartesian)
          enddo

          vec = (/ c1,c2,c3 /)

          vec = matmul(T_to_car,vec)

          cx = vec(1)
          cy = vec(2)
          cz = vec(3)
            
        else

          do ic =1,3
            T_to_car(:,ic) = contravariantVector(ic,x,y,z,cartesian)
          enddo

          vec = (/ c1,c2,c3 /)

          vec = matmul(T_to_car,vec)

          cx = vec(1)
          cy = vec(2)
          cz = vec(3)

        endif

      end subroutine transformVecToCar

c     transformVectorToCurvilinear
c     #################################################################
      subroutine transformVectorToCurvilinear(i,j,k,igx,igy,igz
     .                                       ,cx,cy,cz,covariant
     .                                       ,c1,c2,c3)

c     -----------------------------------------------------------------
c     Interface for transformVec2Curv that reads logical indeces (i,j,k)
c     instead of spatial coordinates.
c     -----------------------------------------------------------------

        implicit none

c     Input variables

        integer(4) :: i,j,k,igx,igy,igz
        real(8)    :: c1,c2,c3,cx,cy,cz
        logical    :: covariant
        

c     Local variables

        integer(4) :: ig,jg,kg
        real(8)    :: x,y,z
        logical    :: cartesian

c     Begin program

        call getCoordinates(i,j,k,igx,igy,igz,ig,jg,kg,x,y,z,cartesian)

        call transformVec2Curv(x,y,z,cartesian,cx,cy,cz
     .                        ,c1,c2,c3,covariant)

      end subroutine transformVectorToCurvilinear

c     transformVec2Curv
c     #################################################################
      subroutine transformVec2Curv(x,y,z,cartesian,cx,cy,cz
     .                            ,c1,c2,c3,covariant)

c     -----------------------------------------------------------------
c     Transforms a Cartesian vector (cx,cy,cz) to curvilinear (c1,c2,c3)
c     at spatial coordinates x,y,z (cartesian if cartesian=.true., logical
c     otherwise).
c     -----------------------------------------------------------------

        implicit none

c     Input variables

        real(8)    :: x,y,z,c1,c2,c3,cx,cy,cz
        logical    :: covariant,cartesian
        

c     Local variables

        integer(4) :: ic
        real(8)    :: T_to_curv(3,3),vec(3),jac

c     Begin program

        jac = jacobian(x,y,z,cartesian)

        if (covariant) then

          do ic =1,3
            T_to_curv(:,ic) = contravariantVector(ic,x,y,z,cartesian)
          enddo

          T_to_curv = jac*transpose(T_to_curv)

          vec = (/ cx,cy,cz /)

          vec = matmul( T_to_curv,vec)

          c1 = vec(1)
          c2 = vec(2)
          c3 = vec(3)
            
        else

          do ic =1,3
            T_to_curv(:,ic) = covariantVector(ic,x,y,z,cartesian)
          enddo

          T_to_curv = jac*transpose(T_to_curv)

          vec = (/ cx,cy,cz /)

          vec = matmul(T_to_curv,vec)

          c1 = vec(1)
          c2 = vec(2)
          c3 = vec(3)

        endif

      end subroutine transformVec2Curv

c     transformFromCurvToCurv
c     #################################################################
      subroutine transformFromCurvToCurv(i,j,k,igx,igy,igz
     .             ,cov1,cov2,cov3,cnv1,cnv2,cnv3,tocnv)
c     -----------------------------------------------------------------
c     Interface for transformCurvToCurv that reads logical indices
c     (i,j,k).
c     -----------------------------------------------------------------

        implicit none

c     Input variables

        integer(4) :: i,j,k,igx,igy,igz
        real(8)    :: cov1,cov2,cov3,cnv1,cnv2,cnv3
        logical    :: tocnv

c     Local variables

        integer(4) :: ig,jg,kg
        real(8)    :: x,y,z
        logical    :: cartesian

c     Begin program

        call getCoordinates(i,j,k,igx,igy,igz,ig,jg,kg,x,y,z,cartesian)

        call transformCurvToCurv(x,y,z,cartesian
     .             ,cov1,cov2,cov3,cnv1,cnv2,cnv3,tocnv)

      end subroutine transformFromCurvToCurv

c     transformCurvToCurv
c     #################################################################
      subroutine transformCurvToCurv(x,y,z,cartesian
     .             ,cov1,cov2,cov3,cnv1,cnv2,cnv3,tocnv)
c     -----------------------------------------------------------------
c     Transforms a curvilinear vector from covariant to contravariant 
c     (tocnv=.true.) or viceversa at spatial coordinates (x,y,z)
c     (cartesian if cartesian=true, logical otherwise).
c     -----------------------------------------------------------------

        implicit none

c     Input variables

        real(8)    :: x,y,z,cov1,cov2,cov3,cnv1,cnv2,cnv3
        logical    :: tocnv,cartesian

c     Local variables

        real(8)    :: gmat(3,3),cov(3),cnv(3)

c     Begin program

        if (tocnv) then
          cov = (/ cov1,cov2,cov3 /)
          gmat= G_super(x,y,z,cartesian)
          cnv = matmul(gmat,cov)
          cnv1 = cnv(1)
          cnv2 = cnv(2)
          cnv3 = cnv(3)
        else
          cnv  = (/ cnv1,cnv2,cnv3 /)
          gmat = G_sub(x,y,z,cartesian)
          cov  = matmul(gmat,cnv)
          cov1 = cov(1)
          cov2 = cov(2)
          cov3 = cov(3)
        endif

      end subroutine transformCurvToCurv

c     volume
c     #################################################################
      function volume(i,j,k,igx,igy,igz) result(vol)

c     -----------------------------------------------------------------
c     Calculates Jacobian of curvilinear coordinate system
c     -----------------------------------------------------------------

        implicit none

c     Input variables

        integer(4) :: i,j,k,igx,igy,igz
        real(8)    :: vol

c     Local variables

        integer(4) :: ig,jg,kg
        real(8)    :: x1,x2,x3,dx1,dx2,dx3,jac
        logical    :: cartesian

c     Begin program

        call getCoordinates(i,j,k,igx,igy,igz,ig,jg,kg,x1,x2,x3
     .                     ,cartesian)

        dx1 = grid_params%dxh(ig)
        dx2 = grid_params%dyh(jg)
        dx3 = grid_params%dzh(kg)

        jac = jacobian(x1,x2,x3,cartesian)
        
        vol = jac*dx1*dx2*dx3

      end function volume

c     vectorNorm
c     ################################################################
      real(8) function vectorNorm(x1,y1,z1,ax,ay,az,covar,cartesian)

c     ---------------------------------------------------------------
c     Finds norm of vector A given its curvilinear components.
c     ---------------------------------------------------------------

      implicit none

c     Call variables

      real(8)    :: x1,y1,z1,ax,ay,az
      logical    :: covar,cartesian

c     Local variables

      real(8)    :: tensor(3,3),cnv(3),cov(3)

c     Begin program

      if (covar) then
        tensor = G_super(x1,y1,z1,cartesian)
        cov = (/ ax,ay,az /)
        cnv = matmul(tensor,cov)
      else
        tensor = G_sub(x1,y1,z1,cartesian)
        cnv = (/ ax,ay,az /)
        cov = matmul(tensor,cnv)
      endif

      vectorNorm = dot_product(cov,cnv)/jacobian(x1,y1,z1,cartesian)

c     End 

      end function vectorNorm

c     scalarProduct
c     ################################################################
      function scalarProduct(x1,y1,z1,cov1,cov2,cov3,cnv1,cnv2,cnv3
     .                      ,cartesian) result (dot)

c     ---------------------------------------------------------------
c     Finds scalar product of two vectors, one covariant and the
c     other contravariant.
c     ---------------------------------------------------------------

      implicit none

c     Call variables

      real(8)    :: dot,cov1,cov2,cov3,cnv1,cnv2,cnv3,x1,y1,z1
      logical    :: cartesian

c     Local variables

      real(8)    :: cnv(3),cov(3)

c     Begin program

      cnv = (/ cnv1,cnv2,cnv3 /)
      cov = (/ cov1,cov2,cov3 /)

      dot = dot_product(cov,cnv)/jacobian(x1,y1,z1,cartesian)

c     End 

      end function scalarProduct

      end module grid_interface

c module grid
c #####################################################################
      module grid

        use grid_interface

        implicit none

        integer(4) :: PER,DIR,NEU,SP,EQU,DEF
        parameter (PER=2,SP=5,DIR=4,NEU=3,EQU=1,DEF=6)

        integer(4),private :: nxx,nyy,nzz

        integer(4) :: mg_ratio,bcond(6)

        real(8),private :: pi

        INTERFACE ASSIGNMENT (=)
          module procedure equateGridStructure
        END INTERFACE

      contains

c     createGrid
c     #################################################################
      subroutine createGrid(nx,ny,nz)

c     -----------------------------------------------------------------
c     Defines logical grid and finds grid quantities
c     -----------------------------------------------------------------

        implicit none

c     Call variables

        integer(4) :: nx,ny,nz

c     Local variables

        integer(4) :: i,n1,n2,n3
        integer(4) :: ngrdx,ngrdy,ngrdz

c     Begin program

        nxx = nx
        nyy = ny
        nzz = nz

        pi = acos(-1d0)

c     Find adequate number of grid levels (for MG)

        n1 = int(dlog(1d0*nxx)/dlog(1d0*mg_ratio)+0.001)
        n2 = int(dlog(1d0*nyy)/dlog(1d0*mg_ratio)+0.001)
        n3 = int(dlog(1d0*nzz)/dlog(1d0*mg_ratio)+0.001)

        ngrdx = max(n1-1,1)
        do i = ngrdx,1,-1
          n1 = nxx/mg_ratio**(i-1)
          if (n1*mg_ratio**(i-1).eq.nxx) exit
        enddo
        ngrdx = i

        ngrdy = max(n2-1,1)
        do i = ngrdy,1,-1
          n2 = nyy/mg_ratio**(i-1)
          if (n2*mg_ratio**(i-1).eq.nyy) exit
        enddo
        ngrdy = i

        ngrdz = max(n3-1,1)
        do i = ngrdz,1,-1
          n3 = nzz/mg_ratio**(i-1)
          if (n3*mg_ratio**(i-1).eq.nzz) exit
        enddo
        ngrdz = i

c     Allocate grid storage structure

        call allocateGridStructure(ngrdx,ngrdy,ngrdz,grid_params)

c     Initialize MG arrays

        grid_params%mg_ratio_x = 1
        grid_params%mg_ratio_y = 1
        grid_params%mg_ratio_z = 1

        grid_params%nxv(1) = nxx
        do i = 2,ngrdx
          grid_params%nxv(i) = grid_params%nxv(i-1) / mg_ratio
          grid_params%mg_ratio_x(i-1) = mg_ratio
        enddo
        do i = ngrdx+1,grid_params%ngrid
          grid_params%nxv(i) = grid_params%nxv(i-1)
        enddo

        grid_params%nyv(1) = nyy
        do i = 2,ngrdy
          grid_params%nyv(i) = grid_params%nyv(i-1) / mg_ratio
          grid_params%mg_ratio_y(i-1) = mg_ratio
        enddo
        do i = ngrdy+1,grid_params%ngrid
          grid_params%nyv(i) = grid_params%nyv(i-1)
        enddo

        grid_params%nzv(1) = nzz
        do i = 2,ngrdz
          grid_params%nzv(i) = grid_params%nzv(i-1) / mg_ratio
          grid_params%mg_ratio_z(i-1) = mg_ratio
        enddo
        do i = ngrdz+1,grid_params%ngrid
          grid_params%nzv(i) = grid_params%nzv(i-1)
        enddo

        grid_params%istartx(1) = 1
        do i = 2,grid_params%ngrid
          grid_params%istartx(i) = grid_params%istartx(i-1)
     .                            +(grid_params%nxv(i-1)+2)
        enddo

        grid_params%istarty(1) = 1
        do i = 2,grid_params%ngrid
          grid_params%istarty(i) = grid_params%istarty(i-1)
     .                            +(grid_params%nyv(i-1)+2)
        enddo

        grid_params%istartz(1) = 1
        do i = 2,grid_params%ngrid
          grid_params%istartz(i) = grid_params%istartz(i-1)
     .                            +(grid_params%nzv(i-1)+2)
        enddo

        grid_params%istartp(1) = 1
        do i = 2,grid_params%ngrid
          grid_params%istartp(i) = grid_params%istartp(i-1)
     .                            +grid_params%nxv(i-1)
     .                            *grid_params%nyv(i-1)
     .                            *grid_params%nzv(i-1)
        enddo

c     Set grid parameters

        grid_params%params = gparams

c     Consistency checks

        if (xmax == 0d0) then
          xmax = 2*pi
          xmin = 0d0
        endif
        if (ymax == 0d0) then
          ymax = 2*pi
          ymin = 0d0
        endif
        if (zmax == 0d0) then
          zmax = 2*pi
          zmin = 0d0
        endif

        call consistencyCheck

c     Define uniform logical grid on ALL grid levels

        call createLogicalGrid(nxx,grid_params%xx,grid_params%dx
     .                        ,grid_params%dxh,grid_params%nxv
cc     .                        ,grid_params%ngrdx,grid_params%istartx
     .                        ,grid_params%ngrid,grid_params%istartx
     .                        ,xmin,xmax,bcond(1),bcond(2))

        call createLogicalGrid(nyy,grid_params%yy,grid_params%dy
     .                        ,grid_params%dyh,grid_params%nyv
cc     .                        ,grid_params%ngrdy,grid_params%istarty
     .                        ,grid_params%ngrid,grid_params%istarty
     .                        ,ymin,ymax,bcond(3),bcond(4))

        call createLogicalGrid(nzz,grid_params%zz,grid_params%dz
     .                        ,grid_params%dzh,grid_params%nzv
cc     .                        ,grid_params%ngrdz,grid_params%istartz
     .                        ,grid_params%ngrid,grid_params%istartz
     .                        ,zmin,zmax,bcond(5),bcond(6))

      end subroutine createGrid

c     allocateGridStructure
c     #################################################################
      subroutine allocateGridStructure(ngridx,ngridy,ngridz,grid_st)
c     -----------------------------------------------------------------
c     Allocates grid structure
c     -----------------------------------------------------------------

        implicit none

c     Call variables

        integer(4)     :: ngridx,ngridy,ngridz
        type(grid_def) :: grid_st

c     Local variables

        integer(4) :: ngrid,nxmg,nymg,nzmg

c     Begin program

        ngrid = max(ngridx,ngridy,ngridz)

        grid_st%ngrdx = ngridx
        grid_st%ngrdy = ngridy
        grid_st%ngrdz = ngridz
        grid_st%ngrid = ngrid

        nxmg=findMGsize(nxx,ngridx,ngrid)
        nymg=findMGsize(nyy,ngridy,ngrid)
        nzmg=findMGsize(nzz,ngridz,ngrid)

        if (.not.associated(grid_st%xx)) then
          allocate(grid_st%xx(nxmg+2*ngrid))
          allocate(grid_st%yy(nymg+2*ngrid))
          allocate(grid_st%zz(nzmg+2*ngrid))
          allocate(grid_st%dx(nxmg+2*ngrid))
          allocate(grid_st%dy(nymg+2*ngrid))
          allocate(grid_st%dz(nzmg+2*ngrid))
          allocate(grid_st%dxh(nxmg+2*ngrid))
          allocate(grid_st%dyh(nymg+2*ngrid))
          allocate(grid_st%dzh(nzmg+2*ngrid))
          allocate(grid_st%nxv(ngrid))
          allocate(grid_st%nyv(ngrid))
          allocate(grid_st%nzv(ngrid))
          allocate(grid_st%ntotv(ngrid))
          allocate(grid_st%istartx(ngrid))
          allocate(grid_st%istarty(ngrid))
          allocate(grid_st%istartz(ngrid))
          allocate(grid_st%istartp(ngrid))
          allocate(grid_st%mg_ratio_x(ngrid))
          allocate(grid_st%mg_ratio_y(ngrid))
          allocate(grid_st%mg_ratio_z(ngrid))
          allocate(grid_st%iline(ngrid))
          allocate(grid_st%jline(ngrid))
          allocate(grid_st%kline(ngrid))
        endif

        grid_st%iline = 0
        grid_st%jline = 0
        grid_st%kline = 0

c     End program

      contains

c     findMGsize
c     #################################################################
      function findMGsize(nn,ngrd,ngrdt) result (nnmg)
      implicit none
c     -----------------------------------------------------------------
c     Finds size for MG vectors, taking into account total grid levels
c     ngrdt, grid levels in the relevant direction ngrd, and the 
c     number of mesh points in the finest grid nn. The formula ensures
c     enough space even if nn=1 in the finest grid. It does NOT include
c     ghost cells (this requires an additional term of 2*ngrdt).
c     -----------------------------------------------------------------

        integer(4) :: nn,nnmg,ngrd,ngrdt

        nnmg = 2*nn + ngrdt - ngrd -1

      end function findMGsize

      end subroutine allocateGridStructure

c     equateGridStructure
c     #################################################################
      subroutine equateGridStructure(grid_st2,grid_st1)

c     -----------------------------------------------------------------
c     Performs grid_st2 = grid_st1
c     -----------------------------------------------------------------

        implicit none

c     Call variables

        type(grid_def),intent(in)  :: grid_st1
        type(grid_def),intent(out) :: grid_st2

c     Local variables

        integer(4)     :: ngridx,ngridy,ngridz

c     Begin program

        ngridx = grid_st1%ngrdx
        ngridy = grid_st1%ngrdy
        ngridz = grid_st1%ngrdz

        call allocateGridStructure(ngridx,ngridy,ngridz,grid_st2)

        grid_st2%iline      = grid_st1%iline
        grid_st2%jline      = grid_st1%jline
        grid_st2%kline      = grid_st1%kline
        grid_st2%xx         = grid_st1%xx        
        grid_st2%yy         = grid_st1%yy        
        grid_st2%zz         = grid_st1%zz        
        grid_st2%dx         = grid_st1%dx        
        grid_st2%dy         = grid_st1%dy        
        grid_st2%dz         = grid_st1%dz        
        grid_st2%dxh        = grid_st1%dxh       
        grid_st2%dyh        = grid_st1%dyh       
        grid_st2%dzh        = grid_st1%dzh       
        grid_st2%nxv        = grid_st1%nxv       
        grid_st2%nyv        = grid_st1%nyv       
        grid_st2%nzv        = grid_st1%nzv       
        grid_st2%ntotv      = grid_st1%ntotv     
        grid_st2%istartx    = grid_st1%istartx   
        grid_st2%istarty    = grid_st1%istarty   
        grid_st2%istartz    = grid_st1%istartz   
        grid_st2%istartp    = grid_st1%istartp   
        grid_st2%mg_ratio_x = grid_st1%mg_ratio_x
        grid_st2%mg_ratio_y = grid_st1%mg_ratio_y
        grid_st2%mg_ratio_z = grid_st1%mg_ratio_z

c     End program

      end subroutine equateGridStructure

c     writeGridStructure
c     #################################################################
      subroutine writeGridStructure(grid_st)

c     -----------------------------------------------------------------
c     Performs grid_st2 = grid_st1
c     -----------------------------------------------------------------

        implicit none

c     Call variables

        type(grid_def) :: grid_st

c     Local variables

c     Begin program

        write (*,*) 'ngrdx',grid_st%ngrdx
        write (*,*) 'ngrdy',grid_st%ngrdy
        write (*,*) 'ngrdz',grid_st%ngrdz
        write (*,*) 'xx',grid_st%xx        
        write (*,*) 'yy',grid_st%yy        
        write (*,*) 'zz',grid_st%zz        
        write (*,*) 'dx',grid_st%dx        
        write (*,*) 'dy',grid_st%dy        
        write (*,*) 'dz',grid_st%dz        
        write (*,*) 'dxh',grid_st%dxh       
        write (*,*) 'dyh',grid_st%dyh       
        write (*,*) 'dzh',grid_st%dzh       
        write (*,*) 'nxv',grid_st%nxv       
        write (*,*) 'nyv',grid_st%nyv       
        write (*,*) 'nzv',grid_st%nzv       
        write (*,*) 'ntotv',grid_st%ntotv     
        write (*,*) 'istartx',grid_st%istartx   
        write (*,*) 'istarty',grid_st%istarty   
        write (*,*) 'istartz',grid_st%istartz   
        write (*,*) 'istartp',grid_st%istartp   
        write (*,*) 'mg_ratio_x',grid_st%mg_ratio_x
        write (*,*) 'mg_ratio_y',grid_st%mg_ratio_y
        write (*,*) 'mg_ratio_z',grid_st%mg_ratio_z

c     End program

      end subroutine writeGridStructure

c     allocateGridMetric
c     #################################################################
      subroutine allocateGridMetric(gmetric)

        implicit none

c     Call variables

        type(MG_grid_metrics) :: gmetric

c     Local variables

        integer(4)      :: igrid,nxp,nyp,nzp

c     Begin program

        if (.not.associated(gmetric%grid)) then
          allocate(gmetric%grid(grid_params%ngrid))
        endif

        do igrid=1,grid_params%ngrid
          if (.not.associated(gmetric%grid(igrid)%jac)) then
            nxp = grid_params%nxv(igrid)+1
            nyp = grid_params%nyv(igrid)+1
            nzp = grid_params%nzv(igrid)+1
            allocate(gmetric%grid(igrid)%jac  (0:nxp,0:nyp,0:nzp))
            allocate(gmetric%grid(igrid)%vol  (0:nxp,0:nyp,0:nzp))
            allocate(gmetric%grid(igrid)%gsub (0:nxp,0:nyp,0:nzp,3,3))
            allocate(gmetric%grid(igrid)%gsup (0:nxp,0:nyp,0:nzp,3,3))
            allocate(gmetric%grid(igrid)%Gamma(0:nxp,0:nyp,0:nzp,3,3,3))
          endif
        enddo

c     End program

      end subroutine allocateGridMetric

c     defineGridMetric
c     #################################################################
      subroutine defineGridMetric

        implicit none

c     Call variables

cc        type(MG_grid_metrics) :: gmetric

c     Local variables

        integer(4)      :: igrid,nxp,nyp,nzp,i,j,k,ig,jg,kg
     .                    ,igx,igy,igz
        real(8)         :: x0,y0,z0
        logical         :: cartsn

c     Begin program

        call allocateGridMetric(gmetric)

        do igrid=1,grid_params%ngrid

          igx = igrid
          igy = igrid
          igz = igrid

          do k = 0,grid_params%nzv(igrid)+1
            do j = 0,grid_params%nyv(igrid)+1
              do i = 0,grid_params%nxv(igrid)+1

                call getCoordinates(i,j,k,igx,igy,igz,ig,jg,kg,x0,y0,z0
     .                             ,cartsn)

                gmetric%grid(igrid)%jac  (i,j,k)
     .                      = jacobian(x0,y0,z0,cartsn)
                gmetric%grid(igrid)%vol  (i,j,k)
     .                      = volume(i,j,k,igx,igy,igz)
                gmetric%grid(igrid)%gsub (i,j,k,:,:)
     .                      = g_sub   (x0,y0,z0,cartsn)
                gmetric%grid(igrid)%gsup (i,j,k,:,:)
     .                      = g_super (x0,y0,z0,cartsn)
                gmetric%grid(igrid)%Gamma(i,j,k,1,:,:)
     .                      = hessian(1,x0,y0,z0,cartsn)
                gmetric%grid(igrid)%Gamma(i,j,k,2,:,:)
     .                      = hessian(2,x0,y0,z0,cartsn)
                gmetric%grid(igrid)%Gamma(i,j,k,3,:,:)
     .                      = hessian(3,x0,y0,z0,cartsn)
              enddo
            enddo
          enddo

        enddo

c     End program

      end subroutine defineGridMetric

c     createLogicalGrid
c     #################################################################
      subroutine createLogicalGrid (nn,xx,dx,dxh,nx,ngrid,istart
     .                             ,lmin,lmax,bcs1,bcs2)

        implicit none

c     Call variables

        integer(4) :: nn,ngrid,nx(ngrid),istart(ngrid),bcs1,bcs2
        real(8)    :: xx(*),dx(*),dxh(*),lmin,lmax

c     Local variables
        
        integer(4) :: ig,i,isig
        real(8)    :: dh,length

c     Begin program

        length = lmax-lmin

c     Periodic

        if (bcs1 == PER .or. bcs2 == PER) then

          do ig = 1,ngrid

            isig = istart(ig)

          !Find cell centers
            dh = length/dfloat(nx(ig))

            xx(1 + isig) = lmin
            do i = 2,nx(ig)+1
              xx(i + isig) = xx(i-1 + isig) + dh
            enddo
            xx(0 + isig) = xx(1+isig) - dh
cc            xx(0 + isig) = xx(nx(ig)+isig)
cc            xx(nx(ig)+1 +isig) = xx(1+isig)

          !Find integer mesh spacings
            do i = 1,nx(ig)
              dx(i + isig) = xx(i+1 + isig) - xx(i + isig)
            enddo
            dx(0       +isig) = dx(nx(ig)+isig)
            dx(nx(ig)+1+isig) = dx(1     +isig)

          !Find half mesh spacings
            do i = 1,nx(ig)+1
              dxh(i + isig) = (dx(i + isig) + dx(i-1 + isig))/2.
            enddo
            dxh(0 + isig) = dxh(nx(ig) + isig)

          enddo

c     Radial (singular point at r=0)

        elseif (bcs1 == SP) then

          do ig = 1,ngrid

            isig = istart(ig)

          !Find cell centers
            dh = length/dfloat(nx(ig))

            xx(0 + isig) = 1d-8*dh
            xx(1 + isig) = dh/2.
            do i = 2,nx(ig)+1
cc            xx(0 + isig) = -dh/2.
cc            do i = 1,nx(ig)+1
              xx(i + isig) = xx(i-1 + isig) + dh
            enddo

          !Find integer mesh spacings
cc            dx(0 + isig) = dh/2.
cc            do i = 1,nx(ig)
            do i = 0,nx(ig)
              dx(i + isig) = xx(i+1 + isig) - xx(i + isig)
            enddo

          !Find half mesh spacings
            dxh(1 + isig) = dh
            do i = 2,nx(ig)
cc            do i = 1,nx(ig)
              dxh(i + isig) = (dx(i + isig) + dx(i-1 + isig))/2.
            enddo
            dxh(0        + isig) = dx(0      + isig)/2.
            dxh(nx(ig)+1 + isig) = dx(nx(ig) + isig)/2.

cc            xx(0 + isig) = 1d-8
          enddo

c     Other

        else

          do ig = 1,ngrid

            isig = istart(ig)

          !Find cell centers
c-ncv            dh = length/dfloat(nx(ig)+1)
            dh = length/dfloat(nx(ig))

c-ncv            xx(0 + isig) = lmin
            xx(0 + isig) = lmin-dh/2.
            do i = 1,nx(ig)+1
              xx(i + isig) = xx(i-1 + isig) + dh
            enddo

          !Find integer mesh spacings
            do i = 0,nx(ig)
              dx(i + isig) = xx(i+1 + isig) - xx(i + isig)
            enddo

          !Find half mesh spacings
            do i = 1,nx(ig)
              dxh(i + isig) = (dx(i + isig) + dx(i-1 + isig))/2.
            enddo
            dxh(0        + isig) = dx(0      + isig)/2.
            dxh(nx(ig)+1 + isig) = dx(nx(ig) + isig)/2.
          enddo

        endif

      end subroutine createLogicalGrid

c     consistencyCheck
c     #################################################################
      subroutine consistencyCheck

c     -----------------------------------------------------------------
c     Checks consistency of grid parameters
c     -----------------------------------------------------------------

        implicit none

c     Input variables

c     Local variables

        real(8) :: major_r

c     Begin program

        select case (coords)
        case ('car')
        case ('scl')
        case ('cyl')
          if (xmin /= 0d0 .and. bcond(1) == SP) then
            write (*,*) 'Error in setup: xmin =/0 is not singular point'
            write (*,*) 'Aborting...'
            stop
          endif
        case ('hel')
          if (xmin /= 0d0 .and. bcond(1) == SP) then
            write (*,*) 'Error in setup: xmin =/0 is not singular point'
            write (*,*) 'Aborting...'
            stop
          endif
        case ('tor')
          if (xmin /= 0d0 .and. bcond(1) == SP) then
            write (*,*) 'Error in setup: xmin =/0 is not singular point'
            write (*,*) 'Aborting...'
            stop
          endif

          major_r = grid_params%params(1)

          if (major_r < xmax) then
            write (*,*) 'Ill-defined toroidal coordinate system'
            write (*,*) 'Major radius < minor radius'
            write (*,*) 'Aborting'
            stop
          endif

        case ('sin')
        case default
          write (*,*) 'Grid not implemented'
          write (*,*) 'Aborting...'
          stop
        end select

      end subroutine consistencyCheck

c     checkGrid
c     #################################################################
      subroutine checkGrid

c     -----------------------------------------------------------------
c     Defines logical grid and finds grid quantities
c     -----------------------------------------------------------------

        implicit none

c     Call variables

c     Local variables

        integer(4) :: i,j,n1,n2,n3
        integer(4) :: igx,igy,igz,isig

c     Begin program

c     Multigrid parameters

cc        igx = grid_params%ngrdx
cc        igy = grid_params%ngrdy
cc        igz = grid_params%ngrdz
        igx = grid_params%ngrid
        igy = grid_params%ngrid
        igz = grid_params%ngrid

        write (*,*)
        write (*,*) 'Coordinate system: ',coords
        write (*,*)
        write (*,*) 'Number of grid levels:'
        write (*,*) 'nx: ',igx,'   ny: ',igy,'   nz: ',igz

        call gridInfo('X',igx,grid_params%istartx,grid_params%nxv
     .               ,grid_params%xx,grid_params%dx,grid_params%dxh)

        call gridInfo('Y',igy,grid_params%istarty,grid_params%nyv
     .               ,grid_params%yy,grid_params%dy,grid_params%dyh)

        call gridInfo('Z',igz,grid_params%istartz,grid_params%nzv
     .               ,grid_params%zz,grid_params%dz,grid_params%dzh)

        call metricTensorCheck(igx,igy,igz
     .             ,grid_params%nxv,grid_params%nyv,grid_params%nzv)

        call hessianCheck(igx,igy,igz
     .             ,grid_params%nxv,grid_params%nyv,grid_params%nzv)

        stop

c     End program

      contains

c     gridInfo
c     #################################################################  
      subroutine gridInfo(char,ig,istart,nxv,xx,dx,dxh)

        implicit  none

c     Call variables

        character*(1) :: char
        integer(4) :: ig,istart(*),nxv(*)
        real(8)    :: xx(*),dx(*),dxh(*)

c     Local variables

        integer(4) :: i,j,isig

c     Begin program

        write (*,*)
        write (*,*) '***************************'
        write (*,*) 'MG grid in ',char,'-axis'
        write (*,*) '***************************'
        do i = 1,ig
          isig = istart(i)
          write (*,*)
          write (*,*) '************* Grid ',i,' **************'
          write (*,*)
          write (*,*) 'Size ',nxv(i)
          write (*,*) 'MG pointer: ',isig
          write (*,*) 'Grid nodes'
          do j = isig,isig+nxv(i)+1
            write (*,10) 'node :',j-isig
     .                  ,'   position: ',xx(j)
     .                  ,'   int dh: ',dx(j)
     .                  ,'   half dh: ',dxh(j)
 10         format (a,i3,a,f6.3,a,f6.3,a,f6.3)
          enddo
        enddo

c     End program

      end subroutine gridInfo

c     hessianCheck
c     #################################################################
      subroutine hessianCheck(igx,igy,igz,nxv,nyv,nzv)

        implicit  none

c     Call variables

        integer(4) :: igx,igy,igz,nxv(*),nyv(*),nzv(*)

c     Local variables

        integer(4) :: i,j,k,i1,j1,k1,ig,jg,kg,ih
        real(8)    :: x1,y1,z1
        real(8)    :: hess(3,3,3),hess_cnv(3,3,3),table(3,3,3)

c     Begin program

        write (*,*) 
        write (*,*) '*************************'
        write (*,*) '     Hessian check'
        write (*,*) '*************************'
        write (*,*)

        do i = 2,nxv(igx)   !Start at 2 to avoid singular points
          do j = 1,nyv(igy)
            do k = 1,nzv(igz)

cc              call getCoordinates(i,j,k,igx,igy,igz,x1,y1,z1,cartesian)

              call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)

              x1 = grid_params%xx(ig)
              y1 = grid_params%yy(jg)
              z1 = grid_params%zz(kg)

              hess(1,:,:) = hessian(1,x1,y1,z1,.false.)
              hess(2,:,:) = hessian(2,x1,y1,z1,.false.)
              hess(3,:,:) = hessian(3,x1,y1,z1,.false.)

              hess_cnv(1,:,:) = hessian_cnv(1,x1,y1,z1)
              hess_cnv(2,:,:) = hessian_cnv(2,x1,y1,z1)
              hess_cnv(3,:,:) = hessian_cnv(3,x1,y1,z1)

              do i1=1,3
                do j1=1,3
                  do k1=1,3
                    table(i1,j1,k1) = hess_cnv(i1,j1,k1)
     .                               -delta(i1,k1)*(hess(1,j1,1)
     .                                             +hess(2,j1,2)
     .                                             +hess(3,j1,3))
     .                               +hess(k1,i1,j1)
                  enddo
                enddo
              enddo

              write (*,5) 'Grid point: (',x1,',',y1,',',z1,')'
 5            format (/,a,f7.3,a,f7.3,a,f7.3,a)

              write (*,*) 
              write (*,*) 'Hessian relation'

              do ih=1,3
                write (*,*) 
                write (*,10) table(ih,1,1:3)
                write (*,10) table(ih,2,1:3)
                write (*,10) table(ih,3,1:3)
 10             format (3f10.3)
              enddo

            enddo
          enddo
        enddo

c     End program

      end subroutine hessianCheck

c     delta
c     #################################################################
      real(8) function delta(i,j)

        integer(4) :: i,j

        delta = 0d0
        if (i == j) delta = 1d0

      end function delta

c     metricTensorCheck
c     #################################################################
      subroutine metricTensorCheck(igx,igy,igz,nxv,nyv,nzv)

        implicit  none

c     Call variables

        integer(4) :: igx,igy,igz,nxv(*),nyv(*),nzv(*)

c     Local variables

        integer(4) :: i,j,k,i1,j1,k1,ig,jg,kg,ih
        real(8)    :: x1,y1,z1,check
        real(8)    :: gup(3,3),gdown(3,3),tensor(3,3)
        logical    :: cartesian

c     Begin program

        check = 0d0

        write (*,*) 
        write (*,*) '*************************'
        write (*,*) '   Metric Tensor check'
        write (*,*) '*************************'
        write (*,*)

        do i = 2,nxv(igx)   !Start at 2 to avoid singular points
          do j = 1,nyv(igy)
            do k = 1,nzv(igz)

              call getCoordinates(i,j,k,igx,igy,igz,ig,jg,kg,x1,y1,z1
     .                           ,cartesian)

cc              call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)
cc
cc              x1 = grid_params%xx(ig)
cc              y1 = grid_params%yy(jg)
cc              z1 = grid_params%zz(kg)

              gup = g_super(x1,y1,z1,cartesian)
              gdown = g_sub(x1,y1,z1,cartesian)

              tensor = matmul(gup,gdown)

              write (*,5) 'Grid point: (',x1,',',y1,',',z1,')'
 5            format (/,a,f7.3,a,f7.3,a,f7.3,a)

              write (*,*) 
              write (*,*) 'Metric tensor product'
              write (*,10) tensor(1,1:3)
              write (*,10) tensor(2,1:3)
              write (*,10) tensor(3,1:3)
 10           format (3f10.3)

            enddo
          enddo
        enddo

c     End program

      end subroutine metricTensorCheck

      end subroutine checkGrid

      end module grid

