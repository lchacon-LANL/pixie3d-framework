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
c     Gets grid packing config
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
c     Sets up grid packing method
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

      integer :: nn(3),dim
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

          if ((.not.bcPER(dim)).and.(.not.g_pack(dim)%pack_both_bdrys))
     .       eta(dim)=find_eta(lmin(dim),lmax(dim)
     .                        ,g_pack(dim)%xp,g_pack(dim)%slope)
        endif
      enddo

      contains

      end subroutine find_pack_slope
      
c     find_eta
c     #################################################################
      function find_eta(xmin,xmax,x0,slope) result(eta)

        implicit none

c     Global variables

        real(8) :: xmin,xmax,x0,slope,eta

c     Local variables

        integer,parameter :: nit=50,neq=4

        integer :: i
        real(8) :: jac(neq,neq),ijac,dxi(nit,neq),res(nit,neq),xi(neq)
        real(8) :: LL,Lm,Lp,alpha
        
        real(8) :: dd,xbar,dx,eps,a1,a2
cc        real(8) :: lambda,ilambda,Tm,Tp,Cm,Cp,Ep,Em
        real(8) :: erfv,exp1,exp2,ff,dff,dexp1d1,dexp1d2,dexp2d1
     .            ,dexp2d2,derfv,derfv1,derfv2,erfv1,erfv2

c     Begin program

        eps = 10*findRoundOff()

        if (slope == 0d0 .or. slope == 1d0) then
          eta = 0d0
          return
        endif

        LL = xmax-xmin
        Lp = xmax-x0
        Lm = LL - Lp

        res = 0d0 ; dxi = 0d0

        if (slope > 1d0 .or. Lp == 0d0 .or. Lm == 0d0) then
          call prep_cubic
        else
          call prep_gaussian
        endif
          
c       Convergence check

        if (i >= nit) then
          write (*,*)
          write (*,*) 'Nonlinear it. in grid-packing routine find_eta'
     .               ,' did not converge'
          write (*,*)
          write (*,*) 'Convergence history:'
          do i=1,nit
            write (*,*) 'Nonlinear it.',i,' residual'
     .                 ,sqrt(sum(res(i,:)**2))
     .                 ,'update',sqrt(sum(dxi(i,:)**2))
          enddo
          write (*,*)
          write (*,*) 'Aborting...'
          stop
        endif

      contains

c     prep_cubic
c     ###############################################################
      subroutine prep_cubic
      
      implicit none

c     ---------------------------------------------------------------
c     Prepares cubic polynomial for grid packing
c     ---------------------------------------------------------------

      i = 0

      if (slope > 1d0) then !Contraction
        
cc          eta = LL**3*(1d0-1d0/slope)/((xmax-x0)**3-(xmin-x0)**3)
        eta = LL**pwr*(1d0-1d0/slope)/((xmax-x0)**pwr-(xmin-x0)**pwr)

      else
                
        xbar= 0.5*(xmax+xmin)
        dx  = xbar - x0
        dd  = (dx/LL)**2

        !Initial guess

        eta = 4d0*(1d0-slope)

        !Newton iteration

        if (dx == 0d0) return

        do i=1,nit
          res(i,1)=(1+2*(eta+slope))**2*(1-slope-0.25*eta) -27.*dd*eta
          ijac = 4*(1+2*(eta+slope)) *(1-slope-0.25*eta)
     .       -0.25*(1+2*(eta+slope))**2                    -27.*dd

          dxi(i,1) = -res(i,1)/ijac
          eta = abs(eta + dxi(i,1)) !Want to find positive root

cc        write (*,*) abs(res(i,1)),dxi(i,1),ijac,eta

          if (abs(res(i,1))<1d-10.or.(abs(dxi(i,1))<eps)) exit
        enddo

      endif
          
      end subroutine prep_cubic

c     prep_gaussian
c     ###############################################################
      subroutine prep_gaussian

      implicit none

c     ---------------------------------------------------------------
c     Prepares Gaussian profile for grid packing
c     ---------------------------------------------------------------

      if (Lp == Lm) then
        Lp = Lp*(1d0+1d-8)
        Lm = Lm*(1d0-1d-8)
      endif

      !Initial condition (assume alpha ~ 0.5 + delta, delta << 1 and solve)
      eta = 1d1
      do i=1,5
        eta =2*sqrt(abs(log(0.5*LL*slope/(2*eta*Lp/sqrt(pi)))))
      enddo
      alpha = 0.5*(1-log(Lm/Lp)/eta/eta)
      a1 = sqrt(4/pi)*eta*Lm/LL
      a2 = sqrt(4/pi)*eta*Lp/LL

      !Newton iteration (variable order: a1,a2,alpha=Ep/L,eta=L/sqrt(sigma))
      do i=1,nit
        erfv1 = sqrt(pi/4.)*erf((1-alpha)*eta)
        erfv2 = sqrt(pi/4.)*erf(   alpha *eta)
        erfv  = sqrt(pi/4.)*erf(          eta)

        derfv1 = exp(-((1-alpha)*eta)**2)
        derfv2 = exp(-((  alpha)*eta)**2)
        derfv  = exp(-eta**2)
        
        exp1    = exp(eta**2*(1-2*alpha))
        dexp1d1 =-2*eta*eta*exp1         !wrt alpha
        dexp1d2 = 2*eta*(1-2*alpha)*exp1 !wrt eta
        
        exp2    = exp(-(alpha*eta)**2)
        dexp2d1 = -2*alpha*eta*eta*exp2  !wrt alpha
        dexp2d2 = -2*eta*alpha*alpha*exp2!wrt eta

        ff  = eta/erfv
        dff = 1./erfv - eta*derfv/(erfv*erfv)
        
        !Residual
        res(i,1) = Lm/LL*eta - a1*erfv1 - a2*erfv + a2*erfv2
        res(i,2) = a1 + a2 - ff
        res(i,3) = a1 - a2*alpha*exp1/(1-alpha)
        res(i,4) = slope - a2*exp2/(1-alpha)

        !Jacobian solution
        jac(1,1) = -erfv1
        jac(1,2) = +erfv2-erfv
        jac(1,3) = a1*eta*derfv1 + a2*eta*derfv2
        jac(1,4) = Lm/LL
     .           - a1*(1-alpha)*derfv1
     .           - a2*derfv
     .           + a2*alpha*derfv2

        jac(2,1) = 1.
        jac(2,2) = 1.
        jac(2,3) = 0d0
        jac(2,4) = -dff

        jac(3,1) = 1.
        jac(3,2) = -alpha*exp1/(1-alpha)
        jac(3,3) = -a2*(1-alpha)**(-2)*exp1
     .             -a2*alpha/(1-alpha)*dexp1d1
        jac(3,4) = -a2*alpha*dexp1d2/(1-alpha)

        jac(4,1) = 0.
        jac(4,2) = -exp2/(1-alpha)
        jac(4,3) = -a2*(1-alpha)**(-2)*exp2 - a2*dexp2d1/(1-alpha)
        jac(4,4) = -a2*dexp2d2/(1-alpha)

        call blockSolve(neq,jac,1,res(i,:),dxi(i,:))

        !Update solution
        a1    = abs(   a1 - dxi(i,1))  !Enforce a1 > 0
        a2    = abs(   a2 - dxi(i,2))  !Enforce a1 > 0
        alpha = abs(alpha - dxi(i,3))  !Enforce alpha > 0
        eta   = abs(  eta - dxi(i,4))  !Enforce eta >= 0d0

cc        write (*,*) i,a1,a2,alpha,eta
        
        !Check convergence
cc        write (*,*) i,dxi(i,:),res(i,:)

        if (  sqrt(sum(res(i,:)**2))<1d-10
     .   .or.(sqrt(sum(dxi(i,:)**2))<eps)) exit
        
      enddo

      end subroutine prep_gaussian

c$$$c     prep_tanh
c$$$c     ###############################################################
c$$$      subroutine prep_tanh
c$$$
c$$$      implicit none
c$$$              
c$$$          LL = xmax-xmin
c$$$          Lp = xmax-x0
c$$$          Lm = LL - Lp
c$$$
c$$$          if (Lp == Lm) then
c$$$            Lp = Lp*(1d0+1d-8)
c$$$            Lm = Lm*(1d0-1d-8)
c$$$          endif
c$$$
c$$$          !Newton iteration
c$$$          lambda = 1d0 ; Ep = Lp
c$$$
c$$$          do i=1,nit
c$$$            Em = LL - Ep
c$$$            ilambda = 1d0/lambda
c$$$
c$$$            Cp = 1d0/cosh(Ep*ilambda)
c$$$            Cp = Cp*Cp
c$$$
c$$$            Cm = 1d0/cosh(Em*ilambda)
c$$$            Cm = Cm*Cm
c$$$
c$$$            Tp = tanh(Ep*ilambda)
c$$$            Tm = tanh(Em*ilambda)
c$$$
c$$$            !Var 1 is Ep, Var 2 is lambda
c$$$            res(i,1) = (Ep + slope*lambda)*Tp - Lp
c$$$            res(i,2) = (Em + slope*lambda)*Tm - Lm
c$$$
c$$$            jac(1,1) = Tp + (Ep+slope*lambda)*Cp*ilambda
c$$$            jac(1,2) =slope*Tp-(Ep+slope*lambda)*Cp*Ep*(ilambda*ilambda)
c$$$            jac(2,1) =-Tm - (Em+slope*lambda)*Cm*ilambda
c$$$            jac(2,2) =slope*Tm-(Em+slope*lambda)*Cm*Em*(ilambda*ilambda)
c$$$
c$$$            call blockSolve(2,jac,1,res(i,:),dxi(i,:))
c$$$
c$$$            Ep = abs(Ep - dxi(i,1))
c$$$            lambda = abs(lambda - dxi(i,2))
c$$$
c$$$cc            write (*,*) i,Ep,lambda,res(i,:),dxi(i,:)
c$$$            if (  sqrt(sum(res(i,:)**2))<1d-10
c$$$     .       .or.(sqrt(sum(dxi(i,:)**2))<eps)) exit
c$$$          enddo
c$$$
c$$$          eta = Ep/lambda
c$$$
c$$$      end subroutine prep_tanh

      end function find_eta

c     scale_coord
c     #################################################################
      function scale_coord(gdef,dim,x1,deriv) result(xlog)

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

        integer,optional :: deriv
        
c     Local variables

        integer :: nn,der

        real(8) :: gp(2),lmin,lmax,TT,bb

        logical :: pack_both_bdrys

c     Begin program

        if (PRESENT(deriv)) then
          der = deriv
        else
          der = 0
        endif
        
        if (.not.gdef%g_pack(dim)%pack) then
          select case(der)
          case(0)
            xlog = x1
          case(1)
            xlog = 1d0
          case(2)
            xlog = 0d0
          end select
          return
        endif

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

        if (bcPER(dim)) then !cosine jacobian
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
        else
          if (gp(2)>1d0.or.(lmin==gp(1)).or.(lmax==gp(1))) then
            xlog = cubic_scl   (dim,x1,eta(dim),lmin,lmax,gp(1),gp(2)
     .                         ,deriv=deriv)
          else
            xlog = gaussian_scl(dim,x1,eta(dim),lmin,lmax,gp(1),gp(2)
     .                         ,deriv=deriv)
          endif
        endif

      contains

c     cubic_scl
c     #################################################################
      function cubic_scl(dim,xi,eta,xmin,xmax,x0,lambda,deriv) result(x)

c     -----------------------------------------------------------------
c     Packs coordinate according to cubic polynomial at x0 (location)
c     with slope lambda.  Note:
c       * lambda <> 1d0 --> dx0 ~ dxi*dlambda
c       * lambda == 0d0 --> No packing
c     -----------------------------------------------------------------

        implicit none

c     Call variables

        integer :: dim
        real(8),INTENT(IN) :: eta,xmin,xmax,x0,lambda,xi
        real(8) :: x

        integer,optional :: deriv
        
c     Local variables

        real(8) :: xi0,xibar,LL,bb,dx

        integer,parameter :: nit=100

        integer :: i,der
        real(8) :: ijac,dxi(nit),res(nit)

c     Begin program

        if (PRESENT(deriv)) then
          der = deriv
        else
          der = 0
        endif
        
        if (eta == 0d0) then

          select case(der)
          case(0)
            x = xi
          case(1)
            x = 1d0
          case(2)
            x = 0d0
          end select
          return

        else

          if (der > 0)
     .         call sstop(0,'cubic_scl','Derivative not implemented')
          
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
c$$$            if (xi < xmin.and.bcSP(dim=dim)) then
c$$$              dx = xmin - xi
c$$$              x =  x0 + lambda*(xmin+dx-xi0) + bb*(xmin+dx-xi0)**3
c$$$              x = 2*xmin - x
c$$$            else
c$$$              x = x0 + lambda*(xi-xi0) + bb*(xi-xi0)**3
c$$$            endif

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
     .                   ' cubic_scl did not converge'
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

      end function cubic_scl

c$$$c     tanh_scl
c$$$c     #################################################################
c$$$      function tanh_scl(dim,xi,eta,xmin,xmax,x0,x0p) result(x)
c$$$
c$$$c     -----------------------------------------------------------------
c$$$c     Packs coordinate using |x|*tanh(x/lambda) prescription at location
c$$$c     x0 with slope x0p. Note:      
c$$$c       * x0p <> 1d0 --> dx0 ~ dxi*x0p
c$$$c       * x0p == 0d0 --> No packing
c$$$c     -----------------------------------------------------------------
c$$$
c$$$        implicit none
c$$$
c$$$c     Global variables
c$$$
c$$$        integer :: dim
c$$$        real(8),INTENT(IN) :: xmin,xmax,x0,eta,xi,x0p
c$$$        real(8) :: x
c$$$
c$$$c     Local variables
c$$$
c$$$        real(8) :: bb,LL,Lp,Lm,Tp,Tm,lambda,xi0,Ep,dx
c$$$
c$$$c     Begin program
c$$$
c$$$        if (x0p == 0d0 .or. eta == 0d0) then
c$$$
c$$$          x = xi
c$$$
c$$$        else
c$$$
c$$$          LL = xmax-xmin
c$$$          Lp = xmax-x0
c$$$          Lm = LL - Lp
c$$$
c$$$          if (Lp == Lm) then
c$$$            Lp = Lp*(1d0+1d-8)
c$$$            Lm = Lm*(1d0-1d-8)
c$$$          endif
c$$$
c$$$          lambda = Lp/((eta+x0p)*tanh(eta))
c$$$
c$$$          Ep = eta*lambda
c$$$
c$$$          xi0 = xmax - Ep
c$$$          bb  = x0p*lambda
c$$$
c$$$          if (xi < xmin.and.bcSP(dim=dim)) then
c$$$            dx = xmin - xi
c$$$            x = x0 + (abs(xmin+dx-xi0) + bb)*tanh((xmin+dx-xi0)/lambda)
c$$$            x = 2*xmin - x
c$$$          else
c$$$            x = x0 + (abs(xi-xi0) + bb)*tanh((xi-xi0)/lambda)
c$$$          endif
c$$$
c$$$        endif
c$$$
c$$$c     End program
c$$$
c$$$      end function tanh_scl

c     gaussian_scl
c     #################################################################
      function gaussian_scl(dim,xi,eta,xmin,xmax,x0,x0p,deriv) result(x)

c     -----------------------------------------------------------------
c     Packs coordinate using two-gaussian prescription at location
c     x0 with slope x0p. Gaussians at centered at boundaries, to
c     ensure that no discontinuity is introduded at ghost cells. Note:      
c       * x0p <> 1d0 --> dx0 ~ dxi*x0p
c       * x0p == 0d0 --> No packing
c     -----------------------------------------------------------------

        implicit none

c     Global variables

        integer :: dim
        real(8),INTENT(IN) :: xmin,xmax,x0,eta,xi,x0p
        real(8) :: x

        integer,optional :: deriv
        
c     Local variables

        integer :: i,der
        
        real(8) :: LL,Lp,Lm,ssigma,alpha,a1,a2,erfv,exp1,exp2,dexp1d1
     .            ,dexp2d1,ff,res,ijac,d1,d2,d1o,d2o

c     Begin program

        if (PRESENT(deriv)) then
          der = deriv
        else
          der = 0
        endif
        
        if (x0p == 0d0 .or. eta == 0d0) then

          select case(der)
          case(0)
            x = xi
          case(1)
            x = 1d0
          case(2)
            x = 0d0
          end select
          
        else

          LL = xmax-xmin
          Lp = xmax-x0
          Lm = LL - Lp

          if (Lp == Lm) then
            Lp = Lp*(1d0+1d-8)
            Lm = Lm*(1d0-1d-8)
          endif

          ssigma = LL/eta

          !solve for alpha, a1, a2
          alpha = 0.5*(1-log(Lm/Lp)/eta/eta) !Initial condition
          do i=1,20
            erfv    = sqrt(pi/4.)*erf(eta)

            exp1    = exp(eta**2*(1-alpha)**2)
            dexp1d1 =-2*(1-alpha)*eta*eta*exp1

            exp2    = exp(eta**2*alpha**2)
            dexp2d1 = 2*alpha*eta*eta*exp2
            
            ff  = eta/erfv

            res = x0p*((1-alpha)*exp2 + alpha*exp1) - ff

            ijac = x0p*(-exp2 + (1-alpha)*dexp2d1
     .                  +exp1 +    alpha *dexp1d1)
            ijac = 1d0/ijac
            
            alpha = alpha - ijac*res

cc            write (*,*) i,"alpha",alpha,"dalpha=",ijac*res
            if (abs(ijac*res) < 1d-12) exit
          enddo

          a2 = x0p*(1-alpha)*exp2
          a1 = a2*alpha/(1-alpha)*exp((1-2*alpha)*eta**2)

          !Find corrections for zero-derivative at domain ends
          d1 = 0d0 ; d2 = 0d0
          do i=1,10
            d1o = d1 ; d2o = d2
            d1 = -a2/a1*(LL+d2)*exp(-((LL+d2)**2-d1**2)/ssigma**2)
            d2 =  a1/a2*(LL-d1)*exp(-((LL-d1)**2-d2**2)/ssigma**2)
            if (abs(d1o-d1) < 1d-10.and.abs(d2o-d2) < 1d-10) exit
          enddo

          !Enforce domain size
c$$$          a1 = a1*erf(LL/ssigma)/(erf((LL-d1)/ssigma)-erf(-d1/ssigma))
c$$$          a2 = a2*erf(LL/ssigma)/(erf((LL+d2)/ssigma)-erf( d2/ssigma))
c$$$          write(*,*) "Gaussian corrections=",d1,d2
c$$$     . ,-a1*2*(xmin-xmin-d1)*(ssigma)**(-2)
c$$$     .         *exp(-((xmin-xmin-d1)/ssigma)**2)
c$$$     .  -a2*2*(xmin-xmax-d2)*(ssigma)**(-2)
c$$$     .         *exp(-((xmin-xmax-d2)/ssigma)**2)  !Derivative at xmin (=0)
c$$$     . ,-a1*2*(xmax-xmin-d1)*(ssigma)**(-2)
c$$$     .         *exp(-((xmax-xmin-d1)/ssigma)**2)
c$$$     .  -a2*2*(xmax-xmax-d2)*(ssigma)**(-2)
c$$$     .         *exp(-((xmax-xmax-d2)/ssigma)**2)  !Derivative at xmax (=0)
c$$$     . ,xmin + 0.5*sqrt(pi)*ssigma*(
c$$$     .                a1*(erf((xmax-xmin-d1)/ssigma)-erf(-d1/ssigma))
c$$$     .               +a2*(erf((LL+d2)/ssigma)-erf((d2)/ssigma)))  !xmax (=1)
c$$$           write(*,*)
          
          !Find coordinate
          select case(der)
          case(0)
            x = xmin + 0.5*sqrt(pi)*ssigma*(
     .                a1*(erf((xi-xmin-d1)/ssigma)-erf(-d1/ssigma))
     .               +a2*(erf((LL+d2)/ssigma)-erf((xmax-xi+d2)/ssigma)))
          case(1)
            x = a1*exp(-((xi-xmin-d1)/ssigma)**2)
     .         +a2*exp(-((xi-xmax-d2)/ssigma)**2)
          case(2)
            x =
     .  -a1*2*(xi-xmin-d1)*(ssigma)**(-2)*exp(-((xi-xmin-d1)/ssigma)**2)
     .  -a2*2*(xi-xmax-d2)*(ssigma)**(-2)*exp(-((xi-xmax-d2)/ssigma)**2)
          end select
        endif

c     End program

      end function gaussian_scl

      end function scale_coord

c     findGlobalUniformGrid
c     #################################################################
      subroutine findGlobalUniformGrid(lmin,lmax,bcs1,bcs2,xg)

c     -----------------------------------------------------------------
c     Finds global grid xg(0:ng+1) in domain (lmin,lmax) with boundary
c     conditions bcs1, bcs2, according to specifications in dpack.
c     The input variable crse_ratio is used to adjust the packing provided
c     in dpack for the fine mesh to coarser grid levels.
c     -----------------------------------------------------------------

        implicit none

c     Call variables

        integer :: bcs1,bcs2

        real(8),INTENT(IN) :: lmin,lmax

        real(8),INTENT(OUT) :: xg(0:)

c     Local variables

        integer :: ng,i
        real(8) :: length,lstart,dh

c     Begin program

        ng = size(xg)-2

        length = lmax-lmin   !Global length, based on global limits

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

c     End program

      end subroutine findGlobalUniformGrid

c     findGlobalPackedGrid
c     #################################################################
      subroutine findGlobalPackedGrid(lmin,lmax,bcs1,bcs2,xp,dx0,dxl,dxr
     .                               ,xg,xi,xil,xir)

c     -----------------------------------------------------------------
c     Finds global grid xg(0:ng+1) in domain (lmin,lmax) with boundary
c     conditions bcs1, bcs2, according to specifications in dpack.
c     The input variable crse_ratio is used to adjust the packing provided
c     in dpack for the fine mesh to coarser grid levels.
c     -----------------------------------------------------------------

        implicit none

c     Call variables

        integer :: bcs1,bcs2

        real(8),INTENT(IN) :: lmin,lmax,dx0,dxl,dxr,xp

        real(8),INTENT(OUT) :: xg(0:),xi,xil,xir

c     Local variables

        integer :: n1,n2,ng,i
        real(8) :: dh,length,lstart,lend
        logical :: pack_l,pack_r,fv_l,fv_r

c     Begin program

        ng = size(xg)-2

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

        dh = length/ng
        lstart = lmin
        lend   = lmax

        pack_l = .false.
        pack_r = .false.

        if (xp < lend .and. xp > lstart) then

          n1 = max(int((xp - lstart)/length*ng),2)

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

          xi = 0d0

        else

          if (xp == lstart) then !Pack to the left
            pack_l = .true.
            if (fv_l) then      !Ensure left boundary is finite-volume type
              call packGrid(ng,lstart+0.5*dx0,lend,dx0,xg(1:ng),xi,fv_r)
            else
              call packGrid(ng,lstart        ,lend,dx0,xg(1:ng),xi,fv_r)
            endif
          else                  !Pack to the right
            pack_r = .true.
            if (fv_r) then      !Ensure right boundary is finite-volume type
             call packGrid(ng,lstart,lend-0.5*dx0,-dx0,xg(1:ng),xi,fv_l)
            else
             call packGrid(ng,lstart,lend        ,-dx0,xg(1:ng),xi,fv_l)
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
          if (pack_l) then
            xg(0) = xg(1) - (xg(2)-xg(1))/xil
          else
            xg(0) = xg(1) - (xg(2)-xg(1))*xil
          endif
        endif

        !Xn ghost cell
        if (bcs2 == PER) then
          xg(ng+1) = xg(ng) + (xg(2) -xg(1)   )
        elseif (abs(bcs2) == SYM) then
          xg(ng+1) = xg(ng) + (xg(ng)-xg(ng-1))
        else
          if (pack_r) then
            xg(ng+1) = xg(ng) + (xg(ng)-xg(ng-1))/xir
          else
            xg(ng+1) = xg(ng) + (xg(ng)-xg(ng-1))*xir
          endif
        endif

c     End program

      end subroutine findGlobalPackedGrid

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

c     are_grid_packing
c     #################################################################
      function are_grid_packing(gdef) result(doing_grid_pack)

c     -----------------------------------------------------------------
c     Sets up grid packing method
c     -----------------------------------------------------------------

      implicit none

      type(grid_mg_def),pointer :: gdef

      logical :: doing_grid_pack

c     Begin program

        doing_grid_pack =  (    gdef%g_pack(1)%pack
     .                      .or.gdef%g_pack(2)%pack
     .                      .or.gdef%g_pack(3)%pack)

c     End program

      end function are_grid_packing

      end module grid_packing
