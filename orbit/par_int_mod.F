c module kernels
c #############################################################################
      module kernels

        use orbit

        integer :: kmax=100      !Maximum number of Fourier modes in Fourier Green's f. form
        real(8) :: t_tol,T_cutoff2,S_cutoff2

      contains

c     tkernel_3d
c     ###############################################################
      function tkernel_3d(ss,ct0,alpha,ex_loop) result(tkernel)

c     ---------------------------------------------------------------
c     3D temperature kernel integral
c     ---------------------------------------------------------------

      implicit none

c     Call variables

      real(8) :: tkernel,ss,ct0,alpha
      logical :: ex_loop

c     Local variables

c     Begin program

      if (ct0 == 0d0) return

      if (alpha == 2d0) then
        tkernel = 0.5/sqrt(pi*ct0)*exp(-0.25*ss**2/ct0) !ct0=chi_par*t0

cc        ex_loop = (abs(ss) > sqrt(s_max**2 + 4*T_cutoff2*ct0))
        ex_loop = (abs(ss) > sqrt(4*T_cutoff2*ct0))
      elseif (alpha == 1d0) then
        tkernel = 1d0/(pi*ct0*(1+(abs(ss)/ct0)**2))

cc        ex_loop = (abs(ss) >  s_max + 2/pi/t_tol*ct0)
        ex_loop = (abs(ss) >  2d0/pi/t_tol*ct0)
      else
        call pstop('tkernel_3d','Fractional order not implemented')
      endif

c     End program

      end function tkernel_3d

c     tkernel_per
c     ###############################################################
      function tkernel_per(ss,ct0,alpha) result(tkernel)

c     ---------------------------------------------------------------
c     2D temperature kernel integral
c     ---------------------------------------------------------------

      implicit none

c     Call variables

      real(8) :: tkernel,ss,ct0,alpha

c     Local variables

      integer :: kk,k0
      real(8) :: twopi_L,dtk,eta,smax,ismax
      logical :: ex_loop,ex_loop2

c     Begin program

      if (ct0 == 0d0) return

      smax = orbit_conf%ode_conf%smax
      ismax = 1d0/smax

      twopi_L = 2d0*pi*ismax

      eta = ct0*(twopi_L)**alpha

      if (eta < 1d0) then
cc        if (my_rank == 0) write (*,*) 'tkernel_per: small eta'

        tkernel = tkernel_3d(ss,ct0,alpha,ex_loop)

        kk = 1
        do 
          tkernel = tkernel + tkernel_3d(ss-kk*smax,ct0,alpha,ex_loop)
     .                      + tkernel_3d(ss+kk*smax,ct0,alpha,ex_loop2)

          if (ex_loop .and. ex_loop2) exit
          kk = kk + 1
        enddo

      else
cc        if (my_rank == 0) write (*,*) 'tkernel_per: large eta'

        !Compute iteration limit
        k0 = ceiling((abs(log(0.5*t_tol))/eta)**(1d0/alpha))

        !Accumulate periodic Green's function
        tkernel = 1d0

        do kk=1,k0
          dtk = exp(-eta*kk**alpha)
          tkernel = tkernel + 2d0*dtk*cos(twopi_L*kk*ss)
        enddo

        tkernel = tkernel*ismax

      endif

c     End program

      end function tkernel_per

c     skernel_3d
c     ###############################################################
      function skernel_3d(ss,ct0,alpha,ex_loop) result(skernel)

c     ---------------------------------------------------------------
c     3D source kernel integral
c     ---------------------------------------------------------------

      implicit none

c     Call variables

      real(8) :: skernel,ss,ct0,alpha
      logical :: ex_loop

c     Local variables

      real(8) :: eta,z

#if defined(pgf90)
      real(8) :: erfc
#endif

c     Begin program

      if (ct0 == 0d0) return

      z = abs(ss)

      if (alpha == 2d0) then
        eta = 0.5*z/sqrt(ct0) !ct0=chi_par/chi_perp*t0
        skernel = (exp(-eta**2)/sqrt(pi)-eta*erfc(eta))/sqrt(ct0)

cc        ex_loop = (z > sqrt(s_max**2 + 4*S_cutoff2*ct0))
        ex_loop = (z > sqrt(4*S_cutoff2*ct0))
cc      elseif (alpha == 1d0) then
cc        skernel = 1d0/(pi*ct0*(1+(abs(ss)/ct0)**2))
cc
cc        slim = s_max + 2/pi/t_tol*ct0
cc
cc        ex_loop = (abs(ss) > slim)
      else
        call pstop('skernel_3d','Fractional order not implemented')
      endif

c     End program

      end function skernel_3d

c     skernel_per
c     ###############################################################
      function skernel_per(ss,ct0,alpha) result(skernel)

c     ---------------------------------------------------------------
c     2D source kernel integral
c     ---------------------------------------------------------------

      implicit none

c     Call variables

      real(8) :: skernel,ss,ct0,alpha

c     Local variables

      integer :: kk,k0
      real(8) :: twopi_L,dtk,eta,smax,ismax
      logical :: ex_loop,ex_loop2

c     Begin program

      if (ct0 == 0d0) return

      smax = orbit_conf%ode_conf%smax
      ismax = 1d0/smax

      twopi_L = 2d0*pi*ismax

      eta = ct0*(twopi_L)**alpha

      if (eta < 1d0/sqrt(t_tol*abs(log(t_tol)))) then

cc        if (my_rank == 0) write (*,*) 'skernel_per: small eta'

        skernel = skernel_3d(ss,ct0,alpha,ex_loop)

        kk = 1
        do
          skernel = skernel + skernel_3d(ss-kk*smax,ct0,alpha,ex_loop)
     .                      + skernel_3d(ss+kk*smax,ct0,alpha,ex_loop2)

          if (ex_loop .and. ex_loop2) exit
          kk = kk + 1
        enddo

      else

cc        if (my_rank == 0) write (*,*) 'skernel_per: large eta'

        !Compute iteration limit
        k0 = ceiling((2d0/t_tol/eta)**(1d0/alpha))
        k0 = min(k0,kmax) !To avoid aliasing errors

        !Accumulate periodic Green's function
        skernel = 1d0

        do kk=1,k0
          eta = ct0*(twopi_L*kk)**alpha
          dtk = (1d0-exp(-eta))/eta

          skernel = skernel + 2d0*dtk*cos(twopi_L*kk*ss)
        enddo

        skernel = skernel*ismax

      endif

c     End program

      end function skernel_per

      end module kernels

c module par_int
c #############################################################################
      module par_int

        use orbit

        use grid

        use timestepping,ONLY: dt,dtold,cnp,cn,cnm,bdf2,time,dtbase

        use kernels

        !Configuration
        logical :: bwgt
        real(8) :: alpha

        !Field line accumulators
cc        integer :: n_lag_ints=0        !Number of lagrangian integrals
        real(8) :: ct0,ct1  !chi*dt
        real(8) :: temp_s,kernel_s,length_s,b0_s

        !Private variables
cc        real(8),private :: orbit_atol=1e-10    !Orbit integrator absolute tolerance
cc     .                    ,orbit_rtol=1e-8     !Orbit integrator relative tolerance
        real(8),private :: orbit_atol=1d-12    !Orbit integrator absolute tolerance
     .                    ,orbit_rtol=1d-12    !Orbit integrator relative tolerance
     .                    ,orbit_dt  =1d-2     !Orbit integrator time step
     .                    ,CHI_LIMIT =1d-4     !Lower bound for chi*dt (delta-function limit)

        real(8),private :: lsode_out(LSODE_ACC_DIM),dx0,dy0,dz0

        logical,private :: init_lag=.false.,l2d,x_per,y_per,z_per
     .                    ,solenoidal

        real(8),pointer,dimension(:,:,:) :: ll,b0

        real(8),pointer,dimension(:,:,:,:) :: sp_coef

        private :: lag_2d,lag_3d,lag_2d_init,lag_3d_init

      contains

c     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c     PRIVATE ROUTINES
c     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ccc     lagrangian_2d_scl
ccc     ##############################################################
cc      subroutine lagrangian_2d_scl(igrid,nx,ny,nz,tt
cc     .                            ,chi,iout,int_f,lsode_f)
cc
cc      implicit none
cc
ccc     --------------------------------------------------------------
ccc     Performs Lagrangian temperature integration according to B
ccc     field (bx,by,bz).
ccc     --------------------------------------------------------------
cc
ccc     Call variables
cc
cc        integer :: nx,ny,nz,igrid,iout
cc
cc        real(8) :: chi
cc
cc        real(8) :: tt (0:nx+1,0:ny+1,0:nz+1)
cc
ccc     Local variables
cc
cc        integer :: i,j,k,ig,jg,kg,nxg,nyg,nzg
cc
cc        real(8) :: x0,y0,z0
cc
cc        real(8),pointer,dimension(:,:,:) :: tg,zeros
cc
cc        logical :: have_lsode_f
cc
ccc     Interfaces
cc
cc        optional :: int_f
cc
cc        INTERFACE
cc          function int_f(fwd,xo,x,yo,y,zo,z,so,ss)
cc          logical :: int_f,fwd
cc          real(8) :: xo,x,yo,y,zo,z,so,ss
cc          end function int_f
cc        END INTERFACE
cc
cc        optional :: lsode_f
cc
cc        !Interface for optional LSODE function evaluation
cc        INTERFACE
cc          subroutine lsode_f(neq,t,y,ydot)
cc          integer :: neq
cc          real(8) :: t,y(neq),ydot(neq)
cc          end subroutine lsode_f
cc        END INTERFACE
cc
ccc     Begin program
cc
cc        if (maxval(tt) == minval(tt)) return
cc
cc        if (chi==0d0) return
cc
cc        ct0 = chi*dt
cc
cc        have_lsode_f = PRESENT(lsode_f)
cc
ccc     Output
cc
cc        if (my_rank == 0) then
cc          if (iout == 1) then
cc            write (*,FMT="(a)",ADVANCE="NO") 'Lagrangian 2D: '
cc          elseif (iout <= 0) then
cc            if (solenoidal) then
cc              write (*,*) 'Solenoidal lagrangian step in 2D...'
cc            else
cc              write (*,*) 'LSODE lagrangian step in 2D...'
cc            endif
cc          endif
cc        endif
cc
ccc     Assemble global temperature field
cc
cc        nxg = grid_params%nxgl(igrid)
cc        nyg = grid_params%nygl(igrid)
cc        nzg = grid_params%nzgl(igrid)
cc
cc        allocate(tg(0:nxg+1,0:nyg+1,0:nzg+1))
cc
cc        call find_global(tt,tg)
cc
ccc     Spline temperature field
cc
cc        call splineFld(tg)
cc
cc        deallocate(tg)
cc
ccc     Orbit integration
cc
cc        if (have_lsode_f) then
cc
cc          do k=1,nz
cc            do j=1,ny
cc              do i=1,nx
cc
cc                call getMGmap(i,j,k,igrid,igrid,igrid,ig,jg,kg)
cc
cc                x0 = grid_params%xx(ig)
cc                y0 = grid_params%yy(jg)
cc                z0 = grid_params%zz(kg)
cccc                call fromLocalToGlobalLimits(i,j,k,ig,jg,kg
cccc     .                                      ,igrid,igrid,igrid)
cccc
cccc                x0 = grid_params%xg(ig)
cccc                y0 = grid_params%yg(jg)
cccc                z0 = grid_params%zg(kg)
cc
cc                temp_s   = 0d0        !Initialize temperature
cc                kernel_s = 0d0        !Initialize kernel integral (for normalization)
cc
cc                call orbit_find(iout-2,x0,y0,z0,.true.,s_max=ll(i,j,k)
cc     .                         ,lsode_f=lsode_f,lsode_sum=lsode_out)
cc                temp_s   = lsode_out(1)
cc                kernel_s = lsode_out(2)
cc
cc                if ((my_rank == 0).and.(iout == 1))
cc     .                      write (*,FMT="(a)",ADVANCE="NO") '+'
cc
cc                tt(i,j,k) = temp_s/kernel_s
cc              enddo
cc            enddo
cc          enddo
cc
cc        else
cc
cc          do k=1,nz
cc            do j=1,ny
cc              do i=1,nx
cc
cc                call getMGmap(i,j,k,igrid,igrid,igrid,ig,jg,kg)
cc
cc                x0 = grid_params%xx(ig)
cc                y0 = grid_params%yy(jg)
cc                z0 = grid_params%zz(kg)
cccc                call fromLocalToGlobalLimits(i,j,k,ig,jg,kg
cccc     .                                      ,igrid,igrid,igrid)
cccc
cccc                x0 = grid_params%xg(ig)
cccc                y0 = grid_params%yg(jg)
cccc                z0 = grid_params%zg(kg)
cc
cc                temp_s   = 0d0        !Initialize temperature
cc                kernel_s = 0d0        !Initialize kernel integral (for normalization)
cc
cc                call orbit_find(iout-2,x0,y0,z0,.true.,s_max=ll(i,j,k)
cc     .                         ,line_int=int_f,lsode_sum=lsode_out)
cc
cc                if ((my_rank == 0).and.(iout == 1))
cc     .                      write (*,FMT="(a)",ADVANCE="NO") '+'
cc
cc                tt(i,j,k) = temp_s/kernel_s
cc              enddo
cc            enddo
cc          enddo
cc
cc        endif
cc
ccc     Free up orbit integral work space
cc
cc        if ((my_rank == 0).and.(iout == 1)) write (*,*)
cc
cc      end subroutine lagrangian_2d_scl

ccc     lagrangian_3d_scl
ccc     ##############################################################
cc      subroutine lagrangian_3d_scl(igrid,nx,ny,nz,tt
cc     .                            ,chi,iout,int_f,lsode_f)
cc
cc      implicit none
cc
ccc     --------------------------------------------------------------
ccc     Performs Lagrangian Green's function temperature integration
ccc     of tt and returns it in tt.
ccc     --------------------------------------------------------------
cc
ccc     Call variables
cc
cc        integer :: nx,ny,nz,igrid,iout
cc
cc        real(8) :: chi
cc
cc        real(8) :: tt(0:nx+1,0:ny+1,0:nz+1)
cc
ccc     Local variables
cc
cc        integer :: i,j,k,ig,jg,kg,nxg,nyg,nzg
cc
cc        real(8) :: x0,y0,z0
cc
cc        real(8),pointer,dimension(:,:,:) :: tg,zeros
cc
cc        logical :: have_lsode_f
cc
ccc     Interfaces
cc
cc        optional :: int_f
cc
cc        INTERFACE
cc          function int_f(fwd,xo,x,yo,y,zo,z,so,ss)
cc          logical :: int_f,fwd
cc          real(8) :: xo,x,yo,y,zo,z,so,ss
cc          end function int_f
cc        END INTERFACE
cc
cc        optional :: lsode_f
cc
cc        !Interface for optional LSODE function evaluation
cc        INTERFACE
cc          subroutine lsode_f(neq,t,y,ydot)
cc          integer :: neq
cc          real(8) :: t,y(neq),ydot(neq)
cc          end subroutine lsode_f
cc        END INTERFACE
cc
ccc     Begin program
cc
cc        if (chi==0d0) return
cc
cc        ct0 = chi*dt
cc
cc        have_lsode_f = PRESENT(lsode_f)
cc
ccc     Output
cc
cc        if (my_rank == 0) then
cc          if (iout == 1) then
cc            write (*,FMT="(a)",ADVANCE="NO") 'Lagrangian 3D: '
cc          elseif (iout <= 0) then
cc            write (*,*) 'Lagrangian step in 3D...'
cc          endif
cc        endif
cc
ccc     Assemble global temperature field
cc
cc        nxg = grid_params%nxgl(igrid)
cc        nyg = grid_params%nygl(igrid)
cc        nzg = grid_params%nzgl(igrid)
cc
cc        allocate(tg(0:nxg+1,0:nyg+1,0:nzg+1))
cc
cc        call find_global(tt,tg)
cc
ccc     Spline temperature field
cc
cc        call splineFld(tg)
cc
cc        deallocate(tg)
cc
ccc     Orbit integration
cc
cc        do k=1,nz
cc          do j=1,ny
cc            do i=1,nx
cc
cc              !Find location on global mesh
cc              call getMGmap(i,j,k,igrid,igrid,igrid,ig,jg,kg)
cc
cc              x0 = grid_params%xx(ig)
cc              y0 = grid_params%yy(jg)
cc              z0 = grid_params%zz(kg)
cccc              call fromLocalToGlobalLimits(i,j,k,ig,jg,kg
cccc     .                                    ,igrid,igrid,igrid)
cccc
cccc              x0 = grid_params%xg(ig)
cccc              y0 = grid_params%yg(jg)
cccc              z0 = grid_params%zg(kg)
cc
cc              temp_s   = 0d0    !Initialize temperature
cc              kernel_s = 0d0    !Initialize kernel integral (for normalization)
cc
cc              call orbit_find(iout-2,x0,y0,z0,.true.
cc     .                       ,line_int=int_f,lsode_f=lsode_f
cc     .                       ,lsode_sum=lsode_out)            !Forward
cc
cc              if (have_lsode_f) then
cc                temp_s   = temp_s   + lsode_out(1)
cc                kernel_s = kernel_s + lsode_out(2)
cc              endif
cc
cc              call orbit_find(iout-2,x0,y0,z0,.false.
cc     .                       ,line_int=int_f,lsode_f=lsode_f
cc     .                       ,lsode_sum=lsode_out) !Backward
cc
cc              if (have_lsode_f) then
cc                temp_s   = temp_s   + lsode_out(1)
cc                kernel_s = kernel_s + lsode_out(2)
cc              endif
cc
cc              if ((my_rank == 0).and.(iout == 1))
cc     .                    write (*,FMT="(a)",ADVANCE="NO") '+'
cc
cc              !Fill local temperature
cc              tt(i,j,k) = temp_s/kernel_s
cc            enddo
cc          enddo
cc        enddo
cc
cc        if ((my_rank == 0).and.(iout == 1)) write (*,*)
cc
cc      end subroutine lagrangian_3d_scl

c     lag_3d
c     ##############################################################
      subroutine lag_3d(igrid,nx,ny,nz,vec_in,vec_out,nints,chi,iout
     .                 ,lsode_f)

      implicit none

c     --------------------------------------------------------------
c     Performs Lagrangian Green's function temperature integration
c     of tt and returns it in tt.
c     --------------------------------------------------------------

c     Call variables

        integer :: nx,ny,nz,igrid,iout,nints

        real(8) :: chi(0:,0:,0:)

        real(8) :: vec_in(0:,0:,0:,:),vec_out(0:,0:,0:,:)

        !Interface for optional LSODE function evaluation
        INTERFACE
          subroutine lsode_f(neq,t,y,ydot)
          integer :: neq
          real(8) :: t,y(neq),ydot(neq)
          end subroutine lsode_f
        END INTERFACE

c     Local variables

        integer :: i,j,k,ig,jg,kg,nxg,nyg,nzg,ieq,neq

        real(8) :: x0,y0,z0,temp_n,gkernel_dt,temp_nm,gkernel_2dt,src_dt
     .            ,pkernel_dt,src_2dt,pkernel_2dt,b0_diag,lb0_diag

        real(8),pointer,dimension(:,:,:,:) :: vg

c     Begin program

        orbit_conf%ode_conf%lsode_neq = 5+nints !Number of LSODE integrals

c     Output

        if (my_rank == 0) then
          if (iout == 1) then
            write (*,FMT="(a)",ADVANCE="NO") 'Lagrangian 3D: '
          elseif (iout <= 0) then
            write (*,*) 'Lagrangian step in 3D...'
          endif
        endif

c     Assemble global temperature field

        nxg = grid_params%nxgl(igrid)
        nyg = grid_params%nygl(igrid)
        nzg = grid_params%nzgl(igrid)

        neq = size(vec_in,4)

        if (neq /= 3) then
          call pstop('lagrangian_3d'
     .              ,'Wrong number of input fields')
        endif

        if (size(vec_out,4) /= 5) then
          call pstop('lag_3d'
     .              ,'Wrong number of output fields')
        endif

        allocate(vg(0:nxg+1,0:nyg+1,0:nzg+1,neq)
     .          ,sp_coef (nxg+2,nyg+2,nzg+2,neq))

        call find_global(vec_in,vg)

c     Spline temperature field

        call splineFlds(vg,fldcoef=sp_coef)

        deallocate(vg)

c     Orbit integration

        b0_diag = 0d0

        do k=1,nz
          do j=1,ny
            do i=1,nx

              if (chi(i,j,k) == 0d0) cycle

              call getMGmap(i,j,k,igrid,igrid,igrid,ig,jg,kg)

              x0 = grid_params%xx(ig)
              y0 = grid_params%yy(jg)
              z0 = grid_params%zz(kg)

              ct0 = chi(i,j,k)*dt
              ct1 = chi(i,j,k)*(dt+dtold)

              call orbit_find(iout-2,x0,y0,z0,.true.,lsode_f=lsode_f
     .                       ,lsode_sum=lsode_out) !Forward

              temp_n      = lsode_out(1)
              gkernel_dt  = lsode_out(2)
              temp_nm     = lsode_out(3)
              gkernel_2dt = lsode_out(4)
              src_dt      = lsode_out(5)
              pkernel_dt  = lsode_out(6)
              src_2dt     = lsode_out(7)
              pkernel_2dt = lsode_out(8)

              call orbit_find(iout-2,x0,y0,z0,.false.,lsode_f=lsode_f
     .                       ,lsode_sum=lsode_out) !Backward

              temp_n      = temp_n      + lsode_out(1)
              gkernel_dt  = gkernel_dt  + lsode_out(2)
              temp_nm     = temp_nm     + lsode_out(3)
              gkernel_2dt = gkernel_2dt + lsode_out(4)
              src_dt      = src_dt      + lsode_out(5)
              pkernel_dt  = pkernel_dt  + lsode_out(6)
              src_2dt     = src_2dt     + lsode_out(7)
              pkernel_2dt = pkernel_2dt + lsode_out(8)

              if ((my_rank == 0).and.(iout == 1))
     .                    write (*,FMT="(a)",ADVANCE="NO") '+'

              vec_out(i,j,k,1) = temp_n !/gkernel_dt
              vec_out(i,j,k,2) = temp_nm!/gkernel_2dt
              vec_out(i,j,k,3) = src_dt !/pkernel_dt 
              vec_out(i,j,k,4) = src_2dt!/pkernel_2dt
              vec_out(i,j,k,5) = pkernel_dt
            enddo
          enddo
        enddo

        if ((my_rank == 0).and.(iout == 1)) write (*,*)

c     Free up orbit integral work space

        deallocate(sp_coef)

      end subroutine lag_3d

c     lag_3d_new
c     ##############################################################
      subroutine lag_3d_new(igrid,nx,ny,nz,vec_in,vec_out,chi,iout)

      implicit none

c     --------------------------------------------------------------
c     Performs Lagrangian Green's function temperature integration
c     of tt and returns it in tt.
c     --------------------------------------------------------------

c     Call variables

        integer :: nx,ny,nz,igrid,iout

        real(8) :: chi(0:,0:,0:)

        real(8) :: vec_in(0:,0:,0:,:),vec_out(0:,0:,0:,:)

        !Interface for optional LSODE function evaluation
        INTERFACE
          subroutine lsode_3d_bdf2_all(neq,t,y,ydot)
          integer :: neq
          real(8) :: t,y(neq),ydot(neq)
          end subroutine lsode_3d_bdf2_all
        END INTERFACE

c     Local variables

        integer :: i,j,k,ig,jg,kg,nxg,nyg,nzg,ieq,neq,nints

        real(8) :: x0,y0,z0,temp_n,gkernel_dt,temp_nm,gkernel_2dt,src_dt
     .            ,pkernel_dt,src_2dt,pkernel_2dt,b0_diag,lb0_diag

        real(8),pointer,dimension(:,:,:,:) :: vg

c     Begin program

        nints = size(vec_out,4)
        orbit_conf%ode_conf%lsode_neq = 5+nints !Number of LSODE integrals

        if (nints > LSODE_ACC_DIM) then
          call pstop('lag_3d_new','Too many LSODE integrals')
        endif

c     Output

        if (my_rank == 0) then
          if (iout == 1) then
            write (*,FMT="(a)",ADVANCE="NO") 'Lagrangian 3D: '
          elseif (iout <= 0) then
            write (*,*) 'Lagrangian step in 3D...'
          endif
        endif

c     Assemble global temperature field

        nxg = grid_params%nxgl(igrid)
        nyg = grid_params%nygl(igrid)
        nzg = grid_params%nzgl(igrid)

        neq = size(vec_in,4)

cc        if (neq /= 3) then
cc          call pstop('lagrangian_3d'
cc     .              ,'Wrong number of input fields')
cc        endif
cc
cc        if (size(vec_out,4) /= 5) then
cc          call pstop('lag_3d'
cc     .              ,'Wrong number of output fields')
cc        endif

        allocate(vg(0:nxg+1,0:nyg+1,0:nzg+1,neq)
     .          ,sp_coef (nxg+2,nyg+2,nzg+2,neq))

        call find_global(vec_in,vg)

c     Spline temperature field

        call splineFlds(vg,fldcoef=sp_coef)

        deallocate(vg)

c     Orbit integration

        b0_diag = 0d0

        do k=1,nz
          do j=1,ny
            do i=1,nx

              if (chi(i,j,k) == 0d0) cycle

              call getMGmap(i,j,k,igrid,igrid,igrid,ig,jg,kg)

              x0 = grid_params%xx(ig)
              y0 = grid_params%yy(jg)
              z0 = grid_params%zz(kg)

              ct0 = chi(i,j,k)*dt
              ct1 = chi(i,j,k)*(dt+dtold)

              call orbit_find(iout-2,x0,y0,z0,.true.
     .                       ,lsode_f=lsode_3d_bdf2_all
     .                       ,lsode_sum=lsode_out) !Forward

              vec_out(i,j,k,:) = lsode_out(1:nints)

              call orbit_find(iout-2,x0,y0,z0,.false.
     .                       ,lsode_f=lsode_3d_bdf2_all
     .                       ,lsode_sum=lsode_out) !Backward

              vec_out(i,j,k,:) = vec_out(i,j,k,:) + lsode_out

              if ((my_rank == 0).and.(iout == 1))
     .                    write (*,FMT="(a)",ADVANCE="NO") '+'
            enddo
          enddo
        enddo

        if ((my_rank == 0).and.(iout == 1)) write (*,*)

c     Free up orbit integral work space

        deallocate(sp_coef)

      end subroutine lag_3d_new

c     lag_2d
c     ##############################################################
      subroutine lag_2d(igrid,nx,ny,nz,vec_in,vec_out,nints,chi,iout
     .                 ,lsode_f)

      implicit none

c     --------------------------------------------------------------
c     Performs Lagrangian integration of vector field vec_in according 
c     to B field (bx,by,bz). Result is stored on vector field vec_out:
c       * vec_in: number of fields to integrate (spline)
c       * vec_out: number of integrals to be performed on vec_in fields.
c     Thus, vec_in and vec_out can have different 4th dimension.
c     --------------------------------------------------------------

c     Call variables

        integer :: nx,ny,nz,igrid,iout,nints

        real(8) :: chi(0:,0:,0:)

        real(8) :: vec_in(0:,0:,0:,:),vec_out(0:,0:,0:,:)

        !Interface for optional LSODE function evaluation
        INTERFACE
          subroutine lsode_f(neq,t,y,ydot)
          integer :: neq
          real(8) :: t,y(neq),ydot(neq)
          end subroutine lsode_f
        END INTERFACE

c     Local variables

        integer :: i,j,k,ig,jg,kg,nxg,nyg,nzg,ieq,neq

        real(8) :: x0,y0,z0,temp_n,gkernel_dt,temp_nm,gkernel_2dt,src_dt
     .            ,pkernel_dt,src_2dt,pkernel_2dt

        real(8),pointer,dimension(:,:,:,:) :: vg

c     Begin program

        orbit_conf%ode_conf%lsode_neq = 3+nints

c     Output

        if (my_rank == 0) then
          if (iout == 1) then
            write (*,FMT="(a)",ADVANCE="NO") 'Lagrangian 2D: '
          elseif (iout <= 0) then
            write (*,*) 'LSODE lagrangian step in 2D...'
          endif
        endif

c     Assemble global temperature field

        nxg = grid_params%nxgl(igrid)
        nyg = grid_params%nygl(igrid)
        nzg = grid_params%nzgl(igrid)

        neq = size(vec_in,4)

        if (neq /= 3) then
          call pstop('lag_2d'
     .              ,'Wrong number of input fields')
        endif

        if (size(vec_out,4) /= 5) then
          call pstop('lag_2d'
     .              ,'Wrong number of output fields')
        endif

        allocate(vg(0:nxg+1,0:nyg+1,0:nzg+1,neq)
     .          ,sp_coef (nxg+2,nyg+2,nzg+2,neq))

        call find_global(vec_in,vg)

c     Spline fields

        call splineFlds(vg,fldcoef=sp_coef)

        deallocate(vg)

c     LSODE orbit integration

        do k=1,nz
          do j=1,ny
            do i=1,nx

              if (chi(i,j,k) == 0d0) cycle

              call getMGmap(i,j,k,igrid,igrid,igrid,ig,jg,kg)

              x0 = grid_params%xx(ig)
              y0 = grid_params%yy(jg)
              z0 = grid_params%zz(kg)

              ct0 = chi(i,j,k)*dt
              ct1 = chi(i,j,k)*(dt+dtold)

              call orbit_find(iout-2,x0,y0,z0,.true.
     .                       ,s_max=ll(i,j,k)
     .                       ,lsode_f=lsode_f,lsode_sum=lsode_out)

              if ((my_rank == 0).and.(iout == 1))
     .                    write (*,FMT="(a)",ADVANCE="NO") '+'

              temp_n      = lsode_out(1)
              temp_nm     = lsode_out(3)
              src_dt      = lsode_out(5)
              src_2dt     = lsode_out(7)

              gkernel_dt  = lsode_out(2)
              gkernel_2dt = lsode_out(4)
              pkernel_dt  = lsode_out(6)
              pkernel_2dt = lsode_out(8)

              !Parallel Xport piece
              vec_out(i,j,k,1) = temp_n !/gkernel_dt
              vec_out(i,j,k,2) = temp_nm!/gkernel_2dt

              !Source
              vec_out(i,j,k,3) = src_dt !/pkernel_dt 
              vec_out(i,j,k,4) = src_2dt!/pkernel_2dt

              !beta integral
              vec_out(i,j,k,5) = pkernel_dt
cc              if (ct0 < CHI_LIMIT) then !Delta function
cc                vec_out(i,j,k,3) = vec_in(i,j,k,3)
cc                vec_out(i,j,k,4) = vec_in(i,j,k,3)
cc              else
cc                vec_out(i,j,k,3) = src_dt !/pkernel_dt 
cc                vec_out(i,j,k,4) = src_2dt!/pkernel_2dt
cc              endif

            enddo
          enddo
        enddo

        if ((my_rank == 0).and.(iout == 1)) write (*,*)

c     Free up orbit integral work space

        deallocate(sp_coef)

      end subroutine lag_2d

c     lag_2d_new
c     ##############################################################
      subroutine lag_2d_new(igrid,nx,ny,nz,vec_in,vec_out,chi,iout)

      implicit none

c     --------------------------------------------------------------
c     Performs Lagrangian integration of vector field vec_in according 
c     to B field (bx,by,bz). Result is stored on vector field vec_out:
c       * vec_in: number of fields to integrate (spline)
c       * vec_out: number of integrals to be performed on vec_in fields.
c     Thus, vec_in and vec_out can have different 4th dimension.
c     --------------------------------------------------------------

c     Call variables

        integer :: nx,ny,nz,igrid,iout

        real(8) :: chi(0:,0:,0:)

        real(8) :: vec_in(0:,0:,0:,:),vec_out(0:,0:,0:,:)

        !Interface for optional LSODE function evaluation
        INTERFACE
          subroutine lsode_2d_bdf2_3eq_all(neq,t,y,ydot)
          integer :: neq
          real(8) :: t,y(neq),ydot(neq)
          end subroutine lsode_2d_bdf2_3eq_all
        END INTERFACE

c     Local variables

        integer :: i,j,k,ig,jg,kg,nxg,nyg,nzg,ieq,neq,nints

        real(8) :: x0,y0,z0,temp_n,gkernel_dt,temp_nm,gkernel_2dt,src_dt
     .            ,pkernel_dt,src_2dt,pkernel_2dt

        real(8),pointer,dimension(:,:,:,:) :: vg

c     Begin program

        nints = size(vec_out,4) !Number of LSODE integrals

        orbit_conf%ode_conf%lsode_neq = 3+nints

        if (nints > LSODE_ACC_DIM) then
          call pstop('lag_2d_new','Too many LSODE integrals')
        endif

c     Output

        if (my_rank == 0) then
          if (iout == 1) then
            write (*,FMT="(a)",ADVANCE="NO") 'Lagrangian 2D: '
          elseif (iout <= 0) then
            write (*,*) 'LSODE lagrangian step in 2D...'
          endif
        endif

c     Assemble global temperature field

        nxg = grid_params%nxgl(igrid)
        nyg = grid_params%nygl(igrid)
        nzg = grid_params%nzgl(igrid)

        neq = size(vec_in,4)

cc        if (neq /= 3) then
cc          call pstop('lag_2d_new'
cc     .              ,'Wrong number of input fields')
cc        endif
cc
cc        if (size(vec_out,4) /= 5) then
cc          call pstop('lag_2d'
cc     .              ,'Wrong number of output fields')
cc        endif

        allocate(vg(0:nxg+1,0:nyg+1,0:nzg+1,neq)
     .          ,sp_coef (nxg+2,nyg+2,nzg+2,neq))

        call find_global(vec_in,vg)

c     Spline fields

        call splineFlds(vg,fldcoef=sp_coef)

        deallocate(vg)

c     LSODE orbit integration

        do k=1,nz
          do j=1,ny
            do i=1,nx

              if (chi(i,j,k) == 0d0) cycle

              call getMGmap(i,j,k,igrid,igrid,igrid,ig,jg,kg)

              x0 = grid_params%xx(ig)
              y0 = grid_params%yy(jg)
              z0 = grid_params%zz(kg)

              ct0 = chi(i,j,k)*dt
              ct1 = chi(i,j,k)*(dt+dtold)

              call orbit_find(iout-2,x0,y0,z0,.true.
     .                       ,s_max=ll(i,j,k)
     .                       ,lsode_f=lsode_2d_bdf2_3eq_all
     .                       ,lsode_sum=lsode_out)

              if ((my_rank == 0).and.(iout == 1))
     .                    write (*,FMT="(a)",ADVANCE="NO") '+'

              vec_out(i,j,k,:) = lsode_out(1:nints)
            enddo
          enddo
        enddo

        if ((my_rank == 0).and.(iout == 1)) write (*,*)

c     Free up orbit integral work space

        deallocate(sp_coef)

      end subroutine lag_2d_new

c     lag_2d_init
c     ##############################################################
      subroutine lag_2d_init(igrid,nx,ny,nz,bb,alp,iout
     .                      ,G_tol,spline_order,B_wgt
     .                      ,psi,chi)

      implicit none

c     --------------------------------------------------------------
c     Initializes 2D Lagrangian step.
c     --------------------------------------------------------------

c     Call variables

        integer :: nx,ny,nz,igrid,iout,spline_order

        real(8),INTENT(IN) :: bb (0:nx+1,0:ny+1,0:nz+1,3)

        real(8) :: G_tol,alp

        logical :: B_wgt

        real(8),INTENT(IN),optional :: psi(0:nx+1,0:ny+1,0:nz+1)
     .                                ,chi(0:nx+1,0:ny+1,0:nz+1)

c     Local variables

        integer :: i,j,k,ig,jg,kg,nxg,nyg,nzg

        real(8) :: x0,y0,z0,beta_err

        real(8),pointer,dimension(:,:,:  ) :: jacg,psig,psil
        real(8),pointer,dimension(:,:,:,:) :: bcar,bg,bcarg,carg

        logical :: have_chi

c     Interfaces

        INTERFACE
          function l_int(fwd,xo,x,yo,y,zo,z,so,ss)
          logical :: l_int,fwd
          real(8) :: xo,x,yo,y,zo,z,so,ss
          end function l_int
        END INTERFACE

c     Begin program

        init_lag = .true.

        l2d = .true.   !Indicate 2D treatment

        solenoidal = PRESENT(psi)

        have_chi = PRESENT(chi)

        alpha = alp

        nxg = grid_params%nxgl(igrid)
        nyg = grid_params%nygl(igrid)
        nzg = grid_params%nzgl(igrid)

        bwgt = B_wgt

c     Compute Green's integral cut-offs

        t_tol = G_tol !Initialize Green's function tolerance

        T_cutoff2 = 1d0
        S_cutoff2 = 1d0
        do i=1,10
          T_cutoff2 = abs(log(  sqrt(pi)*t_tol*sqrt(T_cutoff2)))
          S_cutoff2 = abs(log(2*sqrt(pi)*t_tol*     S_cutoff2 ))
        enddo

c     Find local cartesian components

        allocate(bcar(0:nx+1,0:ny+1,0:nz+1,3))

        bcar = bb

        call transformVector(igrid,igrid,igrid
     .                      ,0,nx+1,0,ny+1,0,nz+1
     .                      ,bcar(:,:,:,1)
     .                      ,bcar(:,:,:,2)
     .                      ,bcar(:,:,:,3),'cnv','car')

c     Fill global arrays

        allocate(bg   (0:nxg+1,0:nyg+1,0:nzg+1,3)
     .          ,bcarg(0:nxg+1,0:nyg+1,0:nzg+1,3)
     .          ,carg (0:nxg+1,0:nyg+1,0:nzg+1,3)
     .          ,jacg (0:nxg+1,0:nyg+1,0:nzg+1))

        call find_global(bb,bg)

        call find_global(bcar,bcarg)

        call find_global(gmetric%grid(igrid)%jac,jacg)
        call find_global(gmetric%grid(igrid)%car,carg)

c     Output

        if (my_rank == 0) then
          write (*,*)
          write (*,*) '2D Lagrangian initialization...'
        endif

c     Computation of periodic orbit quantities (length, beta)

        allocate(ll(0:nx+1,0:ny+1,0:nz+1)
     .          ,b0(0:nx+1,0:ny+1,0:nz+1))

        !Length
        call find_length

        !Beta iteration
        call find_beta

        if ((my_rank == 0).and.(iout == 1)) then
          write (*,*)
          write (*,*)
        endif

c diag ****
cc        if (np == 1) then
cc          write (*,*) 'DIAG -- lag_2d_init, dumping plots'
cc          open(unit=110,file='per-orbits-qtys.bin',form='unformatted'
cc     .        ,status='replace')
cc          call contour(ll(1:nx,1:ny,1),nx,ny,0d0,xmax,0d0,ymax,0,110)
cc          call contour(log10(b0(1:nx,1:ny,1)),nx,ny,0d0,xmax,0d0,ymax,1
cc     .                ,110)
cc          b0 = b0*sum(bcar(:,:,:,:)**2,dim=4)
cc          call contour(b0(1:nx,1:ny,1),nx,ny,0d0,xmax,0d0,ymax,1,110)
cc          if (have_chi) then
cc            b0 = chi*dt/b0*sum(bcar(:,:,:,:)**2,dim=4)
cc          else
cc            b0 = sum(bcar(:,:,:,:)**2,dim=4)
cc          endif
cc          call contour(b0(1:nx,1:ny,1),nx,ny,0d0,xmax,0d0,ymax,1,110)
cc          close(110)
cc          stop
cc        endif
c diag ****

c     Setup periodic integrator

        if (solenoidal) then

          allocate(psil(0:nx +1,0:ny +1,0:nz +1)
     .            ,psig(0:nxg+1,0:nyg+1,0:nzg+1))

          psil = psi
          call find_global(psil,psig)

          bg(:,:,:,1:2) = 0d0
          bg(:,:,:,3  ) = -psig

          call orbit_setup(nxg+2,nyg+2,nzg+2
     .                    ,grid_params%xg
     .                    ,grid_params%yg
     .                    ,grid_params%zg
     .                    ,bg(:,:,:,1)
     .                    ,bg(:,:,:,2)
     .                    ,bg(:,:,:,3)
     .                    ,bcarg(:,:,:,1)
     .                    ,bcarg(:,:,:,2)
     .                    ,bcarg(:,:,:,3)
     .                    ,jacg
     .                    ,carg
     .                    ,bcond
     .                    ,dtime=orbit_dt
     .                    ,atol=orbit_atol  !Field line integration absolute tolerance
     .                    ,rtol=orbit_rtol  !Field line integration relative tolerance
     .                    ,spline_order=spline_order
     .                    ,solen=.true.
     .                    ,B_input_is_A=.true.
     .                    )

          deallocate(psig,psil)

        else

          call orbit_setup(nxg+2,nyg+2,nzg+2
     .                    ,grid_params%xg
     .                    ,grid_params%yg
     .                    ,grid_params%zg
     .                    ,bg(:,:,:,1)
     .                    ,bg(:,:,:,2)
     .                    ,bg(:,:,:,3)
     .                    ,bcarg(:,:,:,1)
     .                    ,bcarg(:,:,:,2)
     .                    ,bcarg(:,:,:,3)
     .                    ,jacg
     .                    ,carg
     .                    ,bcond
     .                    ,dtime=orbit_dt
     .                    ,atol=orbit_atol      !Field line integration absolute tolerance
     .                    ,rtol=orbit_rtol      !Field line integration relative tolerance
     .                    ,spline_order=spline_order
     .                    ,solen=.false.
     .                    ,ds_is_Bds=bwgt       !Use B.ds as integration variable
cc     .                    ,lsode_neq = 3+nints  !3 fields + number of field integrals
     .                    ,lsode_to_smax=.true. !Use arclength as independent variable
     .                    )

        endif

        deallocate(bg,bcarg,jacg,carg,bcar)

c     End program

      contains

c     find_length
c     ##############################################################
      subroutine find_length
c     --------------------------------------------------------------
c     Finds periodic orbit length and stores in ll array
c     --------------------------------------------------------------

c     Setup orbit integrator

        if (solenoidal) then

          if (bwgt) then
            call pstop('lag_2d_init'
     .                ,'Solen version cannot do B.ds yet')
          endif

        else

          call orbit_setup(nxg+2,nyg+2,nzg+2
     .                    ,grid_params%xg
     .                    ,grid_params%yg
     .                    ,grid_params%zg
     .                    ,bg(:,:,:,1)
     .                    ,bg(:,:,:,2)
     .                    ,bg(:,:,:,3)
     .                    ,bcarg(:,:,:,1)
     .                    ,bcarg(:,:,:,2)
     .                    ,bcarg(:,:,:,3)
     .                    ,jacg
     .                    ,carg
     .                    ,bcond
     .                    ,dtime=orbit_dt
     .                    ,atol=orbit_atol       !Field line integration absolute tolerance
     .                    ,rtol=orbit_rtol       !Field line integration relative tolerance
     .                    ,spline_order=spline_order
     .                    ,solen=.false.
     .                    ,ds_is_Bds=bwgt        !Use B.ds as integration variable
     .                    ,lsode_dt_lim=.true.   !Limit LSODE time step by orbit_dt
     .                    ,lsode_neq = 3         !Number of equations (only integrate field line)
     .                    ,lsode_s_indep=.true.  !Use arclength as independent variable
     .                    )

        endif

c     Integrate length

        ll = 0d0
        if (my_rank == 0) write (*,*) '--> LSODE 2D length integral...'

        do k=1,nz
          do j=1,ny
            do i=1,nx

              call getMGmap(i,j,k,igrid,igrid,igrid,ig,jg,kg)

              x0 = grid_params%xx(ig)
              y0 = grid_params%yy(jg)
              z0 = grid_params%zz(kg)

              length_s = 0d0    !Initialize orbit length
              b0_s     = 0d0

              call orbit_find(iout-2,x0,y0,z0,.true.,line_int=l_int)

              if ((my_rank == 0).and.(iout == 1))
     .                    write (*,FMT="(a)",ADVANCE="NO") '+'

              ll(i,j,k) = length_s       !Find int(dlambda)=int(B.ds)
              b0(i,j,k) = b0_s/length_s  !Find int(dlambda/B^2)/int(dlambda)=1/<B^2>
            enddo
          enddo
        enddo

        call orbit_destroy

      end subroutine find_length

c     find_beta
c     ##############################################################
      subroutine find_beta

c     --------------------------------------------------------------
c     Finds periodic orbit beta and stores in b0 array
c     --------------------------------------------------------------

#if defined(FPA)
        use fixed_point_accelerator
#endif

        real(8) :: beta_err,lberr

#if defined(FPA)
        type (fpa_state) :: fpav
#endif
        INTERFACE
          subroutine lsode_2d_beta(neq,t,y,ydot)
          integer :: neq
          real(8) :: t,y(neq),ydot(neq)
          end subroutine lsode_2d_beta
        END INTERFACE

        real(8) :: x(nx*ny*nz),r(nx*ny*nz)

c     Begin program

        if (.not.b_wgt) then
          b0 = 1d0              !Set beta to unity 
          return
        elseif (.not.have_chi) then
          return
        endif

        if (my_rank == 0) write (*,*) '--> LSODE 2D beta integral...'

c     Setup orbits

        call orbit_setup(nxg+2,nyg+2,nzg+2
     .                  ,grid_params%xg
     .                  ,grid_params%yg
     .                  ,grid_params%zg
     .                  ,bg(:,:,:,1)
     .                  ,bg(:,:,:,2)
     .                  ,bg(:,:,:,3)
     .                  ,bcarg(:,:,:,1)
     .                  ,bcarg(:,:,:,2)
     .                  ,bcarg(:,:,:,3)
     .                  ,jacg
     .                  ,carg
     .                  ,bcond
     .                  ,dtime=orbit_dt
     .                  ,atol=orbit_atol       !Field line integration absolute tolerance
     .                  ,rtol=orbit_rtol       !Field line integration relative tolerance
     .                  ,spline_order=spline_order
     .                  ,solen=.false.
     .                  ,ds_is_Bds=bwgt        !Use B.ds as integration variable
     .                  ,lsode_neq = 4         !Number of equations
     .                  ,lsode_to_smax=.true.  !Integrate to smax (implies lsode_s_indep)
     .                  )

c      Compute beta

        !Initial guess (b0 = 1/<B^2> by default)
cc        where (chi*dt*sum(bcar(:,:,:,:)**2,dim=4) < CHI_LIMIT)
cc        where (chi*dt/b0 < CHI_LIMIT)
          b0 = 1d0/sum(bcar(:,:,:,:)**2,dim=4) !Set beta_0=1/B^2 for small B
cc        end where

cc          b0 = 1d0 !Normalization test

        !Iteration
#if defined(FPA)
        !Create fixed-point accelerator
        call fpa_create(fpav,x,maxv=5)

        if (my_rank == 0) then
          write (*,*) '  Using FPA acceleration...'
        endif
#endif

        beta_err = 1d0

        !Picard iteration
        x = reshape(b0(1:nx,1:ny,1:nz),shape(x))

        do while (beta_err > 1d-2)

          beta_err = 0d0

          do k=1,nz
            do j=1,ny
              do i=1,nx

                ct0 = chi(i,j,k)*dt/b0(i,j,k)

cc                if (ct0 < CHI_LIMIT) then  !Delta function
cc                  b0_s = 1d0/sum(bcar(i,j,k,:)**2)
cc                  cycle
cc                else
                  call getMGmap(i,j,k,igrid,igrid,igrid,ig,jg,kg)

                  x0 = grid_params%xx(ig)
                  y0 = grid_params%yy(jg)
                  z0 = grid_params%zz(kg)

                  call orbit_find(iout-2,x0,y0,z0,.true.
     .                           ,s_max=ll(i,j,k)
     .                           ,lsode_f=lsode_2d_beta
     .                           ,lsode_sum=lsode_out)

                  b0_s = lsode_out(1)  !Store LSODE accumulator

                  if ((my_rank == 0).and.(iout == 1))
     .                        write (*,FMT="(a)",ADVANCE="NO") '+'
cc                endif

                b0(i,j,k) = b0_s

              enddo
            enddo
          enddo

          !Compute residual vector
          r = x - reshape(b0(1:nx,1:ny,1:nz),shape(x))

          !Convergence test
          beta_err = sqrt(dot2(r/x,r/x)/nxg/nyg/nzg)

          if (my_rank == 0) write (*,*) '  >> Beta error=',beta_err

#if defined(FPA)
          !Apply fixed-point accelerator
          call fpa_correction(fpav,r,dp=dot2)

          !Update solution
          x = x - r

          b0(1:nx,1:ny,1:nz) = reshape(x,shape(b0(1:nx,1:ny,1:nz)))
#else
          x = reshape(b0(1:nx,1:ny,1:nz),shape(x))
#endif
        enddo

c     Destroy orbit

        call orbit_destroy

#if defined(FPA)
        !Destroy fixed-point accelerator
        call fpa_destroy(fpav)
#endif

      end subroutine find_beta

      end subroutine lag_2d_init

c     lag_3d_init
c     ##############################################################
      subroutine lag_3d_init(igrid,nx,ny,nz,bb,alp,iout,G_tol
     .                      ,spline_order,B_wgt,chi,poinc_file)
      use iosetup

      implicit none

c     --------------------------------------------------------------
c     Initializes 3D Lagrangian step.
c     --------------------------------------------------------------

c     Call variables

        integer :: nx,ny,nz,igrid,spline_order,iout
 
        real(8) :: G_tol,alp
        real(8),dimension(0:,0:,0:,:) :: bb

        logical :: B_wgt

        real(8),INTENT(IN),optional :: chi(0:nx+1,0:ny+1,0:nz+1)

        character(*),optional :: poinc_file

c     Local variables

        integer :: nxg,nyg,nzg,i,j,k,ig,jg,kg
        real(8) :: x0,y0,z0

        real(8),pointer,dimension(:)       :: xx,yy,zz
        real(8),pointer,dimension(:,:,:  ) :: jacg
        real(8),pointer,dimension(:,:,:,:) :: bcar,bg,bcarg,carg

        logical :: have_chi

c     Begin program

        if (my_rank == 0) then
          write (*,*)
          write (*,*) '3D Lagrangian initialization...'
        endif

        init_lag = .true.

        l2d = .false.   !Indicate 3D treatment

        have_chi = PRESENT(chi)

        alpha = alp

        nxg = grid_params%nxgl(igrid)
        nyg = grid_params%nygl(igrid)
        nzg = grid_params%nzgl(igrid)

        bwgt = B_wgt

c     Compute Green's integral cut-offs

        t_tol = G_tol !Initialize Green's function tolerance

        T_cutoff2 = 1d0
        S_cutoff2 = 1d0
        do i=1,10
          T_cutoff2 = abs(log(  sqrt(pi)*t_tol*sqrt(T_cutoff2)))
          S_cutoff2 = abs(log(2*sqrt(pi)*t_tol*     S_cutoff2 ))
        enddo

c     Find global arrays

        if (PRESENT(poinc_file)) then
          call get_global_data_from_file
        else
          call get_global_data
        endif

c     Find initial beta=1/B^2

        allocate(b0(0:nx+1,0:ny+1,0:nz+1))

        call find_beta_3d

c diag ****
cc        if (np == 1.and.b_wgt) then
cc          write (*,*) 'DIAG -- lag_3d_init, dumping plots'
cc          open(unit=110,file='per-orbits-qtys.bin',form='unformatted'
cc     .        ,status='replace')
cc          call contour(log10(b0(1:nx,1:ny,1)),nx,ny,0d0,xmax,0d0,ymax,0
cc     .                ,110)
cc          b0 = b0*sum(bcar(:,:,:,:)**2,dim=4)
cc          call contour(b0(1:nx,1:ny,1),nx,ny,0d0,xmax,0d0,ymax,1,110)
cc          if (have_chi) then
cc            b0 = chi*dt/b0*sum(bcar(:,:,:,:)**2,dim=4)
cc          else
cc            b0 = sum(bcar(:,:,:,:)**2,dim=4)
cc          endif
cc          call contour(b0(1:nx,1:ny,1),nx,ny,0d0,xmax,0d0,ymax,1,110)
cc          close(110)
cccc          stop
cc        endif
c diag ****

c     Orbit initialization on global array

        call orbit_setup(nxg+2,nyg+2,nzg+2
     .                  ,xx
     .                  ,yy
     .                  ,zz
     .                  ,bg(:,:,:,1)
     .                  ,bg(:,:,:,2)
     .                  ,bg(:,:,:,3)
     .                  ,bcarg(:,:,:,1)
     .                  ,bcarg(:,:,:,2)
     .                  ,bcarg(:,:,:,3)
     .                  ,jacg
     .                  ,carg
     .                  ,bcond
     .                  ,dtime=orbit_dt
     .                  ,atol=orbit_atol  !Field line integration absolute tolerance
     .                  ,rtol=orbit_rtol  !Field line integration relative tolerance
     .                  ,spline_order=spline_order
     .                  ,solen=.false.
cc     .                  ,lsode_neq=5 + nints  !5 field eqs + number of field integrals
     .                  ,ds_is_Bds=bwgt   !Use B.ds as integration variable
     .                  )

c     Deallocate memory

        deallocate(bg,bcarg,jacg,carg,bcar,xx,yy,zz)

c     End program

      contains

c     get_global_data
c     ##############################################################
      subroutine get_global_data

c     --------------------------------------------------------------
c     Finds global data (B,geometry) from Poincare file
c     --------------------------------------------------------------

c     Find local cartesian components

        allocate(bcar(0:nx+1,0:ny+1,0:nz+1,3))

        bcar = bb

        call transformVector(igrid,igrid,igrid
     .                      ,0,nx+1,0,ny+1,0,nz+1
     .                      ,bcar(:,:,:,1)
     .                      ,bcar(:,:,:,2)
     .                      ,bcar(:,:,:,3),'cnv','car')

c     Fill global arrays

        allocate(xx(nxg+2),yy(nyg+2),zz(nzg+2))

        allocate(bg   (0:nxg+1,0:nyg+1,0:nzg+1,3)
     .          ,bcarg(0:nxg+1,0:nyg+1,0:nzg+1,3)
     .          ,carg (0:nxg+1,0:nyg+1,0:nzg+1,3)
     .          ,jacg (0:nxg+1,0:nyg+1,0:nzg+1))

        call find_global(bb  ,bg)
        call find_global(bcar,bcarg)

        call find_global(gmetric%grid(igrid)%jac,jacg)
        call find_global(gmetric%grid(igrid)%car,carg)

        xx = grid_params%xg
        yy = grid_params%yg
        zz = grid_params%zg

      end subroutine get_global_data

c     get_global_data_from_file
c     ##############################################################
      subroutine get_global_data_from_file

c     --------------------------------------------------------------
c     Finds global data (B,geometry) from Poincare file
c     --------------------------------------------------------------

        implicit none

        integer :: upixie=12345,ierr

        integer,pointer,dimension(:) :: mupixie

        integer :: ifl,ilmg,ihpg,jlmg,jhpg,klmg,khpg,nx,ny,nz,nfiles
     .            ,bcnd(6)

        logical :: tor

c     Begin program

        if (my_rank == 0) then
          write (*,*) 'Asynchronous input in Lagrangian step'
        endif

c     Check file configuration

        ierr = openFiles(my_rank,np,upixie,trim(poinc_file),mupixie)

        if (ierr /= 0) call pstop('get_global_data_from_file'
     .                           ,'Cannot find files')

c     Read setup info

        nfiles = size(mupixie)

        do ifl=1,nfiles
          read (mupixie(ifl)) nxg,nyg,nzg,bcnd,tor
        enddo

        if (tor .and. (coords /= 'vmc' .and. coords /= 'tor')) then
          call closeFiles(mupixie)
          call pstop('lag_init','Coordinates are not toroidal')
        endif

        if (sum(bcnd-bcond) /= 0) then
          call closeFiles(mupixie)
          call pstop('lag_init','BC conf. does not agree')
        endif

        if (nxd /= nxg .or. nyd /= nyg .or. nzd /= nzg) then
          call closeFiles(mupixie)
          write (*,*) 'Grid size=',nxg,nyg,nzg
          call pstop('lag_init'
     .              ,'Global grid size does not agree')
        endif

c     Allocate global arrays

        allocate(xx(nxg+2),yy(nyg+2),zz(nzg+2))

        allocate(bg   (0:nxg+1,0:nyg+1,0:nzg+1,3)
     .          ,bcarg(0:nxg+1,0:nyg+1,0:nzg+1,3)
     .          ,carg (0:nxg+1,0:nyg+1,0:nzg+1,3)
     .          ,jacg (0:nxg+1,0:nyg+1,0:nzg+1  ))

c     Read global arrays

        do ifl=1,nfiles

          read (mupixie(ifl)) xx
          read (mupixie(ifl)) yy
          read (mupixie(ifl)) zz

          if (nfiles > 1) then
            read (mupixie(ifl)) ilmg,ihpg,jlmg,jhpg,klmg,khpg
          else
            ilmg = 0
            ihpg = nxg+1

            jlmg = 0
            jhpg = nyg+1

            klmg = 0
            khpg = nzg+1
          endif

          read (mupixie(ifl)) bg   (ilmg:ihpg,jlmg:jhpg,klmg:khpg,1)
          read (mupixie(ifl)) bg   (ilmg:ihpg,jlmg:jhpg,klmg:khpg,2)
          read (mupixie(ifl)) bg   (ilmg:ihpg,jlmg:jhpg,klmg:khpg,3)
          read (mupixie(ifl)) bcarg(ilmg:ihpg,jlmg:jhpg,klmg:khpg,1)
          read (mupixie(ifl)) bcarg(ilmg:ihpg,jlmg:jhpg,klmg:khpg,2)
          read (mupixie(ifl)) bcarg(ilmg:ihpg,jlmg:jhpg,klmg:khpg,3)
          read (mupixie(ifl)) jacg (ilmg:ihpg,jlmg:jhpg,klmg:khpg)
          read (mupixie(ifl)) carg (ilmg:ihpg,jlmg:jhpg,klmg:khpg,3)

        enddo

c     Assign local B arrays

        bb = bg(ilomg:ihipg,jlomg:jhipg,klomg:khipg,:)

c     Update grid metrics

        nxl = ihipg-ilomg-1
        nyl = jhipg-jlomg-1
        nzl = khipg-klomg-1

        call defineGridMetric(grid_params
     .           ,xcar=carg(ilomg:ihipg,jlomg:jhipg,klomg:khipg,:))

        gmetric => grid_params%gmetric

c     Close files

        call closeFiles(mupixie)

        deallocate(mupixie)

      end subroutine get_global_data_from_file

c     find_beta_3d
c     ##############################################################
      subroutine find_beta_3d

c     --------------------------------------------------------------
c     Finds periodic orbit beta and stores in b0 array
c     --------------------------------------------------------------

#if defined(FPA)
        use fixed_point_accelerator
#endif

        real(8) :: beta_err,lberr

#if defined(FPA)
        type (fpa_state) :: fpav
#endif
        INTERFACE
          subroutine lsode_3d_beta(neq,t,y,ydot)
          integer :: neq
          real(8) :: t,y(neq),ydot(neq)
          end subroutine lsode_3d_beta
        END INTERFACE

        real(8) :: x(nx*ny*nz),r(nx*ny*nz)

c     Begin program

        if (.not.b_wgt) then
          b0 = 1d0              !Set beta to unity 
          return
        elseif (.not.have_chi) then
          return
        endif

        if (my_rank == 0) write (*,*) '--> LSODE 3D beta integral...'

c     Setup orbits

          !Orbit initialization
          call orbit_setup(nxg+2,nyg+2,nzg+2
     .                    ,xx
     .                    ,yy
     .                    ,zz
     .                    ,bg(:,:,:,1)
     .                    ,bg(:,:,:,2)
     .                    ,bg(:,:,:,3)
     .                    ,bcarg(:,:,:,1)
     .                    ,bcarg(:,:,:,2)
     .                    ,bcarg(:,:,:,3)
     .                    ,jacg
     .                    ,carg
     .                    ,bcond
     .                    ,dtime=orbit_dt
     .                    ,atol=orbit_atol       !Field line integration absolute tolerance
     .                    ,rtol=orbit_rtol       !Field line integration relative tolerance
     .                    ,spline_order=spline_order
     .                    ,solen=.false.
     .                    ,lsode_neq = 4         !Number of equations (only integrate field line)
     .                    ,ds_is_Bds=bwgt        !Use B.ds as integration variable
     .                    ,lsode_s_indep=.true.  !Use arclength as independent variable
     .                    )

c      Compute beta

        !Initial guess (b0 = 1/<B^2> by default)
cc        where (chi*dt*sum(bcar(:,:,:,:)**2,dim=4) < CHI_LIMIT)
cc        where (chi*dt/b0 < CHI_LIMIT)
          b0 = 1d0/sum(bcar(:,:,:,:)**2,dim=4) !Set beta_0=1/B^2 for small B
cc        end where

cc          b0 = 1d0 !Normalization test

        !Iteration
#if defined(FPA)
        !Create fixed-point accelerator
        call fpa_create(fpav,x,maxv=5)

        if (my_rank == 0) then
          write (*,*) '  Using FPA acceleration...'
        endif
#endif

        beta_err = 1d0

        !Picard iteration
        x = reshape(b0(1:nx,1:ny,1:nz),shape(x))

        do while (beta_err > 1d-2)

          beta_err = 0d0

          do k=1,nz
            do j=1,ny
              do i=1,nx

                call getMGmap(i,j,k,igrid,igrid,igrid,ig,jg,kg)

                x0 = grid_params%xx(ig)
                y0 = grid_params%yy(jg)
                z0 = grid_params%zz(kg)

                ct0 = chi(i,j,k)*dt/b0(i,j,k)

cc                length_s  = 0d0    !Initialize orbit length
                b0_s = 0d0

                call orbit_find(iout-2,x0,y0,z0,.true.
     .                         ,lsode_f=lsode_3d_beta
     .                         ,lsode_sum=lsode_out)

                b0_s = lsode_out(1) !Store LSODE accumulator

                call orbit_find(iout-2,x0,y0,z0,.false.
     .                         ,lsode_f=lsode_3d_beta
     .                         ,lsode_sum=lsode_out)

                b0_s = b0_s + lsode_out(1)  !Store LSODE accumulator

                if ((my_rank == 0).and.(iout == 1))
     .                      write (*,FMT="(a)",ADVANCE="NO") '+'

cc                ll(i,j,k) = length_s                  !Find int(dlambda)=int(B.ds)
cc                if (b_wgt) b0(i,j,k) = b0_s/length_s  !Find beta
                b0(i,j,k) = b0_s  !Find beta_0
              enddo
            enddo
          enddo

          !Compute residual vector
          r = x - reshape(b0(1:nx,1:ny,1:nz),shape(x))

          !Convergence test
          beta_err = sqrt(dot2(r/x,r/x)/nxg/nyg/nzg)

          if (my_rank == 0) write (*,*) '  >> Beta error=',beta_err

#if defined(FPA)
          !Apply fixed-point accelerator
          call fpa_correction(fpav,r,dp=dot2)

          !Update solution
          x = x - r

          b0(1:nx,1:ny,1:nz) = reshape(x,shape(b0(1:nx,1:ny,1:nz)))
#else
          x = reshape(b0(1:nx,1:ny,1:nz),shape(x))
#endif
        enddo

c     Destroy orbit

        call orbit_destroy

#if defined(FPA)
        !Destroy fixed-point accelerator
        call fpa_destroy(fpav)
#endif

      end subroutine find_beta_3d

      end subroutine lag_3d_init

c     lag_xport
c     ##############################################################
      subroutine lag_xport(igrid,in_flds,out_flds,chi,iout)

      implicit none

c     --------------------------------------------------------------
c     Performs Lagrangian temperature integration according to B
c     field (bx,by,bz).
c     --------------------------------------------------------------

c     Call variables

        integer :: igrid,iout

        real(8) :: chi(:,:,:),in_flds(:,:,:,:),out_flds(:,:,:,:)

c     Local variables

        integer ::  nx,ny,nz,nints

        !Interface for optional LSODE function evaluation
        INTERFACE
          subroutine lsode_2d_bdf2_3eq(neq,t,y,ydot)
          integer :: neq
          real(8) :: t,y(neq),ydot(neq)
          end subroutine lsode_2d_bdf2_3eq
        END INTERFACE

        INTERFACE
          subroutine lsode_3d_bdf2(neq,t,y,ydot)
          integer :: neq
          real(8) :: t,y(neq),ydot(neq)
          end subroutine lsode_3d_bdf2
        END INTERFACE

c     Begin program

        nints = 8  !Number of LSODE integrals

        nx = grid_params%nxv(igrid)
        ny = grid_params%nyv(igrid)
        nz = grid_params%nzv(igrid)

        if (l2d) then
          call lag_2d(igrid,nx,ny,nz,in_flds,out_flds,nints
     .               ,chi,iout,lsode_f=lsode_2d_bdf2_3eq)
        else
          call lag_3d(igrid,nx,ny,nz,in_flds,out_flds,nints
     .               ,chi,iout,lsode_f=lsode_3d_bdf2)
        endif

      end subroutine lag_xport

c     lag_xport2
c     ##############################################################
      subroutine lag_xport2(igrid,in_flds,out_flds,chi,iout)

      implicit none

c     --------------------------------------------------------------
c     Performs Lagrangian temperature integration according to B
c     field (bx,by,bz).
c     --------------------------------------------------------------

c     Call variables

        integer :: igrid,iout

        real(8) :: chi(:,:,:),in_flds(:,:,:,:),out_flds(:,:,:,:)

c     Local variables

        integer ::  nx,ny,nz,nints

        !Interface for optional LSODE function evaluation
        INTERFACE
          subroutine lsode_2d_bdf2_3eq_dT(neq,t,y,ydot)
          integer :: neq
          real(8) :: t,y(neq),ydot(neq)
          end subroutine lsode_2d_bdf2_3eq_dT
        END INTERFACE

        INTERFACE
          subroutine lsode_3d_bdf2_dT(neq,t,y,ydot)
          integer :: neq
          real(8) :: t,y(neq),ydot(neq)
          end subroutine lsode_3d_bdf2_dT
        END INTERFACE

c     Begin program

        nints = 8  !Number of LSODE integrals

        nx = grid_params%nxv(igrid)
        ny = grid_params%nyv(igrid)
        nz = grid_params%nzv(igrid)

        if (l2d) then
          call lag_2d(igrid,nx,ny,nz,in_flds,out_flds,nints
     .               ,chi,iout,lsode_f=lsode_2d_bdf2_3eq_dT)
        else
          call lag_3d(igrid,nx,ny,nz,in_flds,out_flds,nints
     .               ,chi,iout,lsode_f=lsode_3d_bdf2_dT)
        endif

      end subroutine lag_xport2

c     lag_xport3
c     ##############################################################
      subroutine lag_xport3(igrid,in_flds,out_flds,chi,iout)

      implicit none

c     --------------------------------------------------------------
c     Performs Lagrangian temperature integration according to B
c     field (bx,by,bz).
c     --------------------------------------------------------------

c     Call variables

        integer :: igrid,iout

        real(8) :: chi(:,:,:),in_flds(:,:,:,:),out_flds(:,:,:,:)

c     Local variables

        integer ::  nx,ny,nz

c     Begin program

        nx = grid_params%nxv(igrid)
        ny = grid_params%nyv(igrid)
        nz = grid_params%nzv(igrid)

        if (l2d) then
          call lag_2d_new(igrid,nx,ny,nz,in_flds,out_flds,chi,iout)
        else
          call lag_3d_new(igrid,nx,ny,nz,in_flds,out_flds,chi,iout)
        endif

      end subroutine lag_xport3

c     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c     PUBLIC ROUTINES
c     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

c     lag_update
c     ###############################################################
      subroutine lag_update(igrid,ivar,bcnd,tnm,tn,tnp,chi,iout)

c     ---------------------------------------------------------------
c     Performs Lagrangian update of temperature at time n+1 (returned
c     in tnp) to enforce parallel transport.
c
c     >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
c     MUST BE PRECEDED BY A CALL TO lag_init
c     WORKSPACE FREED  BY A CALL TO lag_destroy
c     >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
c
c     On call sequence:
c       * igrid: grid level
c       * ivar: variable identification
c       * bcnd: boundary conditions
c       * tnm,tn,tnp: 3-level temperature field (for BDF2); returns
c                     lagrangian step result in tnp
c       * chi: parallel transport coeff (array)
c       * iout: output diagnostics level
c     ---------------------------------------------------------------

      implicit none

c     Call variables

      integer :: igrid,iout,bcnd(6),ivar
      real(8),dimension(0:,0:,0:) :: tnm,tn,tnp,chi

c     Local variables

      integer :: nx,ny,nz
      real(8) :: gamma,sgamma,dT_error,b0_diag

      real(8),dimension(:,:,:,:),pointer :: flds,integrals

      real(8),dimension(:,:,:),pointer :: chi_local,rhs,delta_T,zeros

c     Begin program

      if (.not.init_lag) then
        call pstop('lag_update'
     .            ,'Lagrangian integrator uninitialized')
      endif

      nx = size(tnp,1)-2
      ny = size(tnp,2)-2
      nz = size(tnp,3)-2

c     Lagrangian step

      allocate(flds     (0:nx+1,0:ny+1,0:nz+1,3)
     .        ,integrals(0:nx+1,0:ny+1,0:nz+1,5)
     .        ,chi_local(0:nx+1,0:ny+1,0:nz+1)
     .        ,delta_T  (0:nx+1,0:ny+1,0:nz+1)
     .        ,zeros    (0:nx+1,0:ny+1,0:nz+1)
     .        ,rhs      (0:nx+1,0:ny+1,0:nz+1)  )

      zeros = 0d0
      chi_local = 0d0
      flds = 0d0
      integrals = 0d0
      delta_T = 0d0
      rhs = 0d0

c     Find rhs of delta_T iteration
      
      flds(:,:,:,1) = tn
      flds(:,:,:,2) = tnm
      flds(:,:,:,3) = cnp*tnp + cn*tn + cnm*tnm

      where (b0 > 0d0) chi_local = chi/b0  !Divide by beta

      call lag_xport(igrid,flds,integrals,chi_local,iout)

      if (bdf2) then
        sgamma = dt/(dt+dtold)
        gamma = sgamma*sgamma

        where (chi_local /= 0d0) 
          rhs =(       (integrals(:,:,:,1)-tn )  !G(Tn,dt/eps)
     .         - gamma*(integrals(:,:,:,2)-tnm)  !G(Tnm,2dt/eps)
     .         +(       integrals(:,:,:,3)       !P(dT/dt,dt/eps)/gamma/beta
     .          -sgamma*integrals(:,:,:,4))/b0 ) !P(dT/dt,2dt/eps)/beta
     .         /(1d0-sgamma)
        end where
      else
        where (chi_local /= 0d0) 
          rhs = integrals(:,:,:,1)-tn + integrals(:,:,:,3)/b0 !Theta scheme
        end where
      endif

c     Perform delta_T iteration

      delta_T = rhs

      if (.not.bwgt) then
        tnp = (delta_T - cn*tn - cnm*tnm)/cnp
        deallocate(flds,integrals,chi_local,delta_T,rhs)
        return
      endif

cc      if (my_rank == 0) write (*,*) sgamma,cnp,cn,cnm

      if (my_rank ==0) write (*,*) 'Peforming Delta_T iteration...'

      dT_error = 1d0
      do while (dT_error > 1e-6) 
        call setBC(ivar,nx,ny,nz,delta_T,zeros,bcnd
     .            ,igrid,igrid,igrid)

        flds(:,:,:,1) = delta_T
        flds(:,:,:,2) = delta_T
        flds(:,:,:,3) = delta_T

        call lag_xport2(igrid,flds,integrals,chi_local,iout)
      
        if (bdf2) then
          where (chi_local /= 0d0) 
            delta_T = rhs
     .         +(        integrals(:,:,:,1)      !P(dT,dt/eps) 
     .           -sgamma*integrals(:,:,:,2)      !P(dT,2dt/eps)
     .          -(       integrals(:,:,:,3)      !P(dT/B^2,dt/eps)
     .           -sgamma*integrals(:,:,:,4))/b0 )!P(dT/B^2,2dt/eps)
     .          /(1d0-sgamma)
          end where
        else
          where (chi_local /= 0d0) 
            delta_T = rhs + integrals(:,:,:,1)
     .                    - integrals(:,:,:,3)/b0 !Theta scheme
          end where
        endif

        !Compute error
        flds(:,:,:,1) = (delta_T - flds(:,:,:,1))**2

        dT_error = integral(igrid,igrid,igrid,nx,ny,nz,flds(:,:,:,1))
        dT_error = sqrt(dT_error)

        if (my_rank ==0) write (*,*) 'Delta_T error=',dT_error
      enddo

      tnp = (delta_T - cn*tn - cnm*tnm)/cnp

c     Update beta for next time step

      if (bwgt) then
        where (b0 > 0d0) b0 = (1 - integrals(:,:,:,5)/b0)**2
        b0_diag = integral(igrid,igrid,igrid,nx,ny,nz,b0)
        b0_diag = sqrt(b0_diag)

        if (my_rank == 0) then
          write (*,*) '>>>>>Time=',time,' Beta_diff=',b0_diag
        endif

        b0 = integrals(:,:,:,5)
      endif

c     Free workspace

      deallocate(flds,integrals,chi_local,delta_T,rhs)

      end subroutine lag_update

c     lag_residual
c     ###############################################################
      subroutine lag_residual(igrid,tnm,tn,tnp,src,chi,iout)

c     ---------------------------------------------------------------
c     Computes residual of lagrangian integro-differential formulation
c     for anisotropic heat transport.
c
c     >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
c     MUST BE PRECEEDED BY A CALL TO lag_init
c     WORKSPACE FREED  BY A CALL TO lag_destroy
c     >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
c
c     On call sequence:
c       * igrid: grid level
c       * chi: parallel transport coeff
c       * tnm,tn,tnp: 3-level temperature field (for BDF2)
c       * src: on input, source to be propagated by P propagator;
c              on output, returns lagrangian residual
c       * iout: output diagnostics level
c     ---------------------------------------------------------------

      implicit none

c     Call variables

      integer :: igrid,iout
      real(8),dimension(0:,0:,0:) :: tnm,tn,tnp,chi,src

c     Local variables

      integer :: nx,ny,nz
      real(8) :: gamma,sgamma,dT_error,b0_diag

      real(8),dimension(:,:,:,:),pointer :: flds,integrals

      real(8),dimension(:,:,:),pointer :: chi_local,delta_T

c     Begin program

      if (.not.init_lag) then
        call pstop('lag_update'
     .            ,'Lagrangian integrator uninitialized')
      endif

      nx = size(tnp,1)-2
      ny = size(tnp,2)-2
      nz = size(tnp,3)-2

c     Lagrangian step

      allocate(flds     (0:nx+1,0:ny+1,0:nz+1,4)
     .        ,integrals(0:nx+1,0:ny+1,0:nz+1,8)
     .        ,chi_local(0:nx+1,0:ny+1,0:nz+1)
     .        ,delta_T  (0:nx+1,0:ny+1,0:nz+1))

c     Find rhs

      delta_T = cnp*tnp + cn*tn + cnm*tnm

      flds(:,:,:,1) = tn
      flds(:,:,:,2) = tnm
      flds(:,:,:,3) = delta_T
      flds(:,:,:,4) = src

      chi_local = 0d0
      where (b0 > 0d0) chi_local = chi/b0  !Divide by beta

      call lag_xport3(igrid,flds,integrals,chi_local,iout)

      if (bdf2) then
        sgamma = dt/(dt+dtold)
        gamma = sgamma*sgamma

        where (chi_local /= 0d0) 
          src =(       (integrals(:,:,:,1)-tn )  !G(Tn,dt/eps)
     .         - gamma*(integrals(:,:,:,2)-tnm)  !G(Tnm,2dt/eps)
     .         +(       integrals(:,:,:,3)       !P(dT/dt,dt/eps)/gamma/beta
     .          -sgamma*integrals(:,:,:,4))*dt/b0!P(dT/dt,2dt/eps)/beta
     .         +        integrals(:,:,:,5)       !P(dT,dt/eps) 
     .          -sgamma*integrals(:,:,:,6)       !P(dT,2dt/eps)
     .         -(       integrals(:,:,:,7)       !P(dT/B^2,dt/eps)
     .          -sgamma*integrals(:,:,:,8))/b0 ) !P(dT/B^2,2dt/eps)
     .          /(1d0-sgamma)
        end where
      else  !Theta scheme
        where (chi_local /= 0d0) 
          src = integrals(:,:,:,1)-tn + dt*integrals(:,:,:,3)/b0
     .        + integrals(:,:,:,5)       - integrals(:,:,:,7)/b0
        end where
      endif

c     Divide by dt to provide RHS to temporal term

      src = src/dt

ccc     Update beta for next time step
cc
cc      if (bwgt) then
cc        b0 = (b0 - integrals(:,:,:,5))**2
cc        b0_diag = integral(igrid,igrid,igrid,nx,ny,nz,b0)
cc        b0_diag = sqrt(b0_diag)
cc
cc        if (my_rank == 0) then
cc          write (*,*) '>>>>>Time=',time,' Beta_diff=',b0_diag
cc        endif
cc
cc        b0 = integrals(:,:,:,5)
cc      endif

c     Free workspace

      deallocate(flds,integrals,chi_local,delta_T)

      end subroutine lag_residual

c     lag_init
c     ##############################################################
      subroutine lag_init(igr,nx,ny,nz,bb,alp,iout,G_tol,s_order
     .                   ,B_wgt,three_d,psi,chi,poinc_file)

      implicit none

c     --------------------------------------------------------------
c     Initializes Lagrangian step.
c     --------------------------------------------------------------

c     Call variables

      integer :: nx,ny,nz,igr,s_order,iout
 
      real(8) :: G_tol,alp
      real(8),INTENT(IN) :: bb(0:nx+1,0:ny+1,0:nz+1,3)

      logical :: B_wgt,three_d

      real(8),INTENT(IN),optional :: psi(0:nx+1,0:ny+1,0:nz+1)
     .                              ,chi(0:nx+1,0:ny+1,0:nz+1)

      character(*),optional :: poinc_file

c     Local variables

      integer :: nxg,nyg,nzg
      logical :: lag3d

c     Begin program

      nxg = grid_params%nxgl(igr)
      nyg = grid_params%nygl(igr)
      nzg = grid_params%nzgl(igr)

      lag3d = (nxg > 1 .and. nyg > 1 .and. nzg > 1)
     .       .or.(three_d)
     .       .or.PRESENT(poinc_file)

      if (lag3d) then
        call lag_3d_init(igr,nx,ny,nz,bb,alp,iout
     .                  ,G_tol,s_order,B_wgt,chi=chi
     .                  ,poinc_file=poinc_file)
      else
        call lag_2d_init(igr,nx,ny,nz,bb,alp,iout
     .                  ,G_tol,s_order,B_wgt,psi=psi
     .                  ,chi=chi)
      endif

c     End program

      end subroutine lag_init

c     lag_destroy
c     ##############################################################
      subroutine lag_destroy

      implicit none

c     --------------------------------------------------------------
c     Destroys storage for Lagrangian step
c     --------------------------------------------------------------

c     Call variables

c     Local variables

      integer :: alloc_stat

c     Begin program

      call orbit_destroy

      deallocate(ll,stat=alloc_stat)
      deallocate(b0,stat=alloc_stat)

c     End program

      end subroutine lag_destroy

ccc     lagrangian_init_from_file
ccc     ##############################################################
cc      subroutine lagrangian_init_from_file(file,igrid,nx,ny,nz,bb,alp
cc     .                                    ,G_tol,spline_order,nints
cc     .                                    ,B_wgt,chi)
cc
cc      use iosetup
cc
cc      implicit none
cc
ccc     --------------------------------------------------------------
ccc     Initializes Lagrangian step from PIXIE3D's poincare output file.
ccc     Provides bb as output
ccc     --------------------------------------------------------------
cc
ccc     Call variables
cc
cc        character(*) :: file
cc
cc        integer,INTENT(IN)  :: nx,ny,nz,igrid,spline_order,nints
cc
cc        real(8) :: G_tol,alp
cc        real(8),pointer,dimension(:,:,:,:) :: bb
cc
cc        logical :: B_wgt
cc
cc        real(8),INTENT(IN),optional :: chi(0:nx+1,0:ny+1,0:nz+1)
cc
ccc     Local variables
cc
cc        integer :: upixie=12345,ierr,nfiles,i
cc
cc        integer :: nxg,nyg,nzg,bcnd(6)
cc
cc        real(8),pointer,dimension(:)       :: xx,yy,zz
cc        real(8),pointer,dimension(:,:,:)   :: jac,jacg
cc        real(8),pointer,dimension(:,:,:,:) :: bcar,bg,bcarg,carg,car
cc
cc        logical :: tor
cc
cc        logical :: have_chi
cc
cc        integer,pointer,dimension(:) :: mupixie
cc
ccc     Begin program
cc
cc        init_lag = .true.
cc
cc        l2d = .false.   !Indicate 3D treatment
cc
cc        alpha = alp
cc
cc        bwgt = B_wgt
cc
cc        have_chi = PRESENT(chi)
cc
ccc     Compute Green's integral cut-offs
cc
cc        t_tol = G_tol !Initialize Green's function tolerance
cc
cc        T_cutoff2 = 1d0
cc        S_cutoff2 = 1d0
cc        do i=1,10
cc          T_cutoff2 = abs(log(  sqrt(pi)*t_tol*sqrt(T_cutoff2)))
cc          S_cutoff2 = abs(log(2*sqrt(pi)*t_tol*     S_cutoff2 ))
cc        enddo
cc
ccc     Check file configuration
cc
cc        ierr = openFiles(my_rank,np,upixie,trim(file),mupixie)
cc
cc        if (ierr /= 0) call pstop('openGraphicsfiles'
cc     .                           ,'Cannot find files')
cc
ccc     Read global data
cc
cc        nfiles = size(mupixie)
cc
cccc        if (nfiles /= np) then
cc          call read_global_data_async
cccc        else
cccc          call read_global_data_sync
cccc        endif
cc
ccc     Close files
cc
cc        call closeFiles(mupixie)
cc
ccc     Find initial beta=1/B^2
cc
cc        !Output
cc        if (my_rank == 0) then
cc          if (iout == 1) then
cc            write (*,FMT="(a)",ADVANCE="NO") '3D orbit integrals: '
cc          elseif (iout <= 0) then
cc            if (solenoidal) then
cc              write (*,*) 'Solenoidal 3D orbit integrals...'
cc            else
cc              write (*,*) 'LSODE 3D orbit integrals...'
cc            endif
cc          endif
cc        endif
cc
cc        allocate(b0(0:nx+1,0:ny+1,0:nz+1))
cc
cc        call find_beta
cc
ccc     Orbit initialization on global array
cc
cc        call orbit_setup(nxg+2,nyg+2,nzg+2
cc     .                  ,grid_params%xg
cc     .                  ,grid_params%yg
cc     .                  ,grid_params%zg
cc     .                  ,bg(:,:,:,1)
cc     .                  ,bg(:,:,:,2)
cc     .                  ,bg(:,:,:,3)
cc     .                  ,bcarg(:,:,:,1)
cc     .                  ,bcarg(:,:,:,2)
cc     .                  ,bcarg(:,:,:,3)
cc     .                  ,jacg
cc     .                  ,carg
cc     .                  ,bcond
cc     .                  ,dtime=orbit_dt
cc     .                  ,atol=orbit_atol  !Field line integration absolute tolerance
cc     .                  ,rtol=orbit_rtol  !Field line integration relative tolerance
cc     .                  ,spline_order=spline_order
cc     .                  ,solen=.false.
cc     .                  ,lsode_neq=5 + nints  !5 field eqs + number of field integrals
cc     .                  ,ds_is_Bds=bwgt   !Use B.ds as integration variable
cc     .                  )
cc
ccc     Deallocate memory
cc
cc        deallocate(bg,bcarg,jacg,carg,xx,yy,zz)
cc
ccc     End program
cc
cc      contains
cc
ccc     read_global_data_sync
ccc     ##################################################################
cc      subroutine read_global_data_sync
cc
ccc     ------------------------------------------------------------------
ccc     Reads global data from PIXIE3D files in a synchronous manner (i.e.,
ccc     one file per proc).
ccc     ------------------------------------------------------------------
cc
cc        implicit none
cc
ccc     Call variables
cc
ccc     Local variables
cc
cc        integer :: nxl,nyl,nzl,ilmg,ihpg,jlmg,jhpg,klmg,khpg,funit
cc
ccc     Begin program
cc
cc        if (my_rank == 0) then
cc          write (*,*) 'Synchronous input in Lagrangian step'
cc        endif
cc
cc        funit = mupixie(1)
cc
cc        read (funit) nxg,nyg,nzg,bcnd,tor
cc
cc        if (tor .and. (coords /= 'vmc' .and. coords /= 'tor')) then
cc          call closeFiles(mupixie)
cc          call pstop('lagrangian_init','Coordinates are not toroidal')
cc        endif
cc
cc        if (sum(bcnd-bcond) /= 0) then
cc          call closeFiles(mupixie)
cc          call pstop('lagrangian_init','BC conf. does not agree')
cc        endif
cc
cc        if (nxd /= nxg .or. nyd /= nyg .or. nzd /= nzg) then
cc          call closeFiles(mupixie)
cc          write (*,*) 'Grid size=',nxg,nyg,nzg
cc          call pstop('lagrangian_init'
cc     .              ,'Global grid size does not agree')
cc        endif
cc
ccc     Read global position arrays
cc
cc        allocate(xx(nxg+2),yy(nyg+2),zz(nzg+2))
cc
cc        read (funit) xx
cc        read (funit) yy
cc        read (funit) zz
cc
ccc     Local limits
cc
cc        if (np > 1) then
cc          read (funit) ilmg,ihpg,jlmg,jhpg,klmg,khpg
cc        else
cc          ilmg = 0
cc          ihpg = nxg+1
cc
cc          jlmg = 0
cc          jhpg = nyg+1
cc
cc          klmg = 0
cc          khpg = nzg+1
cc        endif
cc
cc        nxl = ihpg-ilmg-1
cc        nyl = jhpg-jlmg-1
cc        nzl = khpg-klmg-1
cc
cc        if (nx /= nxl .or. nyl /= nyl .or. nz /= nzl) then
cc          write (*,*) 'Local grid size=',nxl,nyl,nzl
cc          call pstop('lagrangian_init'
cc     .              ,'Local grid size does not agree')
cc        endif
cc
ccc     Allocate local arrays
cc
cc        allocate(bcar(0:nx+1,0:ny+1,0:nz+1,3)
cc     .          ,jac (0:nx+1,0:ny+1,0:nz+1)
cc     .          ,car (0:nx+1,0:ny+1,0:nz+1,3))
cc
ccc     Fill local arrays
cc
cc        read (funit) bb(:,:,:,1)
cc        read (funit) bb(:,:,:,2)
cc        read (funit) bb(:,:,:,3)
cc        read (funit) bcar(:,:,:,1)
cc        read (funit) bcar(:,:,:,2)
cc        read (funit) bcar(:,:,:,3)
cc        read (funit) jac
cc        read (funit) car
cc
ccc     Allocate global arrays
cc
cc        allocate(bg   (0:nxg+1,0:nyg+1,0:nzg+1,3)
cc     .          ,bcarg(0:nxg+1,0:nyg+1,0:nzg+1,3)
cc     .          ,carg (0:nxg+1,0:nyg+1,0:nzg+1,3)
cc     .          ,jacg (0:nxg+1,0:nyg+1,0:nzg+1))
cc
ccc     Fill global arrays
cc
cc        call find_global(bb,bg)
cc        call find_global(bcar,bcarg)
cc
cc        call find_global(jac,jacg)
cc        call find_global(car,carg)
cc
ccc     Update grid metrics
cc
cc        call defineGridMetric(grid_params,xcar=car)
cc
cc        gmetric => grid_params%gmetric
cc
ccc     Deallocate memory
cc
cc        deallocate(bcar,jac,car)
cc
cc      end subroutine read_global_data_sync
cc
ccc     read_global_data_async
ccc     ##################################################################
cc      subroutine read_global_data_async
cc
ccc     ------------------------------------------------------------------
ccc     Reads global data from PIXIE3D files in an asynchronous manner (i.e.,
ccc     global solution is found in all procs by reading multiple files).
ccc     ------------------------------------------------------------------
cc
cc        implicit none
cc
ccc     Call variables
cc
ccc     Local variables
cc
cc        integer :: ifl,ilmg,ihpg,jlmg,jhpg,klmg,khpg,nx,ny,nz,nfiles
cc
ccc     Begin program
cc
cc        if (my_rank == 0) then
cc          write (*,*) 'Asynchronous input in Lagrangian step'
cc        endif
cc
ccc     Read setup info
cc
cc        nfiles = size(mupixie)
cc
cc        do ifl=1,nfiles
cc          read (mupixie(ifl)) nxg,nyg,nzg,bcnd,tor
cc        enddo
cc
cc        if (tor .and. (coords /= 'vmc' .and. coords /= 'tor')) then
cc          call closeFiles(mupixie)
cc          call pstop('lagrangian_init','Coordinates are not toroidal')
cc        endif
cc
cc        if (sum(bcnd-bcond) /= 0) then
cc          call closeFiles(mupixie)
cc          call pstop('lagrangian_init','BC conf. does not agree')
cc        endif
cc
cc        if (nxd /= nxg .or. nyd /= nyg .or. nzd /= nzg) then
cc          call closeFiles(mupixie)
cc          write (*,*) 'Grid size=',nxg,nyg,nzg
cc          call pstop('lagrangian_init'
cc     .              ,'Global grid size does not agree')
cc        endif
cc
ccc     Allocate global arrays
cc
cc        allocate(xx(nxg+2),yy(nyg+2),zz(nzg+2))
cc
cc        allocate(bg   (0:nxg+1,0:nyg+1,0:nzg+1,3)
cc     .          ,bcarg(0:nxg+1,0:nyg+1,0:nzg+1,3)
cc     .          ,carg (0:nxg+1,0:nyg+1,0:nzg+1,3)
cc     .          ,jacg (0:nxg+1,0:nyg+1,0:nzg+1  ))
cc
ccc     Read global arrays
cc
cc        do ifl=1,nfiles
cc
cc          read (mupixie(ifl)) xx
cc          read (mupixie(ifl)) yy
cc          read (mupixie(ifl)) zz
cc
cc          if (nfiles > 1) then
cc            read (mupixie(ifl)) ilmg,ihpg,jlmg,jhpg,klmg,khpg
cc          else
cc            ilmg = 0
cc            ihpg = nxg+1
cc
cc            jlmg = 0
cc            jhpg = nyg+1
cc
cc            klmg = 0
cc            khpg = nzg+1
cc          endif
cc
cc          read (mupixie(ifl)) bg   (ilmg:ihpg,jlmg:jhpg,klmg:khpg,1)
cc          read (mupixie(ifl)) bg   (ilmg:ihpg,jlmg:jhpg,klmg:khpg,2)
cc          read (mupixie(ifl)) bg   (ilmg:ihpg,jlmg:jhpg,klmg:khpg,3)
cc          read (mupixie(ifl)) bcarg(ilmg:ihpg,jlmg:jhpg,klmg:khpg,1)
cc          read (mupixie(ifl)) bcarg(ilmg:ihpg,jlmg:jhpg,klmg:khpg,2)
cc          read (mupixie(ifl)) bcarg(ilmg:ihpg,jlmg:jhpg,klmg:khpg,3)
cc          read (mupixie(ifl)) jacg (ilmg:ihpg,jlmg:jhpg,klmg:khpg)
cc          read (mupixie(ifl)) carg (ilmg:ihpg,jlmg:jhpg,klmg:khpg,3)
cc
cc        enddo
cc
ccc     Assign local B arrays
cc
cc        bb = bg(ilomg:ihipg,jlomg:jhipg,klomg:khipg,:)
cc
ccc     Update grid metrics
cc
cc        nxl = ihipg-ilomg-1
cc        nyl = jhipg-jlomg-1
cc        nzl = khipg-klomg-1
cc
cc        call defineGridMetric(grid_params
cc     .           ,xcar=carg(ilomg:ihipg,jlomg:jhipg,klomg:khipg,:))
cc
cc        gmetric => grid_params%gmetric
cc
ccc     Deallocate memory
cc
cc      end subroutine read_global_data_async
cc
cc      end subroutine lagrangian_init_from_file

c     per_orbit_init
c     ###############################################################
      subroutine per_orbit_init(x0,x,y0,y,z0,z)

c     ---------------------------------------------------------------
c     Initializes periodic orbit check.
c     ---------------------------------------------------------------

        implicit none

c     Call variables

        real(8) :: x0,x,y0,y,z0,z

c     Local variables

c     Begin program

        orb_init = .false.

c       Determine initial orbit updates

        dx0 = x-x0 !Determines initial direction of orbit
        dy0 = y-y0 !Determines initial direction of orbit
        dz0 = z-z0 !Determines initial direction of orbit

        !Exclude ignorable directions
        if (grid_params%nxgl(1) == 1) dx0 = 0d0
        if (grid_params%nygl(1) == 1) dy0 = 0d0
        if (grid_params%nzgl(1) == 1) dz0 = 0d0

c       Select periodic direction based on largest displacement from surface

        x_per = (abs(dx0) > abs(dy0)) .and. (abs(dx0) > abs(dz0))
        y_per = (abs(dy0) > abs(dx0)) .and. (abs(dy0) > abs(dz0))
        z_per = (abs(dz0) > abs(dy0)) .and. (abs(dz0) > abs(dx0))

c       Set logic to compute lengths using surface of sections with high-order interpolation

        call ss_init(4,x_per,x0,y_per,y0,z_per,z0)

      end subroutine per_orbit_init

c     per_orbit_chk
c     ###############################################################
      subroutine per_orbit_chk(so,ss,xo,x,yo,y,zo,z,ts,xs,ys,zs)

c     ---------------------------------------------------------------
c     Check for end of periodic orbit. Returns ts <> 0d0 if end is
c     found.
c     ---------------------------------------------------------------

        implicit none

c     Call variables

        real(8) :: so,ss,xo,x,yo,y,zo,z,xs,ys,zs,ts

c     Local variables

        real(8) :: xoo,yoo,zoo

c     Begin program

        ts = 0d0

        if (    (x_per.and.(x-xo)*dx0 > 0d0)
     .      .or.(y_per.and.(y-yo)*dy0 > 0d0)
     .      .or.(z_per.and.(z-zo)*dz0 > 0d0)) then

          call sschk_ho(xs,ys,zs,ts)  !--> Returns ts <> 0 if intersects
        endif

      end subroutine per_orbit_chk

ccc     lagrangian_T
ccc     ##############################################################
cc      subroutine lagrangian_T(igrid,tt,chi,iout)
cc
cc      implicit none
cc
ccc     --------------------------------------------------------------
ccc     Performs Lagrangian temperature integration according to B
ccc     field (bx,by,bz).
ccc     --------------------------------------------------------------
cc
ccc     Call variables
cc
cc        integer :: igrid,iout
cc
cc        real(8) :: chi,tt(:,:,:)
cc
ccc     Local variables
cc
cc        integer ::  nx,ny,nz
cc
cc        INTERFACE
cc          function t_int_per(fwd,xo,x,yo,y,zo,z,so,ss)
cc          logical :: t_int_per,fwd
cc          real(8) :: xo,x,yo,y,zo,z,so,ss
cc          end function t_int_per
cc        END INTERFACE
cc
cc        !Interface for optional LSODE function evaluation
cc        INTERFACE
cc          subroutine lsode_3d_T(neq,t,y,ydot)
cc          integer :: neq
cc          real(8) :: t,y(neq),ydot(neq)
cc          end subroutine lsode_3d_T
cc        END INTERFACE
cc
cc        INTERFACE
cc          subroutine lsode_2d_T_3eq(neq,t,y,ydot)
cc          integer :: neq
cc          real(8) :: t,y(neq),ydot(neq)
cc          end subroutine lsode_2d_T_3eq
cc        END INTERFACE
cc
ccc     Begin program
cc
cc        nx = grid_params%nxv(igrid)
cc        ny = grid_params%nyv(igrid)
cc        nz = grid_params%nzv(igrid)
cc
cc        if (l2d) then
cc          if (solenoidal) then
cc            call lagrangian_2d(igrid,nx,ny,nz,tt
cc     .                        ,chi,iout,int_f=t_int_per)
cc          else
cc            call lagrangian_2d(igrid,nx,ny,nz,tt
cc     .                        ,chi,iout,lsode_f=lsode_2d_T_3eq)
cc          endif
cc        else
cc          call lagrangian_3d(igrid,nx,ny,nz,tt
cc     .                      ,chi,iout,lsode_f=lsode_3d_T)
cc        endif
cc
cc      end subroutine lagrangian_T
cc
ccc     lagrangian_S
ccc     ##############################################################
cc      subroutine lagrangian_S(igrid,ss,chi,iout)
cc
cc      implicit none
cc
ccc     --------------------------------------------------------------
ccc     Performs Lagrangian temperature integration according to B
ccc     field (bx,by,bz).
ccc     --------------------------------------------------------------
cc
ccc     Call variables
cc
cc        integer :: igrid,iout
cc
cc        real(8) :: chi,ss(:,:,:)
cc
ccc     Local variables
cc
cc        integer :: nx,ny,nz
cc
cc        INTERFACE
cc          function s_int_per(fwd,xo,x,yo,y,zo,z,so,ss)
cc          logical :: s_int_per,fwd
cc          real(8) :: xo,x,yo,y,zo,z,so,ss
cc          end function s_int_per
cc        END INTERFACE
cc
cc        !Interface for optional LSODE function evaluation
cc        INTERFACE
cc          subroutine lsode_3d_S(neq,t,y,ydot)
cc          integer :: neq
cc          real(8) :: t,y(neq),ydot(neq)
cc          end subroutine lsode_3d_S
cc        END INTERFACE
cc
cc        INTERFACE
cc          subroutine lsode_2d_S_3eq(neq,t,y,ydot)
cc          integer :: neq
cc          real(8) :: t,y(neq),ydot(neq)
cc          end subroutine lsode_2d_S_3eq
cc        END INTERFACE
cc
ccc     Begin program
cc
cc        nx = grid_params%nxv(igrid)
cc        ny = grid_params%nyv(igrid)
cc        nz = grid_params%nzv(igrid)
cc
cc        if (l2d) then
cc          if (solenoidal) then
cc            call lagrangian_2d(igrid,nx,ny,nz,ss
cc     .                        ,chi,iout,int_f=s_int_per)
cc          else
cc            call lagrangian_2d(igrid,nx,ny,nz,ss
cc     .                        ,chi,iout,lsode_f=lsode_2d_S_3eq)
cc          endif
cc        else
cc          call lagrangian_3d(igrid,nx,ny,nz,ss
cc     .                      ,chi,iout,lsode_f=lsode_3d_S)
cc
cc        endif
cc
cc      end subroutine lagrangian_S

      end module par_int

c l_int
c ###################################################################
      function l_int(fwd,xo,x,yo,y,zo,z,so,ss) result(exloop)

c -------------------------------------------------------------------
c     Performs integral of arc length. On input, old and new logical
c     positions and field line length
c -------------------------------------------------------------------

      use par_int

      implicit none

c Call variables

      logical :: exloop,fwd
      real(8) :: xo,x,yo,y,zo,z,so,ss

c Local variables

      integer :: ierr
      real(8) :: xh,yh,zh,sh,bb,bx,by,bz,xs,ys,zs,ts,b2,ds

c Begin program

      if (orb_init) call per_orbit_init(xo,x,yo,y,zo,z)

c Check whether orbit should terminate

      call per_orbit_chk(so,ss,xo,x,yo,y,zo,z,ts,xs,ys,zs) !--> Returns ts <> 0 if intersects

      exloop = (ts /= 0d0)

c Update current position

      if (exloop) then
        x = xs ; y = ys ; z = zs
      else
        ts = ss
      endif

      length_s = ts

c Add B-weighed arc-length

      ds = (ts-so)

      xh = .5*(xo+x)
      yh = .5*(yo+y)
      zh = .5*(zo+z)

      !Cartesian components
      call getB(xh,yh,zh,bx,by,bz,.false.,.true.,ierr)

      b2 = bx**2 + by**2 + bz**2

      b0_s = b0_s + ds/b2       !Integrate dlambda/B^2

      end function l_int

ccc t_int
ccc ###################################################################
cc      function t_int(fwd,xo,x,yo,y,zo,z,so,ss) result(ex_loop)
cc
ccc -------------------------------------------------------------------
ccc     Performs line integral of temperature.
ccc -------------------------------------------------------------------
cc
cc        use par_int
cc
cc        implicit none
cc
ccc Call variables
cc
cc        logical :: ex_loop,fwd
cc        real(8) :: xo,x,yo,y,zo,z,so,ss
cc
ccc Local variables
cc
cc        integer :: ierror
cc        real(8) :: xh,yh,zh,sh,tt,tk,it_max=0d0
cc        save :: it_max
cc
ccc Begin program
cc
ccc Find middle points
cc
cc        xh = 0.5*(x +xo)
cc        yh = 0.5*(y +yo)
cc        zh = 0.5*(z +zo)
cc        sh = 0.5*(ss+so)
cc
ccc Find temperature
cc
cc        tt = evalFld(xh,yh,zh,ierror)
cc
ccc Find kernel
cc
cc        tk = tkernel_3d(sh,ct0,alpha,ex_loop)
cc
ccc Integrate temperature
cc
cc        temp_s   = temp_s   + tt*tk*(ss-so)
cc        kernel_s = kernel_s +    tk*(ss-so)
cc
cc      end function t_int
cc
ccc t_int_per
ccc ###################################################################
cc      function t_int_per(fwd,xo,x,yo,y,zo,z,so,ss) result(exloop)
cc
ccc -------------------------------------------------------------------
ccc     Performs line integral of temperature in periodic domain.
ccc -------------------------------------------------------------------
cc
cc      use par_int
cc
cc      implicit none
cc
ccc Call variables
cc
cc      logical :: exloop,fwd
cc      real(8) :: xo,x,yo,y,zo,z,so,ss
cc
ccc Local variables
cc
cc      integer :: ierror
cc      real(8) :: xh,yh,zh,sh,tt,tk,xs,ys,zs,ts
cc
ccc Begin program
cc
cc      exloop = .false.
cc
ccccc Check whether orbit should terminate (x0,y0,z0 is initial orbit position)
cccc
cccc      if (orb_init) call per_orbit_init(xo,x,yo,y,zo,z)
cccc
ccccc Check whether orbit should terminate
cccc
cccc      call per_orbit_chk(so,ss,xo,x,yo,y,zo,z,ts,xs,ys,zs)  !--> Returns ts <> 0 if intersects
cccc
cccc      exloop = (ts /= 0d0)
cccc
ccccc Update current position
cccc
cccc      if (exloop) then
cccc        x = xs ; y = ys ; z = zs
cccc      else
cc        ts = ss
cccc      endif
cc
ccc Find middle points
cc
cc      xh = 0.5*(x +xo)
cc      yh = 0.5*(y +yo)
cc      zh = 0.5*(z +zo)
cc      sh = 0.5*(ts+so)
cc
ccc Find temperature
cc
cc      tt = evalFld(xh,yh,zh,ierror)
cc
ccc Find kernel
cc
cc      tk = tkernel_per(sh,ct0,alpha)
cc
ccc Integrate temperature
cc
cc      temp_s   = temp_s   + tt*tk*(ts-so)
cc      kernel_s = kernel_s +    tk*(ts-so)
cc
cc      end function t_int_per
cc
ccc s_int
ccc ###################################################################
cc      function s_int(fwd,xo,x,yo,y,zo,z,so,ss) result(ex_loop)
cc
ccc -------------------------------------------------------------------
ccc     Performs line integral of temperature.
ccc -------------------------------------------------------------------
cc
cc        use par_int
cc
cc        implicit none
cc
ccc Call variables
cc
cc        logical :: ex_loop,fwd
cc        real(8) :: xo,x,yo,y,zo,z,so,ss
cc
ccc Local variables
cc
cc        integer :: ierror
cc        real(8) :: xh,yh,zh,sh,tt,tk,it_max=0d0
cc        save :: it_max
cc
ccc Begin program
cc
ccc Find middle points
cc
cc        xh = 0.5*(x +xo)
cc        yh = 0.5*(y +yo)
cc        zh = 0.5*(z +zo)
cc        sh = 0.5*(ss+so)
cc
ccc Find temperature
cc
cc        tt = evalFld(xh,yh,zh,ierror)
cc
ccc Find kernel
cc
cc        tk = skernel_3d(sh,ct0,alpha,ex_loop)
cc
ccc Integrate temperature
cc
cc        temp_s   = temp_s   + tt*tk*(ss-so)
cc        kernel_s = kernel_s +    tk*(ss-so)
cc
cc      end function s_int
cc
ccc s_int_per
ccc ###################################################################
cc      function s_int_per(fwd,xo,x,yo,y,zo,z,so,ss) result(exloop)
cc
ccc -------------------------------------------------------------------
ccc     Performs line integral of temperature in periodic domain.
ccc -------------------------------------------------------------------
cc
cc      use par_int
cc
cc      implicit none
cc
ccc Call variables
cc
cc      logical :: exloop,fwd
cc      real(8) :: xo,x,yo,y,zo,z,so,ss
cc
ccc Local variables
cc
cc      integer :: ierror
cc      real(8) :: xh,yh,zh,sh,tt,tk,xs,ys,zs,ts
cc
ccc Begin program
cc
cc      exloop = .false.
cc
ccccc Check whether orbit should terminate
cccc
cccc      if (orb_init) call per_orbit_init(xo,x,yo,y,zo,z)
cccc
ccccc Check whether orbit should terminate
cccc
cccc      call per_orbit_chk(so,ss,xo,x,yo,y,zo,z,ts,xs,ys,zs)  !--> Returns ts <> 0 if intersects
cccc
cccc      exloop = (ts /= 0d0)
cccc
ccccc Update current position
cccc
cccc      if (exloop) then
cccc        x = xs ; y = ys ; z = zs
cccc      else
cc        ts = ss
cccc      endif
cc
ccc Find middle points
cc
cc      xh = 0.5*(x +xo)
cc      yh = 0.5*(y +yo)
cc      zh = 0.5*(z +zo)
cc      sh = 0.5*(ts+so)
cc
ccc Find temperature
cc
cc      tt = evalFld(xh,yh,zh,ierror)
cc
ccc Find kernel
cc
cc      tk = skernel_per(sh,ct0,alpha)
cc
ccc Integrate temperature
cc
cc      temp_s   = temp_s   + tt*tk*(ts-so)
cc      kernel_s = kernel_s +    tk*(ts-so)
cc
cc      end function s_int_per
cc
ccc beta_int
ccc ###################################################################
cc      function beta_int(fwd,xo,x,yo,y,zo,z,so,ss) result(exloop)
cc
ccc -------------------------------------------------------------------
ccc     Performs integral of arc length. On input, old and new logical
ccc     positions and field line length
ccc -------------------------------------------------------------------
cc
cc      use par_int
cc
cc      implicit none
cc
ccc Call variables
cc
cc      logical :: exloop,fwd
cc      real(8) :: xo,x,yo,y,zo,z,so,ss
cc
ccc Local variables
cc
cc      integer :: ierr
cc      real(8) :: xh,yh,zh,sh,bb,bx,by,bz,xs,ys,zs,ts,b2,ds
cc
ccc Begin program
cc
cc      if (orb_init) call per_orbit_init(xo,x,yo,y,zo,z)
cc
ccc Check whether orbit should terminate
cc
cc      call per_orbit_chk(so,ss,xo,x,yo,y,zo,z,ts,xs,ys,zs) !--> Returns ts <> 0 if intersects
cc
cc      exloop = (ts /= 0d0)
cc
ccc Update current position
cc
cc      if (exloop) then
cc        x = xs ; y = ys ; z = zs
cc      else
cc        ts = ss
cc      endif
cc
ccc Add B-weighed arc-length
cc
cc      ds = (ts-so)
cc
cc      xh = .5*(xo+x)
cc      yh = .5*(yo+y)
cc      zh = .5*(zo+z)
cc
cc      !Cartesian components
cc      call getB(xh,yh,zh,bx,by,bz,.false.,.true.,ierr)
cc
cc      b2 = bx**2 + by**2 + bz**2
cc
cccc      b0_s = b0_s + ds/b2       !Integrate dlambda/B^2
cccc      b0_s = b0_s + ds*skernel_per(0.5*(ts+so),ct0)/b2       !Integrate P(1/B^2,tau)
cc      b0_s = b0_s + ds*skernel_per(0.5*(ts+so),ct0)       !Integrate P(1,tau)
cc
cc      end function beta_int
cc
ccc l_int_3d
ccc ###################################################################
cc      function l_int_3d(fwd,xo,x,yo,y,zo,z,so,ss) result(exloop)
cc
ccc -------------------------------------------------------------------
ccc     Performs integral of arc length. On input, old and new logical
ccc     positions and field line length
ccc -------------------------------------------------------------------
cc
cc      use par_int
cc
cc      implicit none
cc
ccc Call variables
cc
cc      logical :: exloop,fwd
cc      real(8) :: xo,x,yo,y,zo,z,so,ss
cc
ccc Local variables
cc
cc      integer :: ierr
cc      real(8) :: xh,yh,zh,sh,bx,by,bz,tst,b2,ds
cc
ccc Begin program
cc
ccc Add arc length
cc
cc      ds = (ss-so)
cc
cc      length_s = length_s + ds  !Integrate dlambda
cc
cc      tst = skernel_3d(ss,ct0,alpha,exloop)
cc
ccc Add B-weighed arc-length
cc
cc      xh = .5*(xo+x)
cc      yh = .5*(yo+y)
cc      zh = .5*(zo+z)
cc
cc      !Cartesian components
cc      call getB(xh,yh,zh,bx,by,bz,.false.,.true.,ierr)
cc
cc      b2 = bx**2 + by**2 + bz**2
cc
cc      b0_s = b0_s + ds*tst/b2       !Integrate dlambda/B^2
cc
cc      end function l_int_3d

c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
c  LSODE DRIVER ROUTINES
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

c lsode_3d_T
c #########################################################################
      subroutine lsode_3d_T(neq,t,y,ydot)

c -------------------------------------------------------------------------
c     Driver routine for LSODE integration: ydot = ff(y,t)
c -------------------------------------------------------------------------

      use par_int

      implicit none

c Call variables

      integer :: neq
      real(8) :: t,y(neq),ydot(neq)

c Local variables

      integer :: ierror
      real(8) :: x1,y1,z1,s1,tt,it_max=0d0
      save :: it_max

c Begin program

c Call LSODE orbit evaluation function (extracting logical coords x1,y1,z1)

      call lsode_orbit_f(5,t,y(1:5),ydot(1:5),xx1=x1,yy1=y1,zz1=z1)

c Call additional integrators

      s1 = abs(y(5))

      tt = evalFld(x1,y1,z1,ierror)

      ydot(7) = tkernel_3d(s1,ct0,alpha,exit_loop)*ydot(5)  !To normalize

      ydot(6) = tt*ydot(7)

c End program

      end subroutine lsode_3d_T

c lsode_2d_T_3eq
c #########################################################################
      subroutine lsode_2d_T_3eq(neq,t,y,ydot)

c -------------------------------------------------------------------------
c     Driver routine for LSODE integration: ydot = ff(y,t)
c -------------------------------------------------------------------------

      use par_int

      implicit none

c Call variables

      integer :: neq
      real(8) :: t,y(neq),ydot(neq)

c Local variables

      integer :: ierror
      real(8) :: x1,y1,z1,s1,tt,kernel,it_max=0d0
      save :: it_max

c Begin program

c Call LSODE orbit evaluation function (extracting logical coords x1,y1,z1)

      call lsode_orbit_f_3eq(3,t,y(1:3),ydot(1:3),xx1=x1,yy1=y1,zz1=z1)

c Call additional integrators

      tt = evalFld(x1,y1,z1,ierror)
      kernel = tkernel_per(t,ct0,alpha)

      ydot(4) = tt*kernel

      ydot(5) = kernel  !To normalize

c End program

      end subroutine lsode_2d_T_3eq

c lsode_3d_S
c #########################################################################
      subroutine lsode_3d_S(neq,t,y,ydot)

c -------------------------------------------------------------------------
c     Driver routine for LSODE integration: ydot = ff(y,t)
c -------------------------------------------------------------------------

      use par_int

      implicit none

c Call variables

      integer :: neq
      real(8) :: t,y(neq),ydot(neq)

c Local variables

      integer :: ierror
      real(8) :: x1,y1,z1,s1,tt,it_max=0d0
      save :: it_max

c Begin program

c Call LSODE orbit evaluation function (extracting logical coords x1,y1,z1)

      call lsode_orbit_f(5,t,y(1:5),ydot(1:5),xx1=x1,yy1=y1,zz1=z1)

c Call additional integrators

      s1 = abs(y(5))

      tt = evalFld(x1,y1,z1,ierror)

      ydot(7) = skernel_3d(s1,ct0,alpha,exit_loop)*ydot(5)  !To normalize

      ydot(6) = tt*ydot(7)

c End program

      end subroutine lsode_3d_S

c lsode_2d_S_3eq
c #########################################################################
      subroutine lsode_2d_S_3eq(neq,t,y,ydot)

c -------------------------------------------------------------------------
c     Driver routine for LSODE integration: ydot = ff(y,t)
c -------------------------------------------------------------------------

      use par_int

      implicit none

c Call variables

      integer :: neq
      real(8) :: t,y(neq),ydot(neq)

c Local variables

      integer :: ierror
      real(8) :: x1,y1,z1,s1,tt,it_max=0d0
      save :: it_max

c Begin program

c Call LSODE orbit evaluation function (extracting logical coords x1,y1,z1)

      call lsode_orbit_f_3eq(3,t,y(1:3),ydot(1:3),xx1=x1,yy1=y1,zz1=z1)

c Call additional integrators

      tt = evalFld(x1,y1,z1,ierror)
      ydot(5) = skernel_per(t,ct0,alpha)

      ydot(4) = tt*ydot(5)

c End program

      end subroutine lsode_2d_S_3eq

c lsode_3d_bdf2
c #########################################################################
      subroutine lsode_3d_bdf2(neq,t,y,ydot)

c -------------------------------------------------------------------------
c     Driver routine for LSODE integration: ydot = ff(y,t)
c -------------------------------------------------------------------------

      use par_int

      implicit none

c Call variables

      integer :: neq
      real(8) :: t,y(neq),ydot(neq)

c Local variables

      integer :: ierror
      real(8) :: x1,y1,z1,s1,flds(size(sp_coef,4)),B2
      logical :: exl1,exl2,exl3,exl4

c Begin program

      if (neq /= 13) then
        call pstop('lsode_3d_bdf2','Wrong number of equations')
      endif

c Call LSODE orbit evaluation function (extracting logical coords x1,y1,z1)

      call lsode_orbit_f(5,t,y(1:5),ydot(1:5),xx1=x1,yy1=y1,zz1=z1
     .                  ,B2=B2)

      if (.not.bwgt) B2=1d0

c Call additional integrators

      s1 = abs(y(5))

      flds = evalFlds(x1,y1,z1,sp_coef)

      !Tn
      ydot(7) = tkernel_3d(s1,ct0,alpha,exl1)*ydot(5)

      ydot(6) = flds(1)*ydot(7)

      !Tnm
      ydot(9) = tkernel_3d(s1,ct1,alpha,exl2)*ydot(5)

      ydot(8) = flds(2)*ydot(9)

      !Src @ dt
      ydot(11) = skernel_3d(s1,ct0,alpha,exl3)*ydot(5)/B2

      ydot(10) = flds(3)*ydot(11)

      !Src @ dt+dtold
      ydot(13) = skernel_3d(s1,ct1,alpha,exl4)*ydot(5)/B2

      ydot(12) = flds(3)*ydot(13)

c Are we done?

      exit_loop = exl1.and.exl2.and.exl3.and.exl4

c End program

      end subroutine lsode_3d_bdf2

c lsode_3d_bdf2_dT
c #########################################################################
      subroutine lsode_3d_bdf2_dT(neq,t,y,ydot)

c -------------------------------------------------------------------------
c     Driver routine for LSODE integration: ydot = ff(y,t)
c -------------------------------------------------------------------------

      use par_int

      implicit none

c Call variables

      integer :: neq
      real(8) :: t,y(neq),ydot(neq)

c Local variables

      integer :: ierror
      real(8) :: x1,y1,z1,s1,flds(size(sp_coef,4)),B2
      logical :: exl1,exl2,exl3,exl4

c Begin program

      if (neq /= 13) then
        call pstop('lsode_3d_bdf2_dT','Wrong number of equations')
      endif

c Call LSODE orbit evaluation function (extracting logical coords x1,y1,z1)

      call lsode_orbit_f(5,t,y(1:5),ydot(1:5),xx1=x1,yy1=y1,zz1=z1
     .                  ,B2=B2)

      if (.not.bwgt) B2=1d0

c Call additional integrators

      s1 = abs(y(5))

      flds = evalFlds(x1,y1,z1,sp_coef)

      !Src @ dt
      ydot(7) = skernel_3d(s1,ct0,alpha,exl1)*ydot(5)

      ydot(6) = flds(1)*ydot(7)

      !Src * dt+dtold
      ydot(9) = skernel_3d(s1,ct1,alpha,exl2)*ydot(5)

      ydot(8) = flds(2)*ydot(9)

      !Src @ dt
      ydot(11) = skernel_3d(s1,ct0,alpha,exl3)*ydot(5)/B2

      ydot(10) = flds(3)*ydot(11)

      !Src @ dt+dtold
      ydot(13) = skernel_3d(s1,ct1,alpha,exl4)*ydot(5)/B2

      ydot(12) = flds(3)*ydot(13)

c Are we done?

      exit_loop = exl1.and.exl2.and.exl3.and.exl4

c End program

      end subroutine lsode_3d_bdf2_dT

c lsode_3d_bdf2_all
c #########################################################################
      subroutine lsode_3d_bdf2_all(neq,t,y,ydot)

c -------------------------------------------------------------------------
c     Driver routine for LSODE integration: ydot = ff(y,t)
c -------------------------------------------------------------------------

      use par_int

      implicit none

c Call variables

      integer :: neq
      real(8) :: t,y(neq),ydot(neq)

c Local variables

      integer :: ierror
      real(8) :: x1,y1,z1,s1,flds(size(sp_coef,4)),B2,sk,sk2
      logical :: exl1,exl2,exl3,exl4

c Begin program

      if (neq /= 13) then
        call pstop('lsode_3d_bdf2','Wrong number of equations')
      endif

c Call LSODE orbit evaluation function (extracting logical coords x1,y1,z1)

      call lsode_orbit_f(5,t,y(1:5),ydot(1:5),xx1=x1,yy1=y1,zz1=z1
     .                  ,B2=B2)

      if (.not.bwgt) B2=1d0

c Call additional integrators

      s1 = abs(y(5))

      flds = evalFlds(x1,y1,z1,sp_coef)

      !Tn
      ydot(6) = flds(1)*tkernel_3d(s1,ct0,alpha,exl1)

      !Tnm
      ydot(7) = flds(2)*tkernel_3d(s1,ct1,alpha,exl2)

      !Src @ dt
      sk = skernel_3d(s1,ct0,alpha,exl3)
      ydot(8) = flds(4)*sk/B2

      !Src @ dt+dtold
      sk2 = skernel_3d(s1,ct1,alpha,exl4)
      ydot(9) = flds(4)*sk2/B2

      !DT @ dt
      ydot(10) = flds(3)*sk

      !DT @ dt+dtold
      ydot(11) = flds(3)*sk2

      !DT/B^2 @dt
      ydot(12) = flds(3)*sk/B2

      !Src @ dt+dtold
      ydot(13) = flds(3)*sk2/B2

c Are we done?

      exit_loop = exl1.and.exl2.and.exl3.and.exl4

c End program

      end subroutine lsode_3d_bdf2_all

c lsode_2d_bdf2_3eq
c #########################################################################
      subroutine lsode_2d_bdf2_3eq(neq,t,y,ydot)

c -------------------------------------------------------------------------
c     Driver routine for LSODE integration: ydot = ff(y,t)
c -------------------------------------------------------------------------

      use par_int

      implicit none

c Call variables

      integer :: neq
      real(8) :: t,y(neq),ydot(neq)

c Local variables

      integer :: ierror
      real(8) :: x1,y1,z1,s1,it_max=0d0,flds(size(sp_coef,4)),B2
      save :: it_max

c Begin program

      if (neq /= 11) then
        call pstop('lsode_2d_bdf2_3eq','Wrong number of equations')
      endif

c Call LSODE orbit evaluation function (extracting logical coords x1,y1,z1)

      call lsode_orbit_f_3eq(3,t,y(1:3),ydot(1:3),xx1=x1,yy1=y1,zz1=z1
     .                      ,B2=B2)

      if (.not.bwgt) B2=1d0

c Call additional integrators

      flds = evalFlds(x1,y1,z1,sp_coef)

      !Tn
      ydot(5) = tkernel_per(t,ct0,alpha)

      ydot(4) = flds(1)*ydot(5)

      !Tnm
      ydot(7) = tkernel_per(t,ct1,alpha)

      ydot(6) = flds(2)*ydot(7)

      !Src @ dt
      ydot(9) = skernel_per(t,ct0,alpha)/B2

      ydot(8) = flds(3)*ydot(9)

      !Src @ dt+dtold
      ydot(11) = skernel_per(t,ct1,alpha)/B2

      ydot(10) = flds(3)*ydot(11)

c End program

      end subroutine lsode_2d_bdf2_3eq

c lsode_2d_bdf2_3eq_dT
c #########################################################################
      subroutine lsode_2d_bdf2_3eq_dT(neq,t,y,ydot)

c -------------------------------------------------------------------------
c     Driver routine for LSODE integration: ydot = ff(y,t)
c -------------------------------------------------------------------------

      use par_int

      implicit none

c Call variables

      integer :: neq
      real(8) :: t,y(neq),ydot(neq)

c Local variables

      integer :: ierror
      real(8) :: x1,y1,z1,s1,it_max=0d0,flds(size(sp_coef,4)),B2
      save :: it_max

c Begin program

      if (neq /= 11) then
        call pstop('lsode_2d_bdf2_3eq','Wrong number of equations')
      endif

c Call LSODE orbit evaluation function (extracting logical coords x1,y1,z1)

      call lsode_orbit_f_3eq(3,t,y(1:3),ydot(1:3),xx1=x1,yy1=y1,zz1=z1
     .                      ,B2=B2)

      if (.not.bwgt) B2=1d0

c Call additional integrators

      flds = evalFlds(x1,y1,z1,sp_coef)

      !Src @ dt
      ydot(5) = skernel_per(t,ct0,alpha)

      ydot(4) = flds(1)*ydot(5)

      !Src @ dt+dtold
      ydot(7) = skernel_per(t,ct1,alpha)

      ydot(6) = flds(2)*ydot(7)

      !Src @ dt
      ydot(9) = skernel_per(t,ct0,alpha)/B2

      ydot(8) = flds(3)*ydot(9)

      !Src @ dt+dtold
      ydot(11) = skernel_per(t,ct1,alpha)/B2

      ydot(10) = flds(3)*ydot(11)

c End program

      end subroutine lsode_2d_bdf2_3eq_dT

c lsode_2d_bdf2_3eq_all
c #########################################################################
      subroutine lsode_2d_bdf2_3eq_all(neq,t,y,ydot)

c -------------------------------------------------------------------------
c     Driver routine for LSODE integration: ydot = ff(y,t)
c -------------------------------------------------------------------------

      use par_int

      implicit none

c Call variables

      integer :: neq
      real(8) :: t,y(neq),ydot(neq)

c Local variables

      integer :: ierror
      real(8) :: x1,y1,z1,s1,it_max=0d0,flds(size(sp_coef,4)),B2
     .          ,sk,sk2
      save :: it_max

c Begin program

      if (neq /= 11) then
        call pstop('lsode_2d_bdf2_3eq','Wrong number of equations')
      endif

c Call LSODE orbit evaluation function (extracting logical coords x1,y1,z1)

      call lsode_orbit_f_3eq(3,t,y(1:3),ydot(1:3),xx1=x1,yy1=y1,zz1=z1
     .                      ,B2=B2)

      if (.not.bwgt) B2=1d0

c Call additional integrators

      flds = evalFlds(x1,y1,z1,sp_coef)

      !Tn
      ydot(4) = flds(1)*tkernel_per(t,ct0,alpha)

      !Tnm
      ydot(5) = flds(2)*tkernel_per(t,ct1,alpha)

      !Src @ dt
      sk = skernel_per(t,ct0,alpha)
      ydot(6) = flds(4)*sk/B2

      !Src @ dt+dtold
      sk2 = skernel_per(t,ct1,alpha)
      ydot(7) = flds(4)*sk2/B2

      !DT @ dt
      ydot(8) = flds(3)*sk

      !DT @ dt+dtold
      ydot(9) = flds(3)*sk2

      !DT/B^2 @dt
      ydot(10) = flds(3)*sk/B2

      !DTb^2 @ dt+dtold
      ydot(11) = flds(3)*sk2/B2

cc      !Beta
cc      ydot(12) = sk/B2

c End program

      end subroutine lsode_2d_bdf2_3eq_all

c lsode_2d_cn_3eq
c #########################################################################
      subroutine lsode_2d_cn_3eq(neq,t,y,ydot)

c -------------------------------------------------------------------------
c     Driver routine for LSODE integration: ydot = ff(y,t)
c -------------------------------------------------------------------------

      use par_int

      implicit none

c Call variables

      integer :: neq
      real(8) :: t,y(neq),ydot(neq)

c Local variables

      integer :: ierror
      real(8) :: x1,y1,z1,s1,it_max=0d0,flds(size(sp_coef,4)),B2
      save :: it_max

c Begin program

      if (neq /= 7) then
        call pstop('lsode_2d_cn_3eq'
     .            ,'Wrong number of equations')
      endif

c Call LSODE orbit evaluation function (extracting logical coords x1,y1,z1)

      call lsode_orbit_f_3eq(3,t,y(1:3),ydot(1:3),xx1=x1,yy1=y1,zz1=z1
     .                      ,B2=B2)

      if (.not.bwgt) B2=1d0

c Call additional integrators

      flds = evalFlds(x1,y1,z1,sp_coef)

      !Tn
      ydot(5) = tkernel_per(t,ct0,alpha)

      ydot(4) = flds(1)*ydot(5)

      !Src @ dt
      ydot(7) = skernel_per(t,ct0,alpha)/B2

      ydot(6) = flds(3)*ydot(7)

c End program

      end subroutine lsode_2d_cn_3eq

c lsode_2d_beta
c #########################################################################
      subroutine lsode_2d_beta(neq,t,y,ydot)

c -------------------------------------------------------------------------
c     Driver routine for LSODE integration: ydot = ff(y,t)
c -------------------------------------------------------------------------

      use par_int

      implicit none

c Call variables

      integer :: neq
      real(8) :: t,y(neq),ydot(neq)

c Local variables

      integer :: ierror
      real(8) :: B2

c Begin program

      if (neq /= 4) then
        call pstop('lsode_2d_beta','Wrong number of equations')
      endif

c Call LSODE orbit evaluation function (extracting logical coords x1,y1,z1)

      call lsode_orbit_f_3eq(3,t,y(1:3),ydot(1:3),B2=B2)

c Call additional integrators

      ydot(4) = skernel_per(t,ct0,alpha)/B2
cc      ydot(4) = skernel_per(t,ct0,alpha)/sqrt(B2)
cc      ydot(4) = skernel_per(t,ct0,alpha)  !Normalization test

c End program

      end subroutine lsode_2d_beta

c lsode_3d_beta
c #########################################################################
      subroutine lsode_3d_beta(neq,t,y,ydot)

c -------------------------------------------------------------------------
c     Driver routine for LSODE integration: ydot = ff(y,t)
c -------------------------------------------------------------------------

      use par_int

      implicit none

c Call variables

      integer :: neq
      real(8) :: t,y(neq),ydot(neq)

c Local variables

      integer :: ierror
      real(8) :: B2

c Begin program

      if (neq /= 4) then
        call pstop('lsode_3d_beta','Wrong number of equations')
      endif

c Call LSODE orbit evaluation function (extracting logical coords x1,y1,z1)

      call lsode_orbit_f_3eq(3,t,y(1:3),ydot(1:3),B2=B2)

c Call additional integrators

      ydot(4) = skernel_3d(t,ct0,alpha,exit_loop)/B2
cc      ydot(4) = skernel_per(t,ct0)/sqrt(B2)
cc      ydot(4) = skernel_per(t,ct0)  !Normalization test

c End program

      end subroutine lsode_3d_beta
