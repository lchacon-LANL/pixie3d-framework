c correctTimeStepRoutine
c ######################################################################
      subroutine correctTimeStepRoutine(array,imin,imax,jmin,jmax,kmin,
     .                                  kmax,gmits,nwits)

      use variables

      use timeStepping

      use iosetup

      use counters

      implicit none

c Call variables

      integer(4)      :: imin,imax,jmin,jmax,kmin,kmax,nwits,gmits

      type(petsc_var) :: array(imin:imax,jmin:jmax,kmin:kmax)

c Local variables

      type(petsc_array) :: petscarray

      integer(4)        :: ierr,i,j,k,ieq

      real(8) :: mag

c Begin program

      ierr = 0

      call allocatePetscType(petscarray)

      petscarray%array(imin:imax,jmin:jmax,kmin:kmax)
     .         = array(imin:imax,jmin:jmax,kmin:kmax)

c Map old solution

      u_n = petscarray

c diag ******
cc      u_n = u_n - u_np
cc      write (*,*) imin,imax,jmin,jmax,kmin,kmax
cc
cc      do ieq=1,neqd
cc        mag = sqrt(sum(u_n%array_var(ieq)
cc     .                 %array(imin:imax,jmin:jmax,1)**2))
cc        write (*,*) mag
cc      enddo

cc      call writeDerivedType(u_n,6,.true.)

cc      open(unit=110,file='debug.bin',form='unformatted'
cc     .       ,status='replace')
cc      do ieq=1,neqd
cc        call contour(u_n%array_var(ieq)%array(imin:imax,jmin:jmax,1)
cc     .              ,imax-imin+1,jmax-jmin+1,0d0,1d0,0d0,1d0,ieq-1,110)
cc      enddo
cc      close(110)

cc      open(unit=110,file='debug.bin',form='unformatted'
cc     .       ,status='replace')
cc      do ieq=1,neqd
cc        call contour(u_n%array_var(ieq)%array(0:nxd+1,0:nyd+1,0)
cc     .              ,nxd+2,nyd+2,0d0,1d0,0d0,1d0,ieq-1,110)
cc      enddo
cc      close(110)
cc
cc      stop
c diag ******

c Find new time step

      call correctTimeStep(u_n,itime+1,ierr)

c Update counters (only if timeStep is successful)

      itime  = itime + 1

      time   = time  + dt
      tmrst  = tmrst + dt
      nrst   = nrst  + 1

c Deallocate memory

      call deallocatePetscType(petscarray)

      end subroutine correctTimeStepRoutine

c correctTimeStep
c####################################################################
      subroutine correctTimeStep(varray,itm,ierr)

c--------------------------------------------------------------------
c     Correct time step
c--------------------------------------------------------------------

      use timeStepping

      use iosetup

      use variables

      use parameters

      implicit none

c Call variables

      integer(4)       ::  ierr,itm

      type (var_array) :: varray

c Local variables

c Begin program

      call calculate_gammat

      call calculate_dt

      if (sm_flag.eq.1) call calculate_cnfactor

      alpha = 1. - cnfactor

c End program

      contains

c     calculate_gammat
c     #######################################################################
      subroutine calculate_gammat

        use generalOperators

c     Calculation of local growth rate for CN

        integer (4) :: ieq

        real(8)     :: dmag1,dmag2,dpert,mag(neqd)

        real(8)     :: array(0:nxdp,0:nydp,0:nzdp)

c     Begin program

        do ieq=1,neqd

          array = (varray%array_var(ieq)%array
     .            -u_0   %array_var(ieq)%array )
          array = array*array

          dpert = integral(nxd,nyd,nzd,array,1,1,1,.true.)

          dpert = sqrt(dpert)

          array = (u_n%array_var(ieq)%array
     .            -u_0%array_var(ieq)%array )
          array = array*array

          dmag1 = integral(nxd,nyd,nzd,array,1,1,1,.true.)

          array = (varray%array_var(ieq)%array
     .            +u_n   %array_var(ieq)%array
     .         -2.*u_0   %array_var(ieq)%array )
          array = array*array

          dmag2 = integral(nxd,nyd,nzd,array,1,1,1,.true.)

          if (dpert /= 0d0.and.dmag2 /= 0d0) then
            mag(ieq) = .5*dt*sqrt(dmag2)/(dpert-sqrt(dmag1))
          else
            mag(ieq) = 1e30
          endif

        enddo

        gammat = 1./minval(abs(mag))

      end subroutine calculate_gammat

c     calculate_cnfactor
c     #######################################################################
      subroutine calculate_cnfactor

        if (itm.eq.1) then
          cnfactor = .5
        elseif (itm.le.sm_pass+1) then
          cnfactor = .3
        else
cc          write (*,*) 'gamma ',gammat,'dt ',dt
          cnfactor = max(.5 - gammat/12.*dt,0d0)
        endif

      end subroutine calculate_cnfactor

c     calculate_dt
c     #######################################################################
      subroutine calculate_dt


        if (timecorr .or. cnfactor == 1d0) then
          if (ierr.eq.0 .and. (itm.eq.1 .or. cnfactor == 1d0)) then
            call findExplicitDt
          else
            call adapt_dt(dtbase)
          endif
        else
          dt = dtbase
        endif

      end subroutine calculate_dt

c     adapt_dt
c     #######################################################################
      subroutine adapt_dt(dtbase)

        real(8) ::    dtbase
        real(8) ::    coef1,coef2

        coef1 = 0.8             !Time subcycling coefficient
        coef2 = 1.05            !Time recovery   coefficient

        if (ierr.gt.0) then
          dt = dt*coef1
          if (dt < 1d-3*dtbase) then
            write (*,*) 'Time step too small'
            write (*,*) 'Aborting...'
            stop
          endif
          if (ierr.eq.2) write (*,240)
          write (*,400) dt
        else
          if (itm.le.sm_pass+1) then
            dt = dtbase/2.
          elseif (itm.eq.(sm_pass+2)) then
            dt = dtbase
          else
            dt = min(dtbase,dt*coef2)
          endif
        endif

 240    format ('    Too many Newton iterations')
 400    format ('    Subcycling time step... New time step:',f7.4)

      end subroutine adapt_dt

      end subroutine correctTimeStep

c     contour
c     #####################################################################
      subroutine contour(arr,nx,ny,xmin,xmax,ymin,ymax,iopt,nunit)
      implicit none               !For safe fortran
c     ---------------------------------------------------------------------
c     Contours arr in xdraw format
c     Notes:
c      put the next 2 lines in main
c      open(unit=nunit,file='contour.bin',form='unformatted') before
c      close(unit=nunit) after
c     ---------------------------------------------------------------------

c     Call variables

      integer(4) :: nx,ny,iopt,nunit
      real(8)    :: arr(nx,ny),xmin,xmax,ymin,ymax

c     Local variables

      integer(4) :: i,j

c     Begin program

      if(iopt.eq.0) then
        write(nunit) nx-1,ny-1,0
        write(nunit) real(xmin,4),real(xmax,4)
     .              ,real(ymin,4),real(ymax,4) 
      endif
      write(nunit) ((real(arr(i,j),4),i=1,nx),j=1,ny)
cc      write (*,*) ((real(arr(i,j),4),i=1,nx),j=1,ny)
cc      write (*,*)

c     End program

      end subroutine contour
