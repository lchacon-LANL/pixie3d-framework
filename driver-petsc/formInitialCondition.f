c formInitialCondition
c######################################################################
      subroutine formInitialCondition(array,imin,imax,jmin,jmax
     .                                ,kmin,kmax)

c----------------------------------------------------------------------
c     Initializes MG and creates grid
c----------------------------------------------------------------------

      use parameters

      use grid

      use variables

      use timeStepping

      use newtongm

      use constants

      use iosetup

      use icond

      implicit none

c Call variables

      integer(4)      ::imin,imax,jmin,jmax,kmin,kmax

      type(petsc_var) ::array(imin:imax,jmin:jmax,kmin:kmax)

c Local variables

      type(petsc_array) :: petscarray

c Begin program

c Set initial condition

      call setInitialCondition(u_n,u_np)

c Initialize record file

      call initializeRecordFile

c Check time limits

      if (tmax.gt.0d0.and.(tmax-time).le.0d0) then
        write(*,*)
        write(*,*) 'Tmax is less or equal than restarting time'
        write(*,*) 'Aborting'
        stop
      endif

c Initialize counters

cc      gmres_tot = 0
cc      newt_tot  = 0
cc      wh_tot    = 0
cc      ierr      = 0

      itime     = inewtime-1

      nrst      = 0
      tmrst     = 0d0

      dtexp     = 0d0

c Transfer to Petsc format

      petscarray = u_np

      array = petscarray%array(imin:imax,jmin:jmax,kmin:kmax)

      call deallocatePetscType(petscarray)
      call deallocateDerivedType(u_np)

c End program

      end subroutine formInitialCondition


