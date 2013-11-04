c module prj_int
c ######################################################################
      module prj_int

        use ts_base

        logical :: prj_int_ts=.false.

        integer :: prj_int_k=3     !Number of inner steps
     $            ,prj_int_M=10    !Number of steps for outer projection
        real(8) :: prj_int_cnf=0d0 !Theta parameter for inner solver

        logical,private :: debug=.true.

      contains

c     prj_int_update
c     ######################################################################
      function prj_int_update(vnp,k,M,iout) result(ierr)

c     ----------------------------------------------------------------------
c     Performs time advance of solution from vn to vnp using projective
c     integration.
c     ----------------------------------------------------------------------

      implicit none

c     Call variables

      integer :: k,M,ierr,iout

      type(var_array),pointer :: vnp

c     Local variables

      integer :: it
      real(8) :: x(ntotd),h,alpha1,alpha2,cnf_save
      type(var_array),pointer :: vk=>null(),vkp=>null()

c     Begin program

      debug = debug.and.(my_rank == 0).and.(iout > 0)

      cnf_save = cnfactor
      cnfactor = prj_int_cnf

      !Define substep
      h = dt/(k+M+1)

      !Save initial condition
      call equateDerivedType(vk,vnp)

      !Integrate k substeps
      if (debug) write (*,*) 'PROJ. INT.: Integrating',k,' substeps...'
      do it=1,k
        ierr = updateSolution(vk,iout,h)
      enddo

      !Integrate final step to k+1
      if (debug)
     $     write (*,*) 'PROJ. INT.: Integrating projection substep...'
      call equateDerivedType(vkp,vk)
      ierr = updateSolution(vkp,iout,h)

      !Project vk,vkp -> vnp
      if (debug) write (*,*) 'PROJ. INT.: Projecting forward...'
      alpha1 = dble(M+1)
      alpha2 =-dble(M)
      call AXPYDerivedType(alpha1,vkp,alpha2,vk,vnp)

c     End program

      call deallocateDerivedType(vk)
      call deallocateDerivedType(vkp)

      cnfactor = cnf_save

      end function prj_int_update

c     read_PrjInt_Input
c     #################################################################
      subroutine read_PrjInt_Input(file)

c     -----------------------------------------------------------------
c     Initializes grid parameters
c     -----------------------------------------------------------------

      implicit none

c     Call variables

      character(*) :: file

c     Local variables

      integer :: ioerr,uinput=1234

c     Namelists

      !Temporal discretization
      namelist /prj_int/ prj_int_k,prj_int_M,prj_int_cnf,prj_int_ts

c     Begin program

c     Set defaults (see module head)

c     Read initialization parameters

      open(unit=uinput,file=trim(file),status='old')

      read(uinput,prj_int,iostat=ioerr)

      if (ioerr/=0)
     $     call pstop('readGridInput','Problem reading "prj_int"')

      close(unit=uinput)

c     Consistency check

c     Write initialization parameters to standard output

      if (my_rank == 0) then
        write (*,*) '---------------------'
        write (*,*) 'PRJ_INT configuration'
        write (*,*) '---------------------'
        write (*,prj_int)
      endif

c     End program

      end subroutine read_PrjInt_Input

      end module prj_int
