c module constants
c ######################################################################
      module constants

        double precision :: pi
        double precision, allocatable, dimension(:,:,:) :: zeros,ones

      end module constants

c module parameters
c #####################################################################
      module parameters

      integer(4) :: neqd
     .             ,nxd,nyd,nzd
     .             ,nxdp,nydp,nzdp

      integer(4) :: ntotd,ntotdp,ntotd2p,ntimemax

cc      parameter(ntimemax=1000000)

      contains

c     setVectorDimensions
c     #################################################################
      subroutine setVectorDimensions

        nxdp = nxd+1
        nydp = nyd+1
        nzdp = nzd+1

        ntotdp  = nxd*nyd*nzd
        ntotd2p = 2*ntotdp
        ntotd   = neqd*ntotdp

      end subroutine setVectorDimensions

      end module parameters

c module variable_setup
c #####################################################################
      module variable_setup

        use parameters

        implicit none

        type :: var_def
          integer :: bconds(6)          !Boundary conditions
          double precision,pointer,dimension(:,:,:) :: array
          character(20):: descr
        end type var_def

        type :: var_array
          integer :: nvar               !Number of variables
          type (var_def),pointer,dimension(:) :: array_var
        end type var_array

        type (var_array) :: u_0

        INTERFACE ASSIGNMENT (=)
          module procedure equateDerivedType
     .                    ,mapStructureToVector
     .                    ,mapVectorToStructure
        END INTERFACE

        INTERFACE OPERATOR (-)
          module procedure substractDerivedType
        END INTERFACE

      contains

c     allocateDerivedType
c     #################################################################
      subroutine allocateDerivedType(varray)

c     Call variables

        type(var_array)  :: varray

c     Local variables

        integer          :: ieq

c     Begin program

        varray%nvar = neqd

        if (.not.associated(varray%array_var)) then
          allocate(varray%array_var(neqd))
        endif

        if (.not.associated(varray%array_var(1)%array)) then
          do ieq=1,neqd
            allocate(varray%array_var(ieq)%array(0:nxdp,0:nydp,0:nzdp))
cc            varray%array_var(ieq)%array = 0d0
          enddo
        endif

c     End program

      end subroutine allocateDerivedType

c     writeDerivedType
c     #################################################################
      subroutine writeDerivedType(varray,unit,frmt)

c     Call variables

        type(var_array)  :: varray
        integer(4)       :: unit
        logical          :: frmt

c     Local variables

        integer          :: ieq

c     Begin program

        if (frmt) then
          write (unit,*) varray%nvar
          do ieq=1,varray%nvar
            write(unit,*) varray%array_var(ieq)%array
            write(unit,*) varray%array_var(ieq)%bconds
            write(unit,*) varray%array_var(ieq)%descr
          enddo
        else
          write (unit) varray%nvar
          do ieq=1,varray%nvar
            write(unit) varray%array_var(ieq)%array
            write(unit) varray%array_var(ieq)%bconds
            write(unit) varray%array_var(ieq)%descr
          enddo
        endif

c     End program

      end subroutine writeDerivedType

c     readDerivedType
c     #################################################################
      subroutine readDerivedType(varray,unit,format)

c     Call variables

        type(var_array)  :: varray
        integer(4)       :: unit
        logical          :: format

c     Local variables

        integer          :: ieq

c     Begin program

        if (format) then
          read (unit,*) varray%nvar
          do ieq=1,varray%nvar
            read(unit,*) varray%array_var(ieq)%array
            read(unit,*) varray%array_var(ieq)%bconds
            read(unit,*) varray%array_var(ieq)%descr
          enddo
        else
          read (unit) varray%nvar
          do ieq=1,varray%nvar
            read(unit) varray%array_var(ieq)%array
            read(unit) varray%array_var(ieq)%bconds
            read(unit) varray%array_var(ieq)%descr
          enddo
        endif

c     End program

      end subroutine readDerivedType

c     equateDerivedType
c     #################################################################
      subroutine equateDerivedType(varray2,varray1)

c     -----------------------------------------------------------------
c     Performs varray2 = varray1
c     -----------------------------------------------------------------

c     Call variables

        type(var_array),intent(in)  :: varray1
        type(var_array),intent(out) :: varray2

c     Local variables

        integer          :: ieq

c     Begin program

        call allocateDerivedType(varray2)

        varray2%nvar = varray1%nvar
        do ieq=1,varray2%nvar
          varray2%array_var(ieq)%bconds = varray1%array_var(ieq)%bconds
          varray2%array_var(ieq)%array  = varray1%array_var(ieq)%array
          varray2%array_var(ieq)%descr  = varray1%array_var(ieq)%descr
        enddo

c     End program

      end subroutine equateDerivedType

c     initializeDerivedType
c     #################################################################
      subroutine initializeDerivedType(varray)

c     -----------------------------------------------------------------
c     Initializes fields of varray2 using a template (u_0) except for
c     numeric arrays.
c     -----------------------------------------------------------------

c     Call variables

        type(var_array)  :: varray

c     Local variables

        integer          :: ieq

c     Begin program

        call allocateDerivedType(varray)

        varray%nvar = u_0%nvar
        do ieq=1,varray%nvar
          varray%array_var(ieq)%bconds = u_0%array_var(ieq)%bconds
          varray%array_var(ieq)%descr  = u_0%array_var(ieq)%descr
cc          varray%array_var(ieq)%array  = 0d0
        enddo

c     End program

      end subroutine initializeDerivedType

c     substractDerivedType
c     #################################################################
      function substractDerivedType(varray2,varray1) result (varray3)

c     -----------------------------------------------------------------
c     Performs varray3 = varray2 - varray1
c     -----------------------------------------------------------------

c     Call variables

        type(var_array),intent(in)  :: varray1,varray2
        type(var_array)             :: varray3

c     Local variables

        integer          :: ieq

c     Begin program

        call allocateDerivedType(varray3)

        varray3%nvar = varray2%nvar
        do ieq=1,varray3%nvar
          varray3%array_var(ieq)%bconds = varray2%array_var(ieq)%bconds
          varray3%array_var(ieq)%array  = varray2%array_var(ieq)%array
     .                                   -varray1%array_var(ieq)%array
          varray3%array_var(ieq)%descr  = varray2%array_var(ieq)%descr
        enddo

c     End program

      end function substractDerivedType

c     varPack
c     #################################################################
      subroutine varPack(nx,ny,nz,array,bcs,desc,veq,varray)

c     -----------------------------------------------------------------
c     Fills structure fields of varray with pertinent data corresponding
c     to variable veq: boundary conditions, bcs; variable description,
c     desc; numerical data, array).
c     -----------------------------------------------------------------

c     Call variables

        integer          :: nx,ny,nz,veq,bcs(6)
        double precision :: array(0:nx+1,0:ny+1,0:nz+1)
        character(*)     :: desc

        type (var_array) :: varray

c     Local variables

        call allocateDerivedType(varray)

c     Begin program

        varray%array_var(veq)%bconds = bcs
        varray%array_var(veq)%array  = array
        varray%array_var(veq)%descr  = desc

c     End program

      end subroutine varPack

c     varUnPack
c     #################################################################
      subroutine varUnPack(varray,veq,nx,ny,nz,array,bcs,desc)

c     Call variables

        integer          :: nx,ny,nz,veq,bcs(6)
        double precision :: array(0:nx+1,0:ny+1,0:nz+1)
        character(20)    :: desc

        type (var_array) :: varray

c     Local variables

c     Begin program

        bcs   = varray%array_var(veq)%bconds
        array = varray%array_var(veq)%array
        desc  = trim(varray%array_var(veq)%descr)

c     End program

      end subroutine varUnPack

c     mapStructureToVector
c     #################################################################
      subroutine mapStructureToVector(x,varray)

c     -----------------------------------------------------------------
c     Maps solution into Newton vector (without ghost cells)
c     -----------------------------------------------------------------

      implicit none

c     Call variables

      double precision,intent(out) :: x(ntotd)

      type (var_array),intent(in)  :: varray

c     Local variables

      integer*4   i,j,k,ii,ieq,neq

c     Begin program

c     Unpack vector x

      neq = varray%nvar

      do k = 1,nzd
        do j = 1,nyd
          do i = 1,nxd
            do ieq = 1,neq
              ii = ieq + neq*( (i-1) + nxd*(j-1) + nxd*nyd*(k-1) )
              x(ii) = varray%array_var(ieq)%array(i,j,k)
            enddo
          enddo
        enddo
      enddo

c     End program

      end subroutine mapStructureToVector

c     mapVectorToStructure
c     #################################################################
      subroutine mapVectorToStructure(varray,x)

c     -----------------------------------------------------------------
c     Maps Newton vector solution into arrays using varray0 for 
c     boundary conditions.
c     -----------------------------------------------------------------

      implicit none

c     Call variables

      real(8),intent(in)           :: x(ntotd)

      type (var_array),intent(out) :: varray

c     Local variables

      integer          :: i,j,k,ii,ieq,neq

cc      external imposeBoundaryConditions

c     Begin program

c     Initialize varray

      call initializeDerivedType(varray)

c     Unpack vector x

      neq = varray%nvar

      do ieq = 1,neq
        do k = 1,nzd
          do j = 1,nyd
            do i = 1,nxd
              ii = ieq + neq*( (i-1) + nxd*(j-1) + nxd*nyd*(k-1) )
              varray%array_var(ieq)%array(i,j,k) = x(ii)
            enddo
          enddo
        enddo
      enddo

c     Impose boundary conditions (external)

cc      call imposeBoundaryConditions(varray)

c     End program

      end subroutine

      end module variable_setup

c module variables
c #####################################################################
      module variables

        use variable_setup

        type (var_array),target :: u_np,u_n,u_graph

      contains

c     allocateStructures
c     #################################################################
      subroutine allocateStructures

c     Call variables

c     Local variables

c     Begin program

cc        call allocateDerivedType(u_np)
cc        call allocateDerivedType(u_n )
cc        call allocateDerivedType(u_0)
cccc        call allocateDerivedType(utmp)
cc        call allocateDerivedType(u_graph)

c     End program

      end subroutine allocateStructures

      end module variables

c module icond
c ######################################################################
      module icond

        implicit none

        integer :: nh1,nh2,nh3

        logical :: odd,random

        double precision,dimension(10) :: pert

      end module icond

c module timeStepping
c ######################################################################
      module timeStepping

        use parameters

        integer :: numtime,ndstep,nrstep,sm_pass,sm_flag,inewtime

        double precision :: dt,dtbase,cnfactor,cnf_d,alpha,gammat
        double precision :: tmax,dstep,rstep,time,dfreq,dtexp

        double precision,allocatable,dimension(:) :: fold,fsrc

        double precision :: vx_max,vy_max,vz_max,bx_max,by_max,bz_max

        logical          :: timecorr,restart,source

      end module timeStepping

c module newtongm
c ######################################################################
      module newtongm

        integer          ::  maxitnwt,maxitgm,iguess,maxksp
     .                      ,method,global

        double precision ::  tolgm,tolnewt

      end module newtongm

c module precond_setup
c ######################################################################
      module precond_setup

        integer       ::  precpass,nsweep,maxvcyc

        character*(10)::  precon

      end module precond_setup

c module counters
c ######################################################################
      module counters

        integer          :: itnewt,itgmres,itwhistler
        integer          :: gmres_tot,newt_tot,wh_tot

      end module counters

c module iosetup
c ######################################################################
      module iosetup

        integer  ::     ilevel

      end module iosetup
