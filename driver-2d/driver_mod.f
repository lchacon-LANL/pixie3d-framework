c module constants
c ######################################################################
      module constants

        double precision :: pi
        double precision, allocatable, dimension(:,:) :: zeros,ones

      end module constants

c module parameters
c #####################################################################
      module parameters

        integer*4 neqd,nxd,nyd,nxdp,nydp,ntotd,ndiagdp,ntotdp,ntotd2p

      contains

c     setVectorDimensions
c     #################################################################
      subroutine setVectorDimensions

        nxdp = nxd+1
        nydp = nyd+1

        ntotdp  = nxd*nyd
        ntotd2p = 2*ntotdp
        ntotd   = neqd*ntotdp

        ndiagdp=5                !Stencil size

      end subroutine setVectorDimensions

      end module parameters

c module vectorToArrayXfer
c #####################################################################
      module vectorToArrayXfer

      contains

c     mapVectorToArray
c     ###############################################################
      subroutine mapVectorToArray (nx,ny,vector,array,array0,bcs)

c     ---------------------------------------------------------------
c     Maps vector to array filling ghost nodes.
c     ---------------------------------------------------------------

      implicit none

c     Call variables

      integer         :: nx,ny,i,j,bcs(4)
      double precision:: array(0:nx+1,0:ny+1),array0(0:nx+1,0:ny+1)
     .                  ,vector(nx*ny)

c     Local variables

c     Begin program

      do i = 1,nx
        do j = 1,ny
          array(i,j) = vector(i + nx*(j-1))
        enddo
      enddo

      call setBC(array,array0,nx,ny,bcs)

c     End

      end subroutine

c     setBC
c     ###############################################################
      subroutine setBC(array,array0,nx,ny,bcs)

c     ---------------------------------------------------------------
c     Sets adequate boundary conditions on array.
c
c     On input:
c       * array  -> real array with ghost-nodes
c       * array0 -> real array containing Dirichlet BC's
c       * nx,ny  -> domain size
c       * bcs    -> real array of size 4 containing BC setup:
c           + bcs(1) ---> bottom
c           + bcs(2) ---> top
c           + bcs(3) ---> left
c           + bcs(4) ---> right
c         If bcs = 0  --> periodic BC's
c         If bcs = 1  --> dirichlet BC's
c         If bcs = 2  --> natural (Neumann) BC's
c         If bcs = 3  --> second-order derivative = 0 (linear extrapolation)
c     ---------------------------------------------------------------

      implicit none       !For safe fortran

c     Call variables

      integer*4  nx,ny,bcs(4)
      real*8     array(0:nx+1,0:ny+1),array0(0:nx+1,0:ny+1)

c     Local variables

      integer*4  i,j

c     Begin program

c     Bottom

      if (bcs(1).eq.0) then
        array(1:nx,0) = array(1:nx,ny-1)
      elseif (bcs(1).eq.1) then
        array(1:nx,0) = array0(1:nx,0)
      elseif (bcs(1).eq.2) then
        array(1:nx,0) = array(1:nx,1)
cc        array(1:nx,0) = array(1:nx,2)
      elseif (bcs(1).eq.3) then
        array(1:nx,0) = 2*array(1:nx,1) - array(1:nx,2)
      else
        write (*,*) 'Boundary condition type not implemented at bottom'
      endif

c     Top

      if (bcs(2).eq.0) then
        array(1:nx,ny+1) = array(1:nx,2)
      elseif (bcs(2).eq.1) then
        array(1:nx,ny+1) = array0(1:nx,ny+1)
      elseif (bcs(2).eq.2) then
        array(1:nx,ny+1) = array(1:nx,ny)
cc        array(1:nx,ny+1) = array(1:nx,ny-1)
      elseif (bcs(2).eq.3) then
        array(1:nx,ny+1) = 2*array(1:nx,ny) - array(1:nx,ny-1)
      else
        write (*,*) 'Boundary condition type not implemented at top'
      endif

c     Left

      if (bcs(3).eq.0) then
        array(0,:) = array(nx-1,:)
      elseif (bcs(3).eq.1) then
        array(0,:) = array0(0,:)
      elseif (bcs(3).eq.2) then
        array(0,:) = array(1,:)
cc        array(0,:) = array(2,:)
      elseif (bcs(3).eq.3) then
        array(0,:) = 2*array(1,:) - array(2,:)
      else
        write (*,*) 'Boundary condition type not implemented at left'
      endif

c     Right

      if (bcs(4).eq.0) then
        array(nx+1,:) = array(2,:)
      elseif (bcs(4).eq.1) then
        array(nx+1,:) = array0(nx+1,:)
      elseif (bcs(4).eq.2) then
        array(nx+1,:) = array(nx,:)
cc        array(nx+1,:) = array(nx-1,:)
      elseif (bcs(4).eq.3) then
        array(nx+1,:) = 2*array(nx,:) - array(nx-1,:)
      else
        write (*,*) 'Boundary condition type not implemented at right'
      endif

c     Resolve singularity if no BC is Dirichlet

cc      if (      bcs(1) /= 1 .and. bcs(2) /= 1 
cc     .    .and. bcs(3) /= 1 .and. bcs(4) /= 1 ) then
cc        array(1,1) = 0d0
cc      endif

c     End

      end subroutine

      end module vectorToArrayXfer

c module variable_setup
c #####################################################################
      module variable_setup

        use parameters

        use vectorToArrayXfer

        implicit none

        type :: var_def
          integer :: bconds(4)          !Boundary conditions
          double precision,pointer,dimension(:,:) :: array
cold          double precision,dimension(0:nxdp,0:nydp) :: array
          character(20):: descr
        end type var_def

        type :: var_array
          integer :: nvar               !Number of variables
          type (var_def),pointer,dimension(:) :: array_var
        end type var_array

        type (var_array)         :: utmp

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
          do ieq=1,neqd
            allocate(varray%array_var(ieq)%array(0:nxdp,0:nydp))
            varray%array_var(ieq)%array = 0d0
          enddo
        endif

c     End program

      end subroutine allocateDerivedType

c     deallocateDerivedType
c     #################################################################
      subroutine deallocateDerivedType(varray)

c     Call variables

        type(var_array)  :: varray

c     Local variables

        integer          :: ieq

c     Begin program

        if (associated(varray%array_var)) then
          do ieq=1,neqd
            deallocate(varray%array_var(ieq)%array)
          enddo
          deallocate(varray%array_var)
        endif

c     End program

      end subroutine deallocateDerivedType

c     printDerivedType
c     #################################################################
      subroutine writeDerivedType(varray,unit,format)

c     Call variables

        type(var_array)  :: varray
        integer(4)       :: unit
        logical          :: format

c     Local variables

        integer          :: ieq

c     Begin program

        if (format) then
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
      subroutine equateDerivedType(varray1,varray2)

c     -----------------------------------------------------------------
c     Performs varray2 = varray1
c     -----------------------------------------------------------------

c     Call variables

        type(var_array)  :: varray1,varray2

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

c     varPack
c     #################################################################
      subroutine varPack(nx,ny,array,bcs,desc,veq,varray)

c     Call variables

        integer          :: nx,ny,veq,bcs(4)
        double precision :: array(0:nx+1,0:ny+1)
        character(*)     :: desc

        type (var_array) :: varray

c     Local variables

c     Begin program

        varray%array_var(veq)%bconds = bcs
        varray%array_var(veq)%array  = array
        varray%array_var(veq)%descr  = desc

c     End program

      end subroutine varPack

c     varUnPack
c     #################################################################
      subroutine varUnPack(varray,veq,nx,ny,array,bcs,desc)

c     Call variables

        integer          :: nx,ny,veq,bcs(4)
        double precision :: array(0:nx+1,0:ny+1)
        character(20)    :: desc

        type (var_array) :: varray

c     Local variables

c     Begin program

        bcs   = varray%array_var(veq)%bconds
        array = varray%array_var(veq)%array
        desc  = trim(varray%array_var(veq)%descr)

c     End program

      end subroutine varUnPack

c     mapVectorToStructure
c     #################################################################
      subroutine mapVectorToStructure(x,varray,varray0)

c     -----------------------------------------------------------------
c     Maps Newton vector solution into arrays using varray0 for 
c     boundary conditions.
c     -----------------------------------------------------------------

      implicit none

c     Call variables

      real(8)          :: x(ntotd)

      type (var_array) :: varray,varray0

c     Local variables

      integer          :: i,j,ii,ieq,neq

c     Begin program

c     Initialize varray with template

cc      varray = utmp
      call equateDerivedType(utmp,varray)

c     Unpack vector x

      neq = varray%nvar

      do ieq = 1,neq
        do j = 1,nyd
          do i = 1,nxd
            ii = ieq + neq*(i-1) + neq*nxd*(j-1)
            varray%array_var(ieq)%array(i,j) = x(ii)
          enddo
        enddo
      enddo

c     Impose boundary conditions

      call imposeBC(varray,varray0,nxd,nyd)

c     End program

      end subroutine

c     mapStructureToVector
c     #################################################################
      subroutine mapStructureToVector(varray,x)

c     -----------------------------------------------------------------
c     Maps solution into Newton vector (without ghost cells)
c     -----------------------------------------------------------------

      implicit none

c     Call variables

      double precision :: x(ntotd)

      type (var_array) :: varray

c     Local variables

      integer*4   i,j,ii,ieq,neq

c     Begin program

c     Unpack vector x

      neq = varray%nvar

      do j = 1,nyd
        do i = 1,nxd
          do ieq = 1,neq
            ii = ieq + neq*(i-1) + neq*nxd*(j-1)
            x(ii) = varray%array_var(ieq)%array(i,j)
          enddo
        enddo
      enddo

c     End program

      end subroutine mapStructureToVector

c     imposeBC
c     ###############################################################
      subroutine imposeBC (varray,varray0,nx,ny)
      implicit none
c     ---------------------------------------------------------------
c     Sets adequate boundary conditions on array structure varray.
c     ---------------------------------------------------------------

c     Call variables

      integer          :: nx,ny
      type (var_array) :: varray,varray0

c     Local variables

      integer         :: neq,ieq

c     Begin program

      neq = varray%nvar

      do ieq = 1,neq 
cc        call varUnPack(varray ,ieq,nx,ny,array ,bcs,desc)
cc        call varUnPack(varray0,ieq,nx,ny,array0,bcs,desc)

cc        call setBC(array,array0,nx,ny,bcs)

cc        call varPack(nx,ny,array,bcs,desc,ieq,varray)

        call setBC(varray %array_var(ieq)%array
     .            ,varray0%array_var(ieq)%array
     .            ,nx,ny,varray%array_var(ieq)%bconds)

      enddo

c     End

      end subroutine imposeBC

      end module variable_setup

c module variables
c #####################################################################
      module variables

        use variable_setup

cnew        type (var_array) :: u_np,u_n,u_0
        type (var_array),target :: u_np,u_n,u_0

      contains

c     allocatePhysicsVariables
c     #################################################################
      subroutine allocatePhysicsVariables

c     Call variables

c     Local variables

c     Begin program

        call allocateDerivedType(u_np)
        call allocateDerivedType(u_n )
        call allocateDerivedType(u_0)
        call allocateDerivedType(utmp)

c     End program

      end subroutine allocatePhysicsVariables


      end module variables

c module icond
c ######################################################################
      module icond

        integer :: nh,bcond(4)

        double precision,dimension(10) :: pert

      end module icond

c module resistiveWall
c ######################################################################
      module resistiveWall

        use parameters

        !Resistive wall
        logical          :: reswall
        double precision :: tau_wall
        double precision :: drag,u00,u01,du0,tu0_step
        double precision :: maxstrs,uxk,uxkp

        character(4)     :: init

        logical          :: allmodes

        !Vacuum wall
        integer          :: nhvw
        double precision :: eps,phase,vw,w_vw,gain,shift,gr,gi

      end module resistiveWall

c module timeStepping
c ######################################################################
      module timeStepping

        use parameters

        integer :: numtime,ndstep,nrstep,sm_pass,sm_flag,inewtime

        double precision :: dt,dtbase,cnfactor,cnf_d,alpha,gamma
        double precision :: tmax,dstep,rstep,time,dfreq,dtexp

        double precision,allocatable,dimension(:) :: fold,fsrc,cnf
     .                                              ,one_over_dt

        double precision :: u_max,v_max,bx_max,by_max

        logical          :: timecorr,restart,source

      end module timeStepping

c module discret_params
c ######################################################################
      module discret_params

        integer :: advect

        logical :: conserv

      end module discret_params

c module newtongm
c ######################################################################
      module newtongm

        integer          ::  maxitnwt,maxitgm,iguess,kmax,method,global

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


c module graphics
c ######################################################################
      module graphics

        double precision :: diagnostics(20),profiles(20)

        character*(20)   :: diag_desc(20),prof_desc(20)

        integer          :: sel_diag(9),sel_cont(9),nqty,diag_ivar(20)

        integer          :: imin,imax,jmin,jmax

        type :: graph_var_def
          double precision,pointer,dimension(:,:) :: array
          character(20):: descr
        end type graph_var_def

        type (graph_var_def),dimension(20) :: array_graph

      end module graphics

c module iosetup
c ######################################################################
      module iosetup

        character*(13)::  inputfile,contourfile,lineplotfile
     .                   ,profilefile,debugfile

        integer  ::     ucontour,ulineplot,uprofile,udebug

        integer  ::     ilevel,mfile

        logical  ::     plot,debug

      contains

      subroutine defineFiles

        contourfile  = 'contours.bin'
        ucontour     = 20

        lineplotfile = 'timeplots.bin'
        ulineplot    = 8

        profilefile  = 'profiles.bin'
        uprofile     = 15

        debugfile    = 'debug.bin'
        udebug       = 30

        mfile        = 40

      end subroutine defineFiles

c     openGraphicsFiles
c     ##################################################################
      subroutine openGraphicsFiles(restart)

        logical :: restart

        if (.not.restart) then

          open(unit=ulineplot,file=lineplotfile,form='unformatted'
     .        ,status='replace')
          open(unit=ucontour ,file=contourfile ,form='unformatted'
     .        ,status='replace')
          open(unit=uprofile ,file=profilefile ,form='unformatted'
     .        ,status='replace')

        else

          open(unit=ulineplot,file=lineplotfile,form='unformatted'
     .        ,status='old',position='append')
          open(unit=ucontour ,file=contourfile ,form='unformatted'
     .        ,status='old',position='append')
          open(unit=uprofile ,file=profilefile ,form='unformatted'
     .        ,status='old',position='append')

        endif

        if (debug) open(unit=udebug,file=debugfile,form='unformatted')

      end subroutine openGraphicsFiles

c     closeGraphicsFiles
c     ##################################################################
      subroutine closeGraphicsFiles

        close(unit=ulineplot)
        close(unit=uprofile)
        close(unit=ucontour)
        if (debug) close(unit=udebug)

      end subroutine closeGraphicsFiles

      end module iosetup
