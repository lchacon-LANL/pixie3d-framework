
c module graphics
c ######################################################################
      module graphics

        use grid

        use xdraw_graphics

        use HDF5_graphics

      contains

c     initializeGraphics
c     ###############################################################
      subroutine initializeGraphics (igx,igy,igz,bcs)
c     ---------------------------------------------------------------
c     Set graphics files and dumping intervals
c     ---------------------------------------------------------------

        implicit none

c     Call variables

        integer    :: igx,igy,igz,bcs(6)

c     Local variables

        integer    :: nx,ny,nz,isx,isy,isz

        integer    :: i,nplots
        
        character*(30) :: prof_label

c     Begin program

        nx = grid_params%nxv(igx)
        ny = grid_params%nyv(igy)
        nz = grid_params%nzv(igz)

        isx = grid_params%istartx(igx)
        isy = grid_params%istartx(igy)
        isz = grid_params%istartx(igz)

        iggx = igx
        iggy = igy
        iggz = igz

c     Set graphics plotting range

cc        call setGraphicsRange(nx,ny,nz,bcs)

c     Find grid in logical space

        allocate(xl(iming:imaxg),yl(jming:jmaxg),zl(kming:kmaxg))

        xl(:) = grid_params%xx(iming+isx:imaxg+isx)
        yl(:) = grid_params%yy(jming+isy:jmaxg+isy)
        zl(:) = grid_params%zz(kming+isz:kmaxg+isz)

c     Define graphics i/o and initialize graph arrays

cc        call defineGraphics

        call setGraphicsGroups             !External

        call allocateGraphicsVariables

        call defineGraphicsGroups          !External

c     Determine actual number of 3D plots

        do igroup=1,ngroups
          nqty(igroup) = 0
          do i = 1,size(graph(igroup)%array_graph)
            if (len(trim(graph(igroup)%array_graph(i)%descr)) == 0) exit
            nqty(igroup) = nqty(igroup) + 1
          enddo
        enddo

c     Initialize XDRAW

#if defined(petsc)
        if (np == 1) call initializeXDRAW
#else
        call initializeXDRAW
#endif

c     Initialize HDF5

#if defined(hdf5)
        if (hdf_plot) call initializeHDF5
#endif

c     Define application I/O (external)

      call defineGraphicsIO

c     End program

      end subroutine initializeGraphics

c     finalizeGraphics
c     ###############################################################
      subroutine finalizeGraphics 

c     ---------------------------------------------------------------
c     Close graphics files
c     ---------------------------------------------------------------

        implicit none

c     Call variables

c     Local variables

c     Begin program

c     Finalize XDRAW

#if defined(petsc)
        if (np == 1) call finalizeXDRAW
#else
        call finalizeXDRAW
#endif

c     Finalize HDF5

c     Deallocate memory

        call deallocateGraphicsVariables

c     End program

      end subroutine finalizeGraphics

c     setGraphicsRange
c     #################################################################
      subroutine setGraphicsRange(nx,ny,nz,bcond)

c     ---------------------------------------------------------------
c     Find graphics plotting range according to boundary conditions
c     ---------------------------------------------------------------

        implicit none

c     Call variables

        integer    :: nx,ny,nz,bcond(6)

c     Local variables
 
c     Begin program

        if (bcond(1) == PER) then
          iming = 1
          imaxg = nx+1
        elseif (bcSP()) then
cc          iming = 1
          iming = 0
          imaxg = nx+1
        else
          iming = 0
          imaxg = nx+1
        endif

        if (bcond(3) == PER) then
          jming = 1
cc          jming = 0
          jmaxg = ny+1
        else
          jming = 0
          jmaxg = ny+1
        endif

        if (bcond(5) == PER) then
          kming = 0
          kmaxg = nz+1
        else
          kming = 0
          kmaxg = nz+1
        endif

cc          iming = 0
cc          imaxg = nx+1
cc          jming = 0
cc          jmaxg = ny+1
cc          kming = 0
cc          kmaxg = nz+1

cc          iming = ilo
cc          imaxg = ihi
cc          jming = jlo
cc          jmaxg = jhi
cc          kming = klo
cc          kmaxg = khi

c     End program

      end subroutine setGraphicsRange

c     dumpTimeStepPlots
c     #################################################################
      subroutine dumpTimeStepPlots

c     ---------------------------------------------------------------
c     Dumps time plots
c     ---------------------------------------------------------------

        use variables

        implicit none

c     Call variables

c     Local variables

c     Begin program

        if (np == 1) call dumpXDRAWplots

#if defined(hdf5)
        if (hdf_plot) call writeHDF5file
#endif

c     End program

      end subroutine dumpTimeStepPlots

c     allocateGraphicsVariables
c     ##################################################################
      subroutine allocateGraphicsVariables

        implicit none

c     Call variables

cc        integer    :: ngroups

c     Begin program

c     Allocate graph groups

        allocate(graph(ngroups),nqty(ngroups))

c     Allocate XDRAW quantities

        call allocateXDRAWvariables

      end subroutine allocateGraphicsVariables

c     deallocateGraphicsVariables
c     ##################################################################
      subroutine deallocateGraphicsVariables

        implicit none

c     Call variables

c     Local variables

        integer :: alloc_stat

c     Begin program

        deallocate(graph,nqty)

        call deallocateXDRAWvariables

      end subroutine deallocateGraphicsVariables

      end module graphics
