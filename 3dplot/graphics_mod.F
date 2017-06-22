
c module graphics
c ######################################################################
      module graphics

        use grid

        use xdraw_graphics

        use HDF5_graphics

        use ADIOS_graphics

      contains

c     initializeGraphics
c     ###############################################################
      subroutine initializeGraphics
c     ---------------------------------------------------------------
c     Set graphics files and dumping intervals
c     ---------------------------------------------------------------

        implicit none

c     Call variables

c     Local variables

        integer :: igx,igy,igz,bcs(6)

        integer :: nx,ny,nz,isx,isy,isz

        integer :: i,j,k,nplots
        
        character*(30) :: prof_label

c     Begin program

        igx = 1 ; igy = 1 ; igz = 1

        bcs = bcond

        nx = gv%gparams%nxv(igx)
        ny = gv%gparams%nyv(igy)
        nz = gv%gparams%nzv(igz)

        isx = gv%gparams%istartx(igx)
        isy = gv%gparams%istartx(igy)
        isz = gv%gparams%istartx(igz)

        iggx = igx
        iggy = igy
        iggz = igz

c     Set graphics plotting range

cc        call setGraphicsRange(nx,ny,nz,bcs)

c     Find grid in logical space

        allocate(xl(iming:imaxg),yl(jming:jmaxg),zl(kming:kmaxg)
     .          ,xi(iming:imaxg,jming:jmaxg,kming:kmaxg,3))

        xl(:) = gv%gparams%xx(iming+isx:imaxg+isx)
        yl(:) = gv%gparams%yy(jming+isy:jmaxg+isy)
        zl(:) = gv%gparams%zz(kming+isz:kmaxg+isz)

        do k=kming,kmaxg
          do j=jming,jmaxg
            do i=iming,imaxg
              xi(i,j,k,1) = xl(i)
              xi(i,j,k,2) = yl(j)
              xi(i,j,k,3) = zl(k)
            enddo
          enddo
        enddo

c     Define graphics i/o and initialize graph arrays

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
        if (np == 1 .and. xdraw_plot) call initializeXDRAW
#else
        if (xdraw_plot) call initializeXDRAW
#endif

c     Initialize HDF5

#if defined(hdf5)
        if (hdf_plot) call initializeHDF5
#endif

c     Initialize ADIOS

#if defined(adios)
        if (adios_plot) call initializeADIOSPlot
#endif

ccc     Define application I/O (external)
cc
cc      call defineGraphicsIO

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

        deallocate(xl,yl,zl,xi)

c     Finalize XDRAW

#if defined(petsc)
        if (np == 1.and.xdraw_plot) call finalizeXDRAW
#else
        if (xdraw_plot) call finalizeXDRAW
#endif

c     Finalize HDF5

c     Deallocate memory

        call deallocateGraphicsVariables

c     End program

      end subroutine finalizeGraphics

c     setGraphicsRange
c     #################################################################
      subroutine setGraphicsRange

c     ---------------------------------------------------------------
c     Find graphics plotting range according to boundary conditions
c     ---------------------------------------------------------------

        implicit none

c     Call variables

c     Local variables
 
        integer :: igx,igy,igz,nx,ny,nz

c     Begin program

        igx = 1 ; igy = 1 ; igz = 1

        nx = gv%gparams%nxv(igx)
        ny = gv%gparams%nyv(igy)
        nz = gv%gparams%nzv(igz)

        if (bcond(1) == PER) then
cc          iming = 1  !This is when PER domains start at angle 0
cc          imaxg = nx+1
          if (nx == 1) then
            iming = 1
            imaxg = nx+1
          else
            iming = 0
            imaxg = nx
          endif
        elseif (bcSP()) then
          iming = 0
          imaxg = nx+1
        else
#if defined(PLOT_STG)
          iming = 0
          imaxg = nx
#else
          iming = 0
          imaxg = nx+1
#endif
        endif

        if (bcond(3) == PER) then
cc          jming = 1   !This is when PER domains start at angle 0
cc          jmaxg = ny+1
          if (ny == 1) then
            jming = 1
            jmaxg = ny+1
          else
            jming = 0
            jmaxg = ny
          endif
        else
#if defined(PLOT_STG)
          jming = 0
          jmaxg = ny
#else
          jming = 0
          jmaxg = ny+1
#endif
        endif

        if (bcond(5) == PER) then
cc          kming = 1   !This is when PER domains start at angle 0
cc          kmaxg = nz+1
          if (nz == 1) then
            kming = 1
            kmaxg = nz+1
          else
            kming = 0
            kmaxg = nz
          endif
        else
#if defined(PLOT_STG)
          kming = 0
          kmaxg = nz
#else
          kming = 0
          kmaxg = nz+1
#endif
        endif

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

        if (np == 1 .and. xdraw_plot) call dumpXDRAWplots

#if defined(hdf5)
        if (hdf_plot) call writeHDF5file
#endif
#if defined(adios)
        if (adios_plot) call writeADIOSPlotfile
#endif

c     End program

      end subroutine dumpTimeStepPlots

c     allocateGraphicsVariables
c     ##################################################################
      subroutine allocateGraphicsVariables

        implicit none

c     Call variables

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
