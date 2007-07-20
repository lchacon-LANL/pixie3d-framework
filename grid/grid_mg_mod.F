c module grid_mg
c######################################################################
      module grid_mg

        use grid_mpi

        use grid_operations

        integer    :: ngrdx,ngrdy,ngrdz,ngrid,igmax

        integer   ,dimension(:),allocatable ::
     .                     istart,ntotv,istartp,ntotvp
     .                    ,istartb,ntotb,nxv,nyv,nzv,nblock
     .                    ,mg_ratio_x,mg_ratio_y,mg_ratio_z

        type(grid_mg_def),pointer :: MGgrid,MGgrid_sv

      contains

c     allocPointers
c     #################################################################
      subroutine allocPointers(neq,fpointers,mg_grid)

c     -----------------------------------------------------------------
c     Initializes pointers for block MG. In call sequence:
c       * neq (in), integer: number of coupled equations to be solved
c       * mg_grid (in/out), type(grid_mg_def),pointer: contains definition of MG grid
c            level structure.
c       * fpointers (out), logical: indicates if pointers were allocated
c            or not.
c     -----------------------------------------------------------------

      implicit none        !For safe fortran

c     Call variables

      integer             :: neq
      logical,intent(OUT) :: fpointers
      type(grid_mg_def),pointer,optional :: mg_grid

c     Local variables

      integer    :: i

c     Begin program

      fpointers = .false.

      if (PRESENT(mg_grid)) then
        MGgrid=mg_grid
      else
cc        if (.not.associated(MGgrid%xx)) MGgrid=grid_params     !Failsafe grid definition
        if (.not.associated(MGgrid))
     .       call equateGridStructure(MGgrid,grid_params) !Failsafe grid definition
      endif

      ngrdx = MGgrid%ngrdx
      ngrdy = MGgrid%ngrdy
      ngrdz = MGgrid%ngrdz
      ngrid = MGgrid%ngrid

c     Check if MG pointers are allocated

      if (.not.allocated(istart)) then

        allocate(istart (ngrid+1),ntotv (ngrid)
     .          ,istartp(ngrid+1),ntotvp(ngrid)
     .          ,istartb(ngrid+1),ntotb (ngrid)
     .          ,nxv(ngrid),nyv(ngrid),nzv(ngrid),nblock(ngrid)
     .          ,mg_ratio_x(ngrid),mg_ratio_y(ngrid),mg_ratio_z(ngrid))

        mg_ratio_x = MGgrid%mg_ratio_x
        mg_ratio_y = MGgrid%mg_ratio_y
        mg_ratio_z = MGgrid%mg_ratio_z

        nxv = MGgrid%nxv
        nyv = MGgrid%nyv
        nzv = MGgrid%nzv

        istartp(1) = 1
        istart (1) = 1
        istartb(1) = 1
        nblock (1) = mg_ratio_x(1)*mg_ratio_y(1)*mg_ratio_z(1)
        ntotvp (1) = nxv(1)*nyv(1)*nzv(1)
        ntotv  (1) = neq*ntotvp(1)
        ntotb  (1) = nblock(1)*ntotv(1)
        do i = 2,ngrid
          nblock (i) = mg_ratio_x(i)*mg_ratio_y(i)*mg_ratio_z(i)
          ntotvp (i) = nxv(i)*nyv(i)*nzv(i)
          ntotv  (i) = neq*ntotvp(i)
          ntotb  (i) = nblock(i)*ntotv(i)
        enddo
        do i = 2,ngrid+1
          istartp(i) = istartp(i-1) + ntotvp(i-1)
          istart (i) = istart (i-1) + ntotv (i-1)
          istartb(i) = istartb(i-1) + ntotb (i-1)
        enddo

        fpointers = .true.

      endif

c     End program

      end subroutine allocPointers

c     deallocPointers
c     #################################################################
      subroutine deallocPointers(fpointers)

c     -----------------------------------------------------------------
c     Deallocates pointers for block MG if fpointers (in, logical) is
c     true.
c     -----------------------------------------------------------------

      implicit none        !For safe fortran

      logical :: fpointers

      if (fpointers) then
        deallocate(istart,ntotv,istartp,ntotvp
     .            ,istartb,ntotb,nxv,nyv,nzv,nblock
     .            ,mg_ratio_x,mg_ratio_y,mg_ratio_z)
      endif

      end subroutine deallocPointers

c     limits
c     ###############################################################
      subroutine limits(elem,nx,ny,nz,igr,imin,imax,jmin,jmax,kmin,kmax)
      implicit none
c     ---------------------------------------------------------------
c     Finds limits on loops for matvec routines. In call sequence:
c       * elem: if positive, defines lexicographic position of a
c          grid node.
c       * nx,ny,nz (in): grid dimensions
c       * imin,imax,jmin,jmax,kmin,kmax: loop limits.
c          If elem <=0, these sample the whole grid.
c          If elem > 0, these sample stencil around grid node defined
c             by elem.
c     ---------------------------------------------------------------

c     Call variables

      integer    :: elem,nx,ny,nz,imin,imax,jmin,jmax,kmin,kmax,igr

c     Local variables

      integer    :: el1

c     Begin program

      if (elem > 0) then      !Give node coordinates on grid
        el1  = mod(elem,nx*ny)
        if (el1 == 0) el1 = nx*ny

        imin = mod(el1 ,nx)
        if (imin == 0) imin = nx
        imax = imin

        jmin = 1 + (el1 - imin)/nx
        jmax = jmin

        kmin = 1 + (elem - imin - nx*(jmin-1))/(nx*ny)
        kmax = kmin
      elseif (elem <= 0) then !Give grid coordinates
        imin = 1
        imax = nx
        jmin = 1
        jmax = ny
        kmin = 1
        kmax = nz

c       Check for plane/line definitions

        if ( MGgrid%iline(igr) /= 0 ) then
          imin = MGgrid%iline(igr)
          imax = imin
        endif
        if ( MGgrid%jline(igr) /= 0 ) then
          jmin = MGgrid%jline(igr)
          jmax = jmin
        endif
        if ( MGgrid%kline(igr) /= 0 ) then
          kmin = MGgrid%kline(igr)
          kmax = kmin
        endif
      endif

c     End program

      end subroutine limits

      end module grid_mg
