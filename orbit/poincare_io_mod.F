      module poincare_io

        use io

        use file_io

        use grid

        use orbit

        integer :: upixie=12345,ierr

        !Grid configuration
        real(8),private :: pxmin,pxmax,pymin,pymax,pzmin,pzmax
        integer,private :: nxg,nyg,nzg
        integer,private :: bcnd(6)
        logical,private :: tor
        
        type(grid_mg_def),pointer :: g_poinc => null()

        !Global arrays
        real(8),pointer,private,dimension(:,:,:)   :: jacg
        real(8),pointer,private,dimension(:,:,:,:) :: bg,bcarg,carg

        !Spline configuration
        integer :: a_gauge=2   !Gauge (1->A1=0, 2->A2=0)
     .            ,sporder=3   !Order of splines

        !ODE configuration
        integer :: niter=8     !# Picard its. to converge in CN
     .            ,odeorder=1  !Order of ODE integration

        real(8) :: dtime=0.02d0
     .            ,tmax=0d0    !Maximum time
     .            ,smax=0d0    !Maximum arc length (in system normalized length units)
     .            ,rtol=0d-10  !Relative error tolerance in ODE convergence
     .            ,atol=1d-10  !Absolute error tolerance in CN convergence

        logical :: solen=.true.   !Whether to use A or B formulation
cc     .            ,forward=.true. !Whether to go forward or backward in field lines integration
     .            ,picard=.true.  !Whether we use Newton or Picard in CN stepping

      contains

c     read_poinc_file
c     #################################################################
      subroutine read_poinc_file(file)

        implicit none

c     -----------------------------------------------------------------
c     Read poincare file for grid and field data
c     -----------------------------------------------------------------

c     Call variables

        character(*) :: file

c     Local variables

        integer :: nx,ny,nz,i,j,k

        integer :: upixie=12345,ierr

        real(8),allocatable,dimension(:)       :: xx,yy,zz
        real(8),allocatable,dimension(:,:,:)   :: bx,by,bz,jac
        real(8),allocatable,dimension(:,:,:)   :: bxcar,bycar,bzcar
        real(8),allocatable,dimension(:,:,:,:) :: xmap

c     Begin program

        write (*,*)
        write (*,*) '*****************************************'
        write (*,*) 'Opening file ',trim(file),'...(sync)'
        write (*,*) '*****************************************'
        write (*,*)

        open(unit=upixie,file=trim(file),status='unknown'
     .      ,form='unformatted')

        read (upixie) nxg,nyg,nzg,bcnd,tor

        !Add ghost cells
        nx = nxg+2
        ny = nyg+2
        nz = nzg+2

        allocate(xx(nx),yy(ny),zz(nz))
      
        allocate(bx(nx,ny,nz),by(nx,ny,nz),bz(nx,ny,nz))

        allocate(jac  (nx,ny,nz)
     .          ,xmap (nx,ny,nz,3)
     .          ,bxcar(nx,ny,nz)
     .          ,bycar(nx,ny,nz)
     .          ,bzcar(nx,ny,nz))

        read (upixie) xx
        read (upixie) yy
        read (upixie) zz

        read (upixie,iostat=ierr) bx
        if (ierr /= 0) goto 100
        read (upixie,iostat=ierr) by
        if (ierr /= 0) goto 100
        read (upixie,iostat=ierr) bz
        if (ierr /= 0) goto 100

        read (upixie,iostat=ierr) bxcar

        if (ierr /= 0) then  !bxcar not available: assume cylindrical coords
          do k=1,nz
            do j=1,ny
              do i=1,nx
                bxcar(i,j,k) = bx(i,j,k)*cos(yy(j))
     .                       - by(i,j,k)*sin(yy(j))
                bycar(i,j,k) = bx(i,j,k)*sin(yy(j))
     .                       + by(i,j,k)*cos(yy(j))
                bzcar(i,j,k) = bz(i,j,k)

                bx(i,j,k)    = xx(i)*bx(i,j,k)
                bz(i,j,k)    = xx(i)*bz(i,j,k)

                jac(i,j,k)   = xx(i)
                xmap(i,j,k,1)= xx(i)*cos(yy(j))
                xmap(i,j,k,2)= xx(i)*sin(yy(j))
                xmap(i,j,k,3)= zz(k)
              enddo
            enddo
          enddo
        else
          read (upixie,iostat=ierr) bycar
          if (ierr /= 0) goto 100
          read (upixie,iostat=ierr) bzcar
          if (ierr /= 0) goto 100
          read (upixie,iostat=ierr) jac
          if (ierr /= 0) goto 100
          read (upixie,iostat=ierr) xmap
          if (ierr /= 0) goto 100
        endif

        close(upixie)

        !Setup orbit

        call poinc_domain_limits(xx,bcnd(1),pxmin,pxmax)
        call poinc_domain_limits(yy,bcnd(3),pymin,pymax)
        call poinc_domain_limits(zz,bcnd(5),pzmin,pzmax)

        call createGrid(nxg,nyg,nzg,pxmin,pxmax,pymin,pymax,pzmin,pzmax
     .                 ,g_poinc)
cc        call allocateGridStructure(1,1,1,1,1,1,nxg,nyg,nzg,g_poinc)
cc
cc        g_poinc%gxmin = xmin
cc        g_poinc%gxmax = xmax
cc        g_poinc%gymin = ymin
cc        g_poinc%gymax = ymax
cc        g_poinc%gzmin = zmin
cc        g_poinc%gzmax = zmax
cc
cc        g_poinc%xg = xx
cc        g_poinc%yg = yy
cc        g_poinc%zg = zz

        call orbit_setup(g_poinc,1
     .                  ,bx,by,bz,bxcar,bycar,bzcar
     .                  ,jac,xmap,bcnd
     .                  ,tmax=tmax
     .                  ,smax=smax
     .                  ,dtime=dtime
     .                  ,niter=niter
     .                  ,atol=atol
     .                  ,rtol=rtol
     .                  ,order=odeorder
     .                  ,solen=solen
     .                  ,picard=picard
     .                  ,spline_order=sporder
     .                  ,gauge=a_gauge
     .                  ,lsode_neq=5
cc     .                  ,lsode_dt_lim=.true.
     .                  )

c       Deallocate memory

        deallocate(xx,yy,zz)
        deallocate(bx,by,bz)
        deallocate(jac,xmap,bxcar,bycar,bzcar)
        call deallocateGridStructure(g_poinc)

        return

 100    write (*,*) 'Problem reading poincare file; aborting'
        stop

      end subroutine read_poinc_file

c     read_poinc_files
c     ##############################################################
      subroutine read_poinc_files(g_def,file)

        implicit none

c     --------------------------------------------------------------
c     Read poincare file for grid and field data
c     --------------------------------------------------------------

c     Call variables
        
        type(grid_mg_def),pointer :: g_def
        character(*) :: file

c     Local variables

        integer :: nfiles

c     Begin program

c     Check file configuration

        if (.not.is_file(file)) then     !Multiple files
          !Find out number of graphics files
          nfiles = 0
          do while (is_file(trim(file)//'_proc'
     .                    //trim(int2char(nfiles))))
            nfiles=nfiles+1
          enddo
        else
          nfiles = 1
        endif

c     Read POINCARE data

        if (nfiles == 0) then
          call pstop('read_poinc_files'
     .              ,'Cannot find POINCARE files')
        elseif (nfiles /= np) then
          call read_poincare_data_async(file,nfiles)
        else
          call read_poincare_data_sync(file)
        endif

c     Orbit initialization on global array

        call createGrid(nxg,nyg,nzg,pxmin,pxmax,pymin,pymax,pzmin,pzmax
     .                 ,g_def)
cc        call allocateGridStructure(1,1,1,1,1,1,nxg,nyg,nzg,g_poinc)
cc
cc        g_poinc%gxmin = xmin
cc        g_poinc%gxmax = xmax
cc        g_poinc%gymin = ymin
cc        g_poinc%gymax = ymax
cc        g_poinc%gzmin = zmin
cc        g_poinc%gzmax = zmax
cc
cc        g_poinc%xg = xx
cc        g_poinc%yg = yy
cc        g_poinc%zg = zz

        call orbit_setup(g_def,1
     .                  ,bg(:,:,:,1)
     .                  ,bg(:,:,:,2)
     .                  ,bg(:,:,:,3)
     .                  ,bcarg(:,:,:,1)
     .                  ,bcarg(:,:,:,2)
     .                  ,bcarg(:,:,:,3)
     .                  ,jacg
     .                  ,carg
     .                  ,bcnd
     .                  ,tmax=tmax
     .                  ,smax=smax
     .                  ,dtime=dtime
     .                  ,niter=niter
     .                  ,atol=atol
     .                  ,rtol=rtol
     .                  ,order=odeorder
     .                  ,solen=solen
     .                  ,picard=picard
     .                  ,spline_order=sporder
     .                  ,gauge=a_gauge
     .                  ,lsode_neq=5
!#! Marco, 06/05/2013, limited the maximum time step that can be
!performed. It's maximum value will be set to dtime.
     .                  ,lsode_dt_lim=.true.)

c     Deallocate memory

        call dealloc_poincare_data

c     End program

      end subroutine read_poinc_files

c     read_poincare_data_sync
c     ##################################################################
      subroutine read_poincare_data_sync(file)

c     ------------------------------------------------------------------
c     Reads poincare data from PIXIE3D files in a synchronous manner
c     (i.e., one file per proc).
c     ------------------------------------------------------------------

        implicit none

c     Call variables

        character(*) :: file

c     Local variables

        integer :: ilmg,ihpg,jlmg,jhpg,klmg,khpg,nx,ny,nz
     .            ,i,j,k
        real(8),pointer,dimension(:)       :: xx,yy,zz
        real(8),pointer,dimension(:,:,:)   :: jac
        real(8),pointer,dimension(:,:,:,:) :: bb,bcar,car

        character(60) :: poincfile

c     Begin program

        upixie = upixie + my_rank

        if (np > 1) then
          poincfile=trim(file)//'_proc'//trim(int2char(my_rank))
        else
          poincfile=trim(file)
        endif

        write (*,*)
        write (*,*) '***********************************************'
        write (*,*) 'Opening file ',trim(poincfile),'...(sync)'
        write (*,*) '***********************************************'
        write (*,*)

        open(unit=upixie,file=trim(poincfile),status='unknown'
     .      ,form='unformatted')

        read (upixie) nxg,nyg,nzg,bcnd,tor

c     Read poincare position arrays

        allocate(xx(0:nxg+1),yy(0:nyg+1),zz(0:nzg+1))

        read (upixie) xx
        read (upixie) yy
        read (upixie) zz

c     Set computational domain limits

        call poinc_domain_limits(xx,bcnd(1),pxmin,pxmax)
        call poinc_domain_limits(yy,bcnd(3),pymin,pymax)
        call poinc_domain_limits(zz,bcnd(5),pzmin,pzmax)

c     Local limits

        if (np > 1) then
          read (upixie) ilmg,ihpg,jlmg,jhpg,klmg,khpg
        else
          ilmg = 0
          ihpg = nxg+1

          jlmg = 0
          jhpg = nyg+1

          klmg = 0
          khpg = nzg+1
        endif

        nx = ihpg-ilmg - 1
        ny = jhpg-jlmg - 1
        nz = khpg-klmg - 1

c     Allocate local arrays

        allocate(bb  (0:nx+1,0:ny+1,0:nz+1,3)
     .          ,bcar(0:nx+1,0:ny+1,0:nz+1,3)
     .          ,jac (0:nx+1,0:ny+1,0:nz+1)
     .          ,car (0:nx+1,0:ny+1,0:nz+1,3))

c     Fill local arrays

        read (upixie,iostat=ierr) bb(:,:,:,1)
        if (ierr /= 0) goto 100
        read (upixie,iostat=ierr) bb(:,:,:,2)
        if (ierr /= 0) goto 100
        read (upixie,iostat=ierr) bb(:,:,:,3)
        if (ierr /= 0) goto 100

        read (upixie,iostat=ierr) bcar(:,:,:,1)

        if (ierr /= 0) then  !bxcar not available: assume cylindrical coords
          do k=0,nz+1
            do j=0,ny+1
              do i=0,nx+1
                bcar(i,j,k,1) = bb(i,j,k,1)*cos(yy(j))
     .                        - bb(i,j,k,2)*sin(yy(j))
                bcar(i,j,k,2) = bb(i,j,k,1)*sin(yy(j))
     .                        + bb(i,j,k,2)*cos(yy(j))
                bcar(i,j,k,3) = bb(i,j,k,3)

                bb(i,j,k,1)   = xx(i)*bb(i,j,k,1)
                bb(i,j,k,3)   = xx(i)*bb(i,j,k,3)

                jac(i,j,k)    = xx(i)
                car(i,j,k,1)  = xx(i)*cos(yy(j))
                car(i,j,k,2)  = xx(i)*sin(yy(j))
                car(i,j,k,3)  = zz(k)
              enddo
            enddo
          enddo
        else
          read (upixie,iostat=ierr) bcar(:,:,:,2)
          if (ierr /= 0) goto 100
          read (upixie,iostat=ierr) bcar(:,:,:,3)
          if (ierr /= 0) goto 100
          read (upixie,iostat=ierr) jac
          if (ierr /= 0) goto 100
          read (upixie,iostat=ierr) car
          if (ierr /= 0) goto 100
        endif

        close(upixie)

c     Allocate global arrays

        allocate(bg   (0:nxg+1,0:nyg+1,0:nzg+1,3)
     .          ,bcarg(0:nxg+1,0:nyg+1,0:nzg+1,3)
     .          ,carg (0:nxg+1,0:nyg+1,0:nzg+1,3)
     .          ,jacg (0:nxg+1,0:nyg+1,0:nzg+1))

c     Fill global arrays

        call find_global(bb,bg)
        call find_global(bcar,bcarg)

        call find_global(jac,jacg)
        call find_global(car,carg)

c     Deallocate memory

        deallocate(bb,bcar,jac,car,xx,yy,zz)

        return

 100    call pstop('read_poincare_data_sync'
     .            ,'Problem reading poincare file; aborting')

      end subroutine read_poincare_data_sync

c     read_poincare_data_async
c     ##################################################################
      subroutine read_poincare_data_async(file,nfiles)

c     ------------------------------------------------------------------
c     Reads poincare data from PIXIE3D files in an asynchronous manner
c     (i.e., poincare data is read in all procs by reading multiple files).
c     ------------------------------------------------------------------

        implicit none

c     Call variables

        integer :: nfiles

        character(*) :: file

c     Local variables

        integer :: ifl,ilmg,ihpg,jlmg,jhpg,klmg,khpg,i,j,k,ierr
        integer,pointer,dimension(:) :: mupixie
        real(8),pointer,dimension(:) :: xx,yy,zz
        character(60),allocatable,dimension(:) :: mpixie

c     Begin program

c     Setup files

        allocate(mupixie(nfiles),mpixie(nfiles))

        if (nfiles == 1) then
          mupixie(1) = upixie
          mpixie (1) = trim(file)
        else
          do ifl=1,nfiles
            mupixie(ifl) = upixie+ifl
            mpixie (ifl) = trim(file)//'_proc'//trim(int2char(ifl-1))
            write (*,*) '  => ',trim(mpixie(ifl))
     $                 ,', unit=',mupixie(ifl)
          enddo
        endif

c     Read setup info

        if (my_rank == 0) then
          write (*,*)
          write (*,*) '*****************************************'
        endif
#if defined(use_mpi) || defined(petsc)
        call MPI_Barrier(MPI_COMM_WORLD,mpierr)
#endif

        do ifl=1,nfiles

          write (*,'(a,i4,3a)') ' Proc ',my_rank,' opening file '
     .                ,trim(mpixie(ifl)),'...(async)'

          open(mupixie(ifl),file=trim(mpixie(ifl))
     $        ,form='unformatted',status='old')

          read (mupixie(ifl)) nxg,nyg,nzg,bcnd,tor

        enddo

cc        if (tor .and. (coords /= 'vmc' .and. coords /= 'tor')) then
cc          call closeFiles(mupixie)
cc          call pstop('get_global_data_from_file'
cc     .              ,'Coordinates are not toroidal')
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
cc          call pstop('get_global_data_from_file'
cc     .              ,'Global grid size does not agree')
cc        endif

c     Allocate global arrays

        allocate(xx(0:nxg+1),yy(0:nyg+1),zz(0:nzg+1))

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

          read (mupixie(ifl),iostat=ierr)
     .         bg   (ilmg:ihpg,jlmg:jhpg,klmg:khpg,1)
          if (ierr /= 0) goto 100
          read (mupixie(ifl),iostat=ierr)
     .         bg   (ilmg:ihpg,jlmg:jhpg,klmg:khpg,2)
          if (ierr /= 0) goto 100
          read (mupixie(ifl),iostat=ierr)
     .         bg   (ilmg:ihpg,jlmg:jhpg,klmg:khpg,3)
          if (ierr /= 0) goto 100

          read (mupixie(ifl),iostat=ierr)
     .         bcarg(ilmg:ihpg,jlmg:jhpg,klmg:khpg,1)

          if (ierr /= 0) then  !bxcar not available: assume cylindrical coords
            do k=klmg,khpg
              do j=jlmg,jhpg
                do i=ilmg,ihpg
                  bcarg(i,j,k,1) = bg(i,j,k,1)*cos(yy(j))
     .                           - bg(i,j,k,2)*sin(yy(j))
                  bcarg(i,j,k,2) = bg(i,j,k,1)*sin(yy(j))
     .                           + bg(i,j,k,2)*cos(yy(j))
                  bcarg(i,j,k,3) = bg(i,j,k,3)

                  bg(i,j,k,1)    = xx(i)*bg(i,j,k,1)
                  bg(i,j,k,3)    = xx(i)*bg(i,j,k,3)

                  jacg(i,j,k)    = xx(i)
                  carg(i,j,k,1)  = xx(i)*cos(yy(j))
                  carg(i,j,k,2)  = xx(i)*sin(yy(j))
                  carg(i,j,k,3)  = zz(k)
                enddo
              enddo
            enddo
          else
            read (mupixie(ifl),iostat=ierr)
     .           bcarg(ilmg:ihpg,jlmg:jhpg,klmg:khpg,2)
            if (ierr /= 0) goto 100
            read (mupixie(ifl),iostat=ierr)
     .           bcarg(ilmg:ihpg,jlmg:jhpg,klmg:khpg,3)
            if (ierr /= 0) goto 100
            read (mupixie(ifl),iostat=ierr)
     .           jacg (ilmg:ihpg,jlmg:jhpg,klmg:khpg)
            if (ierr /= 0) goto 100
            read (mupixie(ifl),iostat=ierr)
     .           carg (ilmg:ihpg,jlmg:jhpg,klmg:khpg,:)
            if (ierr /= 0) goto 100
          endif

          close(mupixie(ifl))

        enddo

c     Set computational domain limits

        call poinc_domain_limits(xx,bcnd(1),pxmin,pxmax)
        call poinc_domain_limits(yy,bcnd(3),pymin,pymax)
        call poinc_domain_limits(zz,bcnd(5),pzmin,pzmax)

        deallocate(xx,yy,zz)

c     Deallocate memory

        deallocate(mupixie)
        deallocate(mpixie)
        return

 100    call pstop('read_poincare_data_async'
     .            ,'Problem reading poincare file; aborting')

      end subroutine read_poincare_data_async

c     dealloc_poincare_data
c     ##################################################################
      subroutine dealloc_poincare_data

        deallocate(bg,bcarg,jacg,carg)
        call deallocateGridStructure(g_poinc)

      end subroutine dealloc_poincare_data

c     poinc_domain_limits
c     ##################################################################
      subroutine poinc_domain_limits(x,bc,xmin,xmax)

c     ------------------------------------------------------------------
c     Set domain limits according to dimension vector x and boundary
c     condition bc.
c     ------------------------------------------------------------------

        implicit none

c     Call variables

        integer :: bc
        real(8) :: x(:),xmin,xmax

c     Local Variables

        integer :: nn

c     Begin program

        nn = size(x)

        if (bc == PER) then     !Ghost cell is at x(1)
cc          x = x - x(2)  !Set angle = 0 at x(2)
cc          if (x(2) /= 0d0) then !PIXIE3D data, not SPECYL
cc            x = x - 0.5*(x(1)+x(2)) !Set angle = 0 at face
cc          endif
          xmin = 0d0 ; xmax = x(nn)-x(2)
cc          xmin = -0.5*(x(2)-x(1)) ; xmax = x(nn)-xmin
        else
          xmin = 5d-1*(x(1   )+x(2 ))
          xmax = 5d-1*(x(nn-1)+x(nn))
        endif

c     End program

      end subroutine poinc_domain_limits

c     get_poinc_domain_config
c     ##################################################################
      subroutine get_poinc_domain_ptys(xmin,xmax,ymin,ymax,zmin,zmax
     .                                ,bcnds,nx_g,ny_g,nz_g,toroidal)

c     ------------------------------------------------------------------
c     Extract domain properties from POINCARE file
c     ------------------------------------------------------------------

        implicit none

c     Call variables

        integer :: bcnds(6)
        real(8) :: xmin,xmax,ymin,ymax,zmin,zmax
        integer,optional :: nx_g,ny_g,nz_g
        logical,optional :: toroidal

c     Local Variables

c     Begin program

        xmin = pxmin
        xmax = pxmax
        ymin = pymin
        ymax = pymax
        zmin = pzmin
        zmax = pzmax

        bcnds = bcnd

        if (PRESENT(nx_g)) nx_g = nxg
        if (PRESENT(ny_g)) ny_g = nyg
        if (PRESENT(nz_g)) nz_g = nzg

        if (PRESENT(toroidal)) toroidal = tor

c     End program

      end subroutine get_poinc_domain_ptys

      end module poincare_io
