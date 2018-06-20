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
        logical,private :: tor,specyl_read=.false.

        !Global arrays
        real(8),pointer,private,dimension(:)       :: xx,yy,zz
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

c     read_poinc_files
c     ##############################################################
      subroutine read_poinc_files(file)

        implicit none

c     --------------------------------------------------------------
c     Read poincare file for grid and field data
c     --------------------------------------------------------------

c     Call variables

        character(*) :: file

c     Local variables

        integer :: nfiles,nx,ny,nz

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
          call read_poincare_data_async(file)
        else
          call read_poincare_data_sync(file)
        endif

c     Orbit initialization on global array

        !Add ghost cells
        nx = nxg+2
        ny = nyg+2
        nz = nzg+2

        call orbit_setup(nx,ny,nz,xx,yy,zz
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

        integer :: ilmg,ihpg,jlmg,jhpg,klmg,khpg,nx,ny,nz,i,j,k
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

cc        if (tor .and. (coords /= 'vmc' .and. coords /= 'tor')) then
cc          call pstop('lagrangian_init','Coordinates are not toroidal')
cc        endif
cc
cc        if (sum(bcnd-bcond) /= 0) then
cc          call pstop('lagrangian_init','BC conf. does not agree')
cc        endif

c     Read poincare position arrays

        allocate(xx(0:nxg+1),yy(0:nyg+1),zz(0:nzg+1))

        read (upixie) xx
        read (upixie) yy
        read (upixie) zz

        !Shift angles for PIXIE3D convention
        yy = yy - 0.5*(yy(0)+yy(1))
        zz = zz - 0.5*(zz(0)+zz(1))

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
          specyl_read = .true.
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
          specyl_read = .false.
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

        deallocate(bb,bcar,jac,car)

        return

 100    call pstop('read_poincare_data_sync'
     .            ,'Problem reading poincare file; aborting')

      end subroutine read_poincare_data_sync

c     read_poincare_data_async
c     ##################################################################
      subroutine read_poincare_data_async(file)

c     ------------------------------------------------------------------
c     Reads poincare data from PIXIE3D files in an asynchronous manner
c     (i.e., poincare data is read in all procs by reading multiple files).
c     ------------------------------------------------------------------

        implicit none

c     Call variables

        character(*) :: file

c     Local variables

        integer :: ifl,ilmg,ihpg,jlmg,jhpg,klmg,khpg,i,j,k,nfiles,ierr
        integer,pointer,dimension(:) :: mupixie
        character(60),allocatable,dimension(:) :: mpixie

c     Begin program

c     Setup files

        if (is_file(file)) then !Single file
          nfiles = 1
        else
          !Find out number of graphics files
          nfiles = 0
          do while (is_file(trim(file)//'_proc'
     .                    //trim(int2char(nfiles))))
            nfiles=nfiles+1
          enddo
        endif

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

#if defined(use_mpi)
        call MPI_Barrier(MPI_COMM_WORLD,mpierr)
#endif
        if (my_rank == 0) then
          write (*,*)
          write (*,*) '*****************************************'
        endif

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

        if (sum(bcnd-bcond) /= 0) then
          call closeFiles(mupixie)
          call pstop('get_global_data_from_file'
     .              ,'BC conf. does not agree')
        endif

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

          !Shift angles for PIXIE3D convention
          yy = yy - 0.5*(yy(0)+yy(1))
          zz = zz - 0.5*(zz(0)+zz(1))

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
            specyl_read = .true.
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
            specyl_read = .false.
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

c     Set domain limits

        call poinc_domain_limits(xx,bcnd(1),pxmin,pxmax)
        call poinc_domain_limits(yy,bcnd(3),pymin,pymax)
        call poinc_domain_limits(zz,bcnd(5),pzmin,pzmax)

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

        deallocate(bg,bcarg,jacg,carg,xx,yy,zz)

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

        if (bc == PER) then !Redefine periodic limits consistently with PIXIE3D
          xmin = 0d0 ; xmax = x(nn)-x(2) 
        else
          xmin = 5d-1*(x(1   )+x(2 ))
          xmax = 5d-1*(x(nn-1)+x(nn))
        endif

c     End program

      end subroutine poinc_domain_limits

c     get_poinc_domain_ptys
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
