program parareal_c2f

!===========================================================================
! PERFORMS COARSE-TO-FINE SOLUTION PROLONGATION FOR PARAREAL ALGORITHM
!===========================================================================

  use var_io

  use grid

  implicit none

  character(Len=1024) :: ifile,ofile
  integer :: n_args,ntimelevels,itlevel
  integer :: ierr,perr,grid_prolongation_factor,order
  logical :: debug=.true.,smooth

  !namelist
  namelist /c2f/ ntimelevels,ifile,ofile,grid_prolongation_factor,order,smooth

  !PIXIE3D variables
  type :: state
     type(var_array),pointer :: vf=>null(),vc=>null()
  end type state
  type(var_array),pointer :: vdum=>null()

  type(state),pointer,dimension(:) :: tslice

  integer :: i,j,k,ig,jg,kg,ivar,igx,igy,igz
  integer :: it,itm,send_buf(3),rec_buf(3)
  real(8) :: dt,tt

  integer,allocatable,dimension(:) :: itime
  real(8),allocatable,dimension(:) ::  time

  !Interpolation
  real(8) :: xp,yp,zp,interp
  real(8),allocatable,dimension(:) :: xx,yy,zz

  integer :: kx,ky,kz,nxs,nys,nzs,dim,flg

  real(8), dimension(:),allocatable:: tx,ty,tz,work
  real(8), dimension(:,:,:),allocatable:: bcoef

  real(8) :: db2val,db3val
  external   db2val,db3val

  real(8),pointer,dimension(:,:,:) :: array

! Begin program

#if defined(petsc)
  call PetscInitialize(PETSC_NULL_CHARACTER,perr)
  call initMPI(MPI_COMM_WORLD,np,my_rank)
#endif

  debug = debug.and.(my_rank == 0)

  igx = 1
  igy = 1
  igz = 1

  !Read application grid setup configuration
  call readGridInput('c2f.in')

  !Defaults
  order = 3
  grid_prolongation_factor = 1
  ntimelevels = 1
  ifile = ''
  ofile = ''
  smooth = .false.

  open(unit=uinput,file='c2f.in',status='old')
  read(uinput,c2f,iostat=ierr)
  close(uinput)

!!$#if !defined(vec_pot)
!!$  call setAuxVarDims(8,23)  !8 scalars, 22 vectors
!!$#else
!!$  call setAuxVarDims(8,26)  !8 scalars, 25 vectors
!!$#endif

  !Allocate temporal variables
  allocate(tslice(0:ntimelevels) &
          ,itime (0:ntimelevels) &
          ,time  (0:ntimelevels))

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !Open COARSE file and read in configuration!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if (debug) write (*,*) "Reading COARSE (input) file ",trim(ifile)

  call openRestartFileForRead(file=ifile)

  if (NDEPS > 0 .and. NDEPV > 0 .and. NAUXS > 0 .and. NAUXV > 0) then
    call setDepVarDims(NDEPS,NDEPV)
    call setAuxVarDims(NAUXS,NAUXV)
  else
    call pstop("pit_c2f","No variable info available")
  endif

  !Set global limits
#if defined(petsc)
  send_buf = (/ ihig,jhig,khig /)
  call MPI_Allreduce(send_buf,rec_buf,3 &
                    ,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,mpierr)

  nxd = rec_buf(1)
  nyd = rec_buf(2)
  nzd = rec_buf(3)
#else
  nxd = nxl
  nyd = nyl
  nzd = nzl
#endif

  if (debug) then
    write (*,*)
    write (*,*) 'Coarse (reference) problem dimensions'
    write (*,*) 'NEQ=',neqd
    write (*,*) 'NX=',nxd
    write (*,*) 'NY=',nyd
    write (*,*) 'NZ=',nzd
  endif

  !Set vector dimensions and allocate variables
  call allocateGlobalVar(gv)

  do itlevel=0,ntimelevels
     call allocateDerivedType(tslice(itlevel)%vc)
  enddo

  call allocateDerivedType(vdum)

  !Read first file until last time slice
  ierr = 0
  it = 0
  do 
    ierr = readRecordFile(urecord,itm,tt,dt,vdum)

    if (ierr /= 0) then
      exit
    else

      if (it > ntimelevels) then
        do itlevel=2,ntimelevels  !Equilibrium (it=0) is untouched
          itime(itlevel-1) = itime(itlevel)
          time (itlevel-1) =  time(itlevel)
          call equateDerivedType(tslice(itlevel-1)%vc,tslice(itlevel)%vc)
        enddo
      endif

      itime(min(ntimelevels,it)) = itm
      time (min(ntimelevels,it)) = tt
      call equateDerivedType(tslice(min(ntimelevels,it))%vc,vdum)

      it = it + 1

    endif
  enddo

  call closeRestartFileForRead

  ntimelevels = min(ntimelevels,it-1)

  if (ntimelevels == 0) then
    call pstop('parareal_c2f','Could not read any time level info')
  endif

  if (debug) write (*,*) 'Number of available time levels=',ntimelevels

  call deallocateDerivedType(vdum)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !Spline setup of reference (coarse) solution!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  nxs = nxl+2
  nys = nyl+2
  nzs = nzl+2

  allocate(xx(nxs),yy(nys),zz(nzs))

  call getMGmap(1,1,1,igx,igy,igz,ig,jg,kg)

  xx(1:nxs) = grid_params%xx(ig-1:ig+nxl)
  yy(1:nys) = grid_params%yy(jg-1:jg+nyl)
  zz(1:nzs) = grid_params%zz(kg-1:kg+nzl)

  flg = 0

  kx = min(order+1,nxs-1)
  ky = min(order+1,nys-1)
  kz = min(order+1,nzs-1)

  dim = nxs*nys*nzs + max(2*kx*(nxs+1),2*ky*(nys+1),2*kz*(nzs+1))

  allocate(tx(nxs+kx))
  allocate(ty(nys+ky))
  allocate(tz(nzs+kz))
  allocate(work(dim))
  allocate(bcoef(nxs,nys,nzs))

  !!!!!!!!!!!!!!!!!!!!!!
  !Perform prolongation!
  !!!!!!!!!!!!!!!!!!!!!!

  nullify(grid_params)
  call destroyGrid(gv%gparams)
      
  if (debug) write (*,*) "Performing PARAREAL coarse-to-fine prolongation"

  if (nxl > 1) nxl = nxl*grid_prolongation_factor
  if (nyl > 1) nyl = nyl*grid_prolongation_factor
  if (nzl > 1) nzl = nzl*grid_prolongation_factor

  if (nxd > 1) nxd = nxd*grid_prolongation_factor
  if (nyd > 1) nyd = nyd*grid_prolongation_factor
  if (nzd > 1) nzd = nzd*grid_prolongation_factor

  if (debug) then
    write (*,*) 'Fine (target) problem dimensions'
    write (*,*) 'NEQ=',neqd
    write (*,*) 'NX=',nxd
    write (*,*) 'NY=',nyd
    write (*,*) 'NZ=',nzd
  endif

  !Set vector dimensions and allocate variables
  call createGrid(nxd,nyd,nzd,gv%gparams)

  call setVectorDimensions

  do itlevel=0,ntimelevels
    call allocateDerivedType(tslice(itlevel)%vf)
  enddo

  !Prolongate
  if (grid_prolongation_factor == 1) then
    do itlevel=0,ntimelevels
      call equateDerivedType(tslice(itlevel)%vf,tslice(itlevel)%vc)
    enddo
  else
    do itlevel=0,ntimelevels

      do ivar=1,tslice(itlevel)%vf%nvar

        tslice(itlevel)%vf%array_var(ivar)%bconds &
             = tslice(itlevel)%vc%array_var(ivar)%bconds
        tslice(itlevel)%vf%array_var(ivar)%descr  &
             = tslice(itlevel)%vc%array_var(ivar)%descr

        tslice(itlevel)%vf%array_var(ivar)%bc_dep_list &
             = tslice(itlevel)%vc%array_var(ivar)%bc_dep_list
        tslice(itlevel)%vf%array_var(ivar)%dom_dep_list &
             = tslice(itlevel)%vc%array_var(ivar)%dom_dep_list

        !Spline reference component
        array => tslice(itlevel)%vc%array_var(ivar)%array

        call db3ink(xx,nxs,yy,nys,zz,nzs,array &
                   ,nxs,nys,kx,ky,kz,tx,ty,tz,bcoef,work,flg)

        !Interpolate to solution grid and compute L2-norm
        array => tslice(itlevel)%vf%array_var(ivar)%array

        do k = 0,nzl+1
          do j = 0,nyl+1
            do i = 0,nxl+1
              call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)

              xp = grid_params%xx(ig)
              yp = grid_params%yy(jg)
              zp = grid_params%zz(kg)

              array(i,j,k) = db3val(xp,yp,zp,0,0,0,tx,ty,tz,nxs,nys,nzs &
                                   ,kx,ky,kz,bcoef,work)
            enddo
          enddo
        enddo

        !Apply binomial filter
        if (smooth) then
          array(1:nxl,1:nyl,1:nzl) = 0.25*(array(0:nxl-1,1:nyl,1:nzl) &
                                          +array(2:nxl+1,1:nyl,1:nzl) &
                                      +2d0*array(1:nxl  ,1:nyl,1:nzl)) 
          array(1:nxl,1:nyl,1:nzl) = 0.25*(array(1:nxl,0:nyl-1,1:nzl) &
                                          +array(1:nxl,2:nyl+1,1:nzl) &
                                      +2d0*array(1:nxl,1:nyl  ,1:nzl))
          array(1:nxl,1:nyl,1:nzl) = 0.25*(array(1:nxl,1:nyl,0:nzl-1) &
                                          +array(1:nxl,1:nyl,2:nzl+1) &
                                      +2d0*array(1:nxl,1:nyl,1:nzl  ))
	endif

      enddo
    enddo

  endif

  !!!!!!!!!!!!!!!!!!!!!!!!!!
  !Write FINE (output) file!
  !!!!!!!!!!!!!!!!!!!!!!!!!!

  if (debug) write (*,*) "Writing output ",trim(ofile)
  do itlevel=0,ntimelevels
    if (debug) write (*,*) 'itime=',itime(itlevel),'; time=',time(itlevel)
    call writeRecordFile(ofile,itime(itlevel),time(itlevel),dt,tslice(itlevel)%vf,init=(itlevel==0))
  enddo

  !!!!!!!!!!!!!!!!!!!!!!
  !Deallocate variables!
  !!!!!!!!!!!!!!!!!!!!!!

  do itlevel=0,size(tslice)-1
    call deallocateDerivedType(tslice(itlevel)%vc)
    call deallocateDerivedType(tslice(itlevel)%vf)
  enddo
  deallocate(tslice,itime,time)
  deallocate(tx,ty,tz,work,bcoef,xx,yy,zz)
  nullify(array)

  call deallocateStructures
  call deallocateGlobalVar(gv)

#if defined(petsc)
  call PetscFinalize(mpierr)
#endif

end program parareal_c2f 

