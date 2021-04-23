program parareal_advance

!===========================================================================
! Computes PARAREAL new state from fine and coarse solutions according to
! prescription:
!      u_new = u_fine - u_coarse_old + u_coarse_new
!===========================================================================

  use var_io

  implicit none

  character(Len=1024) :: ifile1,ifile2,ifile3,ofile
  integer :: n_args,ntimelevels,itlevel
  integer :: ierr
  logical :: debug=.true.

  !namelist
  namelist /advance/ ntimelevels,ifile1,ifile2,ifile3,ofile

  !PIXIE3D variables
  type :: state
     type(var_array),pointer :: v1=>null(),v2=>null(),v3=>null()
  end type state
  type(var_array),pointer :: vdum=>null()

  type(state),pointer,dimension(:) :: tslice

  integer :: nxl1,nyl1,nzl1,it,itime,send_buf(3),rec_buf(3)
  real(8) :: dt,time

  integer,allocatable,dimension(:) :: itime1,itime2,itime3
  real(8),allocatable,dimension(:) ::  time1, time2, time3

! Begin program

!!$  n_args = COMMAND_ARGUMENT_COUNT()
!!$
!!$  if (n_args < 4) then
!!$     call pstop('parareal_advance','Wrong call sequence')
!!$  endif
!!$
!!$  if (debug) write (*,*) "Reading command arguments"
!!$
!!$  call get_command_argument(1, ifile1)
!!$  call get_command_argument(2, ifile2)
!!$  call get_command_argument(3, ifile3)
!!$  call get_command_argument(4, ofile)

#if defined(petsc)
  call MPI_Init(mpierr)
  call MPI_Comm_rank(MPI_COMM_WORLD,my_rank,mpierr)
  call MPI_Comm_size(MPI_COMM_WORLD,np     ,mpierr)
#if defined(adios)
  call adios_init_noxml (adios_err)
  !allocate 100MB buffer for ADIOS
  call adios_allocate_buffer (100,adios_err)
#endif
#endif

  debug = debug.and.(my_rank == 0)

  !Read configuration
  open(unit=uinput,file='advance.in',status='old')
  read(uinput,advance,iostat=ierr)
  close(uinput)

  !Open first file and read in configuration
  if (debug) write (*,*) "Reading file1 ",trim(ifile1)

  call openRestartFileForRead(file=ifile1)

  if (NDEPS > 0 .and. NDEPV > 0 .and. NAUXS > 0 .and. NAUXV > 0) then
    call setDepVarDims(NDEPS,NDEPV)
    call setAuxVarDims(NAUXS,NAUXV)
  else
    call pstop("pit_advance","No variable info available")
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
     write (*,*) 'Problem dimensions'
     write (*,*) 'NEQ=',neqd
     write (*,*) 'NX=',nxd
     write (*,*) 'NY=',nyd
     write (*,*) 'NZ=',nzd
  endif

  !Set vector dimensions and allocate variables
  call setVectorDimensions

  allocate(tslice(0:ntimelevels) &
          ,itime1(0:ntimelevels) &
          ,itime2(0:ntimelevels) &
          ,itime3(0:ntimelevels) &
          ,time1 (0:ntimelevels) &
          ,time2 (0:ntimelevels) &
          ,time3 (0:ntimelevels))

  do itlevel=0,ntimelevels
     call allocateDerivedType(tslice(itlevel)%v1)
     call allocateDerivedType(tslice(itlevel)%v2)
     call allocateDerivedType(tslice(itlevel)%v3)
  enddo

  call allocateDerivedType(vdum)

  !Read first file until last time slice
  it = 0
  do 
     ierr = readRecordFile(urecord,itime,time,dt,vdum)

     if (ierr /= 0) then
        exit
     else

        if (it > ntimelevels) then
           do itlevel=2,ntimelevels
              itime1(itlevel-1) = itime1(itlevel)
              time1 (itlevel-1) =  time1(itlevel)
              call equateDerivedType(tslice(itlevel-1)%v1,tslice(itlevel)%v1)
           enddo
        endif

        itime1(min(ntimelevels,it)) = itime
        time1 (min(ntimelevels,it)) = time
        call equateDerivedType(tslice(min(ntimelevels,it))%v1,vdum)

        it = it + 1
     endif
  enddo

  call closeRestartFileForRead

  ntimelevels = min(ntimelevels,it-1)

  if (debug) write (*,*) 'Number of available time levels=',ntimelevels

  !Read second file until last time slice
  if (debug) write (*,*) "Reading file2 ",trim(ifile2)
  nxl1 = nxl
  nyl1 = nyl
  nzl1 = nzl
  call openRestartFileForRead(file=ifile2)
  if (nxl1 /= nxl .or. nyl1 /= nyl .or. nzl1 /= nzl) then
    call pstop("parareal_advance","grid levels to not agree in v2")
  endif

  it = 0
  do 
     ierr = readRecordFile(urecord,itime,time,dt,vdum)

     if (ierr /= 0) then
        exit
     else
        if (it > ntimelevels) then
           do itlevel=2,ntimelevels
              itime2(itlevel-1) = itime2(itlevel)
              time2 (itlevel-1) =  time2(itlevel)
              call equateDerivedType(tslice(itlevel-1)%v2,tslice(itlevel)%v2)
           enddo
        endif

        itime2(min(ntimelevels,it)) = itime
        time2 (min(ntimelevels,it)) = time
        call equateDerivedType(tslice(min(ntimelevels,it))%v2,vdum)

        it = it + 1
     endif
  enddo

  call closeRestartFileForRead

  !Read third file until last time slice
  if (debug) write (*,*) "Reading file3 ",trim(ifile3)

  nxl1 = nxl
  nyl1 = nyl
  nzl1 = nzl

  call openRestartFileForRead(file=ifile3)
  if (nxl1 /= nxl .or. nyl1 /= nyl .or. nzl1 /= nzl) then
    call pstop("parareal_advance","grid levels to not agree in v3")
  endif

  it = 0
  do 
     ierr = readRecordFile(urecord,itime,time,dt,vdum)

     if (ierr /= 0) then
        exit
     else
        if (it > ntimelevels) then
           do itlevel=2,ntimelevels
              itime3(itlevel-1) = itime3(itlevel)
              time3 (itlevel-1) =  time3(itlevel)
              call equateDerivedType(tslice(itlevel-1)%v3,tslice(itlevel)%v3)
           enddo
        endif

        itime3(min(ntimelevels,it)) = itime
        time3 (min(ntimelevels,it)) = time
        call equateDerivedType(tslice(min(ntimelevels,it))%v3,vdum)

        it = it + 1
     endif
  enddo

  call closeRestartFileForRead

  !Consistency check
  if (debug) then
    do itlevel=0,ntimelevels
      write (*,*) 'Time levels=',time1(itlevel),time2(itlevel),time3(itlevel)
    enddo
  endif
  if (   (time1(ntimelevels) /= time2(ntimelevels)) &
     .or.(time2(ntimelevels) /= time3(ntimelevels))) then
    call pstop("parareal_advance","time levels to not agree")
  endif

  !Perform advance
  if (debug) write (*,*) "Performing PARAREAL advance"

  do itlevel=1,ntimelevels
    !vdum = v1-v2
    call AXPYDerivedType(1d0,tslice(itlevel)%v1,-1d0,tslice(itlevel)%v2,vdum)

    !v1 = vdum+v3
    call AXPYDerivedType(1d0,tslice(itlevel)%v3,1d0,vdum,tslice(itlevel)%v1)
  enddo

  !Write output file
  if (debug) write (*,*) "Writing output ",trim(ofile)

  do itlevel=0,ntimelevels
    if (debug) write (*,*) 'itime=',itime1(itlevel),'; time=',time1(itlevel)
    call writeRecordFile(ofile,itime1(itlevel),time1(itlevel),dt,tslice(itlevel)%v1,init=(itlevel==0))
  enddo

  !Deallocate variables
  call deallocateDerivedType(vdum)

  do itlevel=0,size(tslice)-1
    call deallocateDerivedType(tslice(itlevel)%v1)
    call deallocateDerivedType(tslice(itlevel)%v2)
    call deallocateDerivedType(tslice(itlevel)%v3)
  enddo

  deallocate(tslice,itime1,itime2,itime3,time1,time2,time3)

#if defined(petsc)
  call MPI_Finalize(mpierr)
#if defined(adios)
  call adios_finalize(my_rank,adios_err)
#endif
#endif

end program parareal_advance 
