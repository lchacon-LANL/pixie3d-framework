
c module iosetup
c ######################################################################
      module iosetup

        use io

        use variables

#if defined(adios)
        use adios_vars
ccnorb        use adios_read
#endif

        integer        :: ilevel

        integer        :: urecord=10000,usave=20000,uinput=25
        character*(50) :: recordfile='record.bin',recordsave,inputfile

        integer        :: nfiles,ifile
        integer       ,allocatable,dimension(:) :: murecord,musave
        character*(50),allocatable,dimension(:) :: mrecfile,mrecsave

        logical        :: debug,clean=.false.,tst_flg=.false.

      contains

#if defined(adios)
c     openRecordFileAdiosForRead(ierr)
c     #################################################################
      subroutine openRecordFileAdiosForRead(ierr)
        integer, intent(out) :: ierr

        integer :: gcnt, vcnt, acnt, nullidx 
        character (len=20), dimension(1)  :: gnamelist  ! list of groups (=1)

        call MPI_Comm_dup (MPI_COMM_WORLD, adios_world_comm, ierr)
        call adios_fopen (adios_fh, adiosfile, adios_world_comm, 
     .                     gcnt, ierr)
        call adios_inq_file (adios_fh,vcnt,acnt,
     .                    adios_tstart,adios_ntsteps,gnamelist,ierr)

        adios_tstep = adios_tstart + adios_ntsteps - 1
        call adios_gopen (adios_fh, adios_gh, gnamelist(1), 
     .                     adios_varcount, adios_attrcount, ierr)
        call adios_inq_group (adios_gh, adios_varlist, 
     .                       adios_attrlist, ierr)

        ! count /var/ variables
        adios_nvar=0
        do i=1,adios_varcount
            if (index(adios_varlist(i),"/var/") == 1) then
                adios_nvar = adios_nvar+1
        !        ! get rid of trailing char(0) in names
        !        nullidx=index(adios_varlist(i),char(0))
        !        if (nullidx == 0) nullidx=len(trim(adios_varlist(i)))+1
        !        adios_varnames(adios_nvar) = 
     .  !                  adios_varlist(i)(1:nullidx-1)
            endif
        enddo
        
        adios_tread = 0  ! next time read first timestep

        ! debug
        write (*,"(a,i0,a,a)"), "rank=",my_rank,": opened ", adiosfile
        if (my_rank == 0) then
            write (*,"(a, i5)") "Number of vars: ", adios_varcount
            write (*,*), "Variables:"
            do i=1,adios_varcount
                write (*,"(i5, a, a)") i,")  ", trim(adios_varlist(i))
            enddo
        endif

      end subroutine openRecordFileAdiosForRead

c     closeRecordFileAdiosForRead(ierr)
c     #################################################################
      subroutine closeRecordFileAdiosForRead(ierr)
        integer, intent(out) :: ierr
        call adios_gclose(adios_gh, ierr)
        call adios_fclose(adios_fh, ierr)
        write (*,*), "Closed ", adiosfile
      end subroutine closeRecordFileAdiosForRead
#endif

c     readRecordFile
c     #################################################################
      subroutine readRecordFile(unit,itime,time,dt,varray,ierr)

c     -----------------------------------------------------------------
c     Reads record file
c     -----------------------------------------------------------------

      implicit none

c     Call variables

      integer    :: ierr,itime(1),tmp(1),unit
      real(8)    :: time(1),dt(1)
      type(var_array),pointer:: varray

#if defined(adios)
      integer*8       :: start(1), readcount(1)  ! dimensions are 64bit in adios
      integer*8       :: m ! actually read bytes
      logical, save   :: firstread = .true.
      ierr = 0
      !firstread = (adios_tread == 0)
      write (*,"(a,i0,a,i2)"), "rank=",my_rank,
     .       " adios read time=", adios_tread
      if (adios_tread < adios_ntsteps) then
          ! read 1 item from the 1D arrays, from position adios_tread
          start = adios_tread 
          readcount = 1
          write (*,"(a,i0,a,i2,a,i2)"), "rank=",my_rank,
     .           ": read in time. start=",start,"count=",readcount
          call adios_read_var (adios_gh, "time", 
     .                     start, readcount, time, m)
          start = adios_tread 
          write (*,"(a,i0,a,i2,a,i2)"), "rank=",my_rank,
     .           ": read in itime. start=",start,"count=",readcount
          call adios_read_var (adios_gh, "itime", 
     .                     start, readcount, itime, m)
          start = adios_tread 
          write (*,"(a,i0,a,i2,a,i2)"), "rank=",my_rank,
     .           ": read in dt. start=",start,"count=",readcount
          call adios_read_var (adios_gh, "dt",  
     .                     start, readcount, dt, m)

          if (firstread) then
                start = adios_tread 
                call adios_read_var (adios_gh, "nvar", 
     .                              start, readcount, tmp, m)
                varray%nvar = tmp(1)
                ! consistency check
                if (adios_nvar.ne.varray%nvar) then

                    write (*,"(a,i0,a,i2,a,i2)"), "rank=",my_rank,
     .                     " ERROR: adios_nvar=",adios_nvar,
     .                     " != varray%nvar=",varray%nvar
                endif
          endif
          write (*,"(a,i0,a,i2,a,f6.2,a,i2)"), "rank=",my_rank,
     .         ": read itime=",itime," time=",time," nvar=",varray%nvar
          call readDerivedTypeAdios(adios_gh,varray,adios_tread,
     .           firstread,ierr)

          adios_tread = adios_tread + 1
      else  ! EOF
         ierr = -2
      endif

      firstread = .false.
      write (*,*), "Adios read returns ierr",ierr

#else
c     Begin program

      ierr = 0

      read (unit,iostat=ierr,end=100) time
      if (ierr /= 0) goto 200

      read (unit,iostat=ierr,end=100) itime
      if (ierr /= 0) goto 200

      read (unit,iostat=ierr,end=100) dt
      if (ierr /= 0) goto 200

      call readDerivedType(varray,unit,.false.,ierr)

      return

 100  ierr=-2  !EOF
cc      write (*,*) 'readRecordFile: EOF in file unit',unit
      return

 200  ierr=-1  !Error, but not EOF
      write (*,*) 'readRecordFile: Error reading file unit',unit

c     End
#endif

      end subroutine readRecordFile

c     writeRecordFile
c     #################################################################
      subroutine writeRecordFile(unit,itime,time,dt,varray)

c     -----------------------------------------------------------------
c     Writes record file
c     -----------------------------------------------------------------

      implicit none

c     Call variables

      integer    :: itime,unit
      real(8)    :: time,dt
      type(var_array),pointer:: varray

#if defined(adios)
c     ADIOS variables
      integer*8    :: handle, totalsize, groupsize
      integer      :: err
      character(2) :: mode='a'//char(0)
      logical      :: addl_write  ! true: write names and bconds


c     for calculating group size
      integer          :: xsize, ysize, zsize
      integer          :: xloghost, xhighost
      integer          :: yloghost, yhighost
      integer          :: zloghost, zhighost


c     Offsets to include ghosts
c     1. internal processes: no ghost cells written out
        xloghost = 1
        yloghost = 1
        zloghost = 1
        xhighost = -1
        yhighost = -1
        zhighost = -1
      
c     2. processes that has some face 
c       x=0 face
        if (ilomg == 0) xloghost = 0

c       x=1 face
        if (ihipg == nxd+1) xhighost = 0

c       y=0 face
        if (jlomg == 0) yloghost = 0

c       y=1 face
        if (jhipg == nyd+1) yhighost = 0

c       z=0 face
        if (klomg == 0) zloghost = 0

c       z=1 face
        if (khipg == nzd+1) zhighost = 0

c       calculate size actually written for each dimension
        xsize = ihip+xhighost-ilom-xloghost+1  
        ysize = jhip+yhighost-jlom-yloghost+1  
        zsize = khip+zhighost-klom-zloghost+1  

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c     Create/Append adios file
      addl_write=.false.
      if (adios_open_num.eq.0) then
          mode = 'w'//char(0)
          if (my_rank == 0) then
c             Process 0 (arbitrary choice) 
c             in the first timestep writes 
c             the boundary conditions and variable names
              addl_write = .true.
          endif
      else 
          mode = 'a'//char(0)
      endif

      adios_open_num = adios_open_num+1
      groupsize = 10 * 4 +     ! 12 integers
     .  2 * 8 +                ! 2 doubles + nvar*double 3D arrays
     .  varray%nvar * 8 * xsize*ysize*zsize

      if (addl_write) then
         groupsize = groupsize + 4 + 4 + 4 + ! namelen
     .     varray%nvar * 6 * 4 +  ! bconds
     .     varray%nvar * 20       ! descr (var names)
      endif

c      write (*,"('rank= ',i0,' ihip=',i0,' ilom=',i0,' jhip=',
c     . i0,' jlom=',i0,' khip=',i0,' klom=',i0)") my_rank, ihip, ilom, 
c     . jhip, jlom, khip, klom

      call adios_open (handle, "record", adiosfile, mode,
     .                 adios_world_comm, err)

c      print *, "ADIOS: adios_group_size =", groupsize
      call adios_group_size(handle, groupsize, totalsize, err)
c      print *, "ADIOS: write time", time

c     Write all variables to buffer
      if (my_rank == 0) then
          call adios_write (handle, "time", time, err)
          call adios_write (handle, "itime", itime, err)
          call adios_write (handle, "dt", dt, err)
      endif
      call writeDerivedTypeAdios(handle, varray, addl_write)

c     Close adios file which starts I/O writing from the buffer
c      print *, "ADIOS: close", time

      call adios_close (handle, err)

#else

      write (unit) time
      write (unit) itime
      write (unit) dt

      call writeDerivedType(varray,unit,.false.)
#endif

#if !defined(absoft_) && !defined(absoft)
      call flush(unit)
#endif

c     End

      end subroutine writeRecordFile

c     closeRecordFile
c     #################################################################
      subroutine closeRecordFile

c     -----------------------------------------------------------------
c     Closes record file
c     -----------------------------------------------------------------

      implicit none

c     Begin program

      close(urecord)

c     End

      end subroutine closeRecordFile

      end module iosetup
