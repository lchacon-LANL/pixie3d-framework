c module graphics_variables
c ######################################################################
      module graphics_variables

        use grid

        use iosetup

        integer(4)     :: ngraph,ngroups
        parameter        (ngraph=30)

        real(8), allocatable, dimension(:) :: xl,yl,zl

        type :: graph_var_def
          double precision,pointer,dimension(:,:,:) :: array
          character(20)                             :: descr
        end type graph_var_def

        type :: graph_group
          type (graph_var_def),dimension(ngraph) :: array_graph
          logical                                :: cartesian
          character(20)                          :: descr
        end type graph_group

        type (graph_group),pointer,dimension(:) :: graph

        integer(4),allocatable,dimension(:) :: nqty

        logical        :: plot

        integer(4)     :: iming,imaxg,jming,jmaxg,kming,kmaxg,igroup
     .                   ,iggx,iggy,iggz,iig,jjg,kkg

        integer(4)     :: sel_diag(9),sel_graph(9)

        integer(4)     :: ndplot,iplot,jplot,kplot
        real(8)        :: dplot,tmplot

        integer(4) :: nfiles,ifile
        integer(4),allocatable,dimension(:)     :: murecord
        character*(20),allocatable,dimension(:) :: mrecfile

      end module graphics_variables

c module XDRAW
c #####################################################################
      module XDRAW

        use graphics_variables

        use variables

        character*(20) :: debugfile,axes_label(2)=(/'x','y'/)
     .                   ,prof_axis_label = 'x'
        character*(20),allocatable,dimension(:) ::
     .                    graphfile,drawgraph,profilefile,drawprof
        character*(20),allocatable,dimension(:,:) :: prof_desc
        integer(4)    ,allocatable,dimension(:,:) :: prof_ivar,prof_log
        logical       ,allocatable,dimension(:)   :: prof_spline

        integer(4)     :: udebug,prof_axis=1

        integer(4),allocatable,dimension(:) :: ugraph,uprofile

        integer(4),allocatable,dimension(:,:) :: sel_gr

      contains

c     initializeXDRAW
c     #################################################################
      subroutine initializeXDRAW

c     ---------------------------------------------------------------
c     Initializes XDRAW graphics
c     ---------------------------------------------------------------

        implicit none

c     Call variables

c     Local variables
 
        integer(4) :: nplots

        character*(30) :: prof_label

c     Begin program

c     Open XDRAW graphics files

        call openXDRAWGraphicsFiles

c     Initialize profile description array

        do igroup=1,ngroups
          where (graph(igroup)%array_graph%descr /= '')
            prof_desc(igroup,:) = 'Prof '
     .                            // graph(igroup)%array_graph%descr 
          elsewhere
            prof_desc(igroup,:) = ''
          end where
        enddo

c     Initialize contour files

        do igroup=1,ngroups
          call contour_step(nqty(igroup),ugraph(igroup)
     .                     ,graph(igroup)%array_graph
     .                     ,graph(igroup)%cartesian,0)
        enddo

c     Create draw*.in files

c       Contours

        do igroup=1,ngroups
          nplots = count (sel_gr(igroup,:) /= 0)

          call createDrawInMfile(nplots,graphfile(igroup)
     .                          ,graph(igroup)%descr
     .                          ,'time',trim(axes_label(1))
     .                          ,trim(axes_label(2)),sel_gr(igroup,:)
     .                          ,graph(igroup)%array_graph%descr
     .                          ,'-L50',drawgraph(igroup))
        enddo

c       Profiles (selects same profiles as contour plots except vector plots)

        do igroup=1,ngroups
          nplots = count (sel_gr(igroup,:) /= 0)
          prof_label = 'Profiles '//graph(igroup)%descr
          call createDrawInGfile(nplots,profilefile(igroup),prof_label
     .          ,trim(prof_axis_label),abs(sel_gr(igroup,:))
     .          ,prof_desc(igroup,:),prof_ivar(igroup,:)
     .          ,prof_log(igroup,:),drawprof(igroup)
     .          ,prof_spline(igroup))
        enddo

c     End program

      end subroutine initializeXDRAW

c     allocateXDRAWvariables
c     ##################################################################
      subroutine allocateXDRAWvariables(ngroups)

        implicit none

c     Call variables

        integer(4) :: ngroups

c     Begin program

        allocate(graphfile(ngroups)
     .          ,drawgraph(ngroups)
     .          ,profilefile(ngroups)
     .          ,drawprof(ngroups))

        allocate(ugraph(ngroups),uprofile(ngroups))

        allocate(prof_desc(ngroups,ngraph),sel_gr(ngroups,9)
     .          ,prof_ivar(ngroups,ngraph),prof_log(ngroups,ngraph)
     .          ,prof_spline(ngroups))

      end subroutine allocateXDRAWvariables

c     deallocateXDRAWvariables
c     ##################################################################
      subroutine deallocateXDRAWvariables

        implicit none

c     Call variables

c     Begin program

        deallocate(graphfile,drawgraph,profilefile,drawprof)

        deallocate(ugraph,uprofile)

        deallocate(prof_desc,sel_gr,prof_ivar,prof_log,prof_spline)

      end subroutine deallocateXDRAWvariables

c     dumpXDRAWPlots
c     #################################################################
      subroutine dumpXDRAWPlots

c     ---------------------------------------------------------------
c     Dumps time plots
c     ---------------------------------------------------------------

        use variables

        implicit none

c     Call variables

c     Local variables
 
c     Begin program

        do igroup=1,ngroups
          call profplot(nqty(igroup),uprofile(igroup)
     .                 ,graph(igroup)%array_graph
     .                 ,graph(igroup)%cartesian,prof_axis,.false.)

          call contour_step(nqty(igroup),ugraph(igroup)
     .                     ,graph(igroup)%array_graph
     .                     ,graph(igroup)%cartesian,1)
        enddo

c     End program

      end subroutine dumpXDRAWPlots

c     profplot
c     ###############################################################
      subroutine profplot(nqty,unit,array_graph,cartesian,comp,average)

c     ---------------------------------------------------------------
c     Does time plots of averaged profiles
c     ---------------------------------------------------------------

        implicit none

c     Call variables

        integer(4) :: nqty,unit,comp

        logical    :: cartesian,average

        type (graph_var_def),dimension(ngraph) :: array_graph

c     Local variables

        integer(4) :: i,j,k,iimin,iimax,jjmin,jjmax,kkmin,kkmax,igr
        integer(4) :: l,llmin,llmax
        real(8 )   :: car(3),coord
        real(8),allocatable,dimension(:,:) :: profiles

c     Begin program

        iimin = iming
        iimax = imaxg

        jjmin = jming
        jjmax = jmaxg

        kkmax = 1
        kkmin = 1

        select case (comp)
        case(1)
          if (.not.average) then
            jjmin = max(2*jmaxg/3,1)
            jjmax = max(2*jmaxg/3,1)
          endif

          llmin = iming
          llmax = imaxg
        case(2)
          if (.not.average) then
            iimin = max(2*imaxg/3,1)
            iimax = max(2*imaxg/3,1)
          endif

          llmin = jming
          llmax = jmaxg
        case(3)

          write (*,*) 'Profile direction not implemented'
          write (*,*) 'Aborting..'
          stop

        end select

        allocate(profiles(llmin:llmax,nqty))

c     Average in periodic direction(s)

        profiles = 0d0

        do igr =1,nqty
          do i=iimin,iimax
            do j=jjmin,jjmax
              do k=kkmin,kkmax
                select case(comp)
                case(1)
                  l = i
                case(2)
                  l = j
                case(3)
                  l = k
                end select
                profiles(l,igr) = profiles(l,igr)
     .                        + array_graph(igr)%array(i,j,k)
              enddo
            enddo
          enddo
          select case(comp)
          case(1)
            profiles(:,igr) = profiles(:,igr)/(jjmax-jjmin+1)
     .                                       /(kkmax-kkmin+1)
          case(2)
            profiles(:,igr) = profiles(:,igr)/(iimax-iimin+1)
     .                                       /(kkmax-kkmin+1)
          case(3)
            profiles(:,igr) = profiles(:,igr)/(jjmax-jjmin+1)
     .                                       /(iimax-iimin+1)
          end select
        enddo

c     Dump plots

        do l=llmin,llmax
          if (cartesian) then
            select case(comp)
            case(1)
              call getCartesianCoordinates(l,1,1,iggx,iggy,iggz
     .                                ,iig,jjg,kkg,car(1),car(2),car(3))
            case(2)
              call getCartesianCoordinates(1,l,1,iggx,iggy,iggz
     .                                ,iig,jjg,kkg,car(1),car(2),car(3))
            case(3)
              call getCartesianCoordinates(1,1,l,iggx,iggy,iggz
     .                                ,iig,jjg,kkg,car(1),car(2),car(3))
            end select
            write(unit) real(car(comp)),real(profiles(l,1:nqty))
          else
            select case(comp)
            case(1)
              coord = xl(l)
            case(2)
              coord = yl(l)
            case(3)
              coord = zl(l)
            end select
            write(unit) real(coord),real(profiles(l,1:nqty))
          endif
        enddo
        write(unit)

        deallocate (profiles)

c     End

      end subroutine profplot

c     contour
c     #####################################################################
      subroutine contour(arr,nx,ny,xmin,xmax,ymin,ymax,iopt,nunit)
      implicit none               !For safe fortran
c     ---------------------------------------------------------------------
c     Contours 2D array in xdraw format. In call:
c       * arr: 2D array to be plotted
c       * nx,ny: dimensions of array
c       * xmin,xmax,ymin,ymax: 2D domain limits
c       * iopt: whether to initialize xdraw plot (iopt=0) or not.
c       * nunit: integer file identifier.
c     ---------------------------------------------------------------------

c     Call variables

      integer(4) :: nx,ny,iopt,nunit
      real(8)    :: arr(nx,ny),xmin,xmax,ymin,ymax

c     Local variables

      integer(4) :: i,j

c     Begin program

      if(iopt == 0) then
        write(nunit) nx-1,ny-1,0
        write(nunit) real(xmin,4),real(xmax,4)
     .              ,real(ymin,4),real(ymax,4) 
      endif
      write(nunit) ((real(arr(i,j),4),i=1,nx),j=1,ny)

c     End program

      end subroutine contour

c     contour_step
c     #####################################################################
      subroutine contour_step(nqty,nunit,array_graph,cartsn,iopt)

c     ---------------------------------------------------------------------
c     Dumps graph arrays in array_graph in xdraw format.
c     ---------------------------------------------------------------------

        implicit none

c     Call variables

        integer(4) :: nqty,iopt,nunit

        logical    :: cartsn

        type (graph_var_def),dimension(ngraph) :: array_graph

c     Local variables

        integer(4) :: nsblk,nusblk
        integer(4) :: i,j,k,igr
        real(8),allocatable,dimension(:,:,:,:) :: car

c     Begin program

        nsblk  = 1              ! Number of structured blocks
        nusblk = 0              ! Number of unestructured blocks

c     X-Y plots

        k = max((kmaxg-kming+1)/3,1)
cc        k = kming

        if (kmaxg -kming /= 1) then
          write (*,*) 'Z-slice contour at k=',k
        endif

        if(iopt.eq.0) then

c         Write initialization headings

          write(nunit) nsblk,nusblk,nqty
          write(nunit) (imaxg-iming),(jmaxg-jming)

          if (cartsn) then      !Write cartesian coordinates

            allocate (car(iming:imaxg,jming:jmaxg,kming:kmaxg,3))

            do j=jming,jmaxg
              do i=iming,imaxg
                call getCartesianCoordinates(i,j,k,iggx,iggy,iggz
     .              ,iig,jjg,kkg,car(i,j,k,1),car(i,j,k,2),car(i,j,k,3))
              enddo
            enddo
 
            write(nunit)
     .           ((real(car(i,j,k,1),4),i=iming,imaxg),j=jming,jmaxg)
     .          ,((real(car(i,j,k,2),4),i=iming,imaxg),j=jming,jmaxg)

            deallocate(car)

          else             !Write logical coordinates

            write(nunit) ((real(xl(i),4),i=iming,imaxg),j=jming,jmaxg)
     .                  ,((real(yl(j),4),i=iming,imaxg),j=jming,jmaxg)

          endif

        else

c         Write arrays

          do igr=1,nqty
            write(nunit) ((real(array_graph(igr)%array(i,j,k),4)
     .                    ,i=iming,imaxg),j=jming,jmaxg)
          enddo

        endif

c     X-Z plots

cc        j = jming
cc
cc        if(iopt.eq.0) then
cc
cc          allocate (car(iming:imaxg,jming:jmaxg,kming:kmaxg,3))
cc
ccc         Write initialization headings
cc
cc          write(nunit) nsblk,nusblk,nqty
cc          write(nunit) (imaxg-iming),(kmaxg-kming)
cc
cc          do k=kming,kmaxg
cc            do i=iming,imaxg
cc              car(i,j,k,:) = x_xi(xl(i),yl(j),zl(k))
cc            enddo
cc          enddo
cc
cc          write(nunit) ((real(car(i,j,k,1),4),i=iming,imaxg),k=kming,kmaxg)
cc     .                ,((real(car(i,j,k,3),4),i=iming,imaxg),k=kming,kmaxg)
cc
cc          deallocate(car)
cc
cc        else
cc
cc          do igr=1,nqty
cc            write(nunit) ((real(array_graph(igr)%array(i,j,k),4)
cc     .                    ,i=iming,imaxg),k=kming,kmaxg)
cc          enddo
cc
cc        endif

c     End program

      end subroutine contour_step

c     createDrawInGfile
c     ################################################################
      subroutine createDrawInGfile(nvar,binf,graphlabel,ivar
     .                       ,dvariables,descr,ivariables,logvar
     .                       ,drawfile,spline)

c     ----------------------------------------------------------------
c     Creates the G-type (x-y curves) draw.in file for xdraw graphics.
c     On call sequence:
c       * nvar : number of variables in consideration for plotting
c       * binf : name of binary file to be read
c       * graphlabel : character array of arbitrary length
c       * ivar : independent variable axis name
c       * dvariables : integer array of size (nvar) identifying the 
c                      dependent variables (although they can also be 
c                      independent variables)
c       * descr      : character array containing dependent variables
c                      descriptions.
c       * ivariables : integer array specifying independent variables 
c                      for all the dependent variables.
c                      These include any of the dependent
c                      variable indeces plus 0 (indicating the true
c                      independent variable).
c       * logvar     : integer array specifying which dependent
c                      variables will be plotted in log scale.
c       * drawfile   : character variable specifying the name of the
c                      draw.in file.
c       * spline     : logical variable that specifies whether 
c                      xdraw should do spline interpolations or not.
c     ----------------------------------------------------------------

        implicit none

c     Call variables

        integer     :: nvar,ivariables(ngraph),dvariables(nvar)
     .                ,logvar(ngraph)
        character(*):: binf,drawfile,graphlabel,descr(ngraph),ivar
        logical     :: spline

c     Local variables

        integer      :: i,j,mfile
        character(50):: label
        logical      :: flag

c     Begin program

        mfile = 111

        open (unit=mfile,file=drawfile,status='unknown')

        write (mfile,'(1a,/)') 
     .     'Type (G=Graph, C=Contour, M=TS contour):      G'
        write (mfile,'(a)') 'filename(s)'
        write (mfile,'(a,/)') trim(binf)

        write (mfile,'(2a,/)') 'graph title: ',trim(graphlabel)
        write (mfile,'(a)') 'variable names'
        write (mfile,'(2a)') ' 0       ',trim(ivar)
        do i=1,size(descr)
          if (len(trim(descr(i))) == 0) exit
          write (mfile,'(i2,2a)') i  ,'       ',trim(descr(i))
        enddo
        write (mfile,'(/,a)') 'ix       iy        window title'

        do i=1,size(descr)
          if (len(trim(descr(i))) == 0) exit
          flag = .false.
          do j = 1,nvar
            if (dvariables(j).eq.i) then
              flag = .true.
              exit
            endif
          enddo
          if (flag) then
            if (spline) then
              if(logvar(i) == 1) then
                write (mfile,'(i2,a,i2,2a)') ivariables(i),'      '
     .                                ,i,'L&        ',trim(descr(i))
              else
                write (mfile,'(i2,a,i2,2a)') ivariables(i),'      '
     .                                ,i,'&        ',trim(descr(i))
              endif
            elseif(logvar(i) == 1) then
              write (mfile,'(i2,a,i2,2a)')  ivariables(i),'      '
     .                                ,i,'L        ',trim(descr(i))
            else
              write (mfile,'(i2,a,i2,2a)')  ivariables(i),'      '
     .                                ,i,'         ',trim(descr(i))
            endif
          else
            if(logvar(i) == 1) then
              write (mfile,'(a,i2,a,i2,2a)') ';',ivariables(i),'      '
     .                                ,i,'L        ',trim(descr(i))

            else
              write (mfile,'(a,i2,a,i2,2a)') ';',ivariables(i),'      '
     .                                ,i,'         ',trim(descr(i))
            endif
          endif
        enddo

        write (mfile,*)

        close (mfile)

c     End program

      end subroutine createDrawInGfile

c     createDrawInMfile
c     ################################################################
      subroutine createDrawInMfile(nvar,binf,graphlabel,ivart,ivarx
     .                       ,ivary,variables,descr,options,drawfile)

c     ----------------------------------------------------------------
c     Creates drawfile for "M-type" xdraw graphics.
c     ----------------------------------------------------------------

        implicit none

c     Call variables

        integer     :: nvar,variables(nvar)
        character(*):: binf,drawfile,graphlabel,options,descr(ngraph)
        character(*):: ivart,ivarx,ivary

c     Local variables

        integer      :: i,j,mfile
        character(50):: label
        logical      :: flag

c     Begin program

        mfile = 111

        open (unit=mfile,file=drawfile,status='unknown')

        write (mfile,'(1a,/)') 
     .     'Type (G=Graph, C=Contour, M=TS contour):      M'
        write (mfile,'(a)') 'filename(s)'
        write (mfile,'(a,/)') trim(binf)

        write (mfile,'(2a,/)') 'comment:      ',graphlabel
        write (mfile,'(a)') 'independent variable names'
        write (mfile,'(2a)')   't       ',trim(ivart)
        write (mfile,'(2a)')   'x       ',trim(ivarx)
        write (mfile,'(2a,/)') 'y       ',trim(ivary)
        write (mfile,'(a)') 'dependent variable names'
        do i=1,size(descr)
          if (len(trim(descr(i))) == 0) exit
          write (mfile,'(i2,2a)') i  ,'       ',descr(i)
        enddo
        write (mfile,'(/,a)') 'iqty       options     window title'

        do i=1,size(variables)
          if (variables(i) > 0) then
            write (mfile,'(i2,4a)') variables(i)
     .         ,'         ',options
     .         ,'         ',trim(descr(variables(i)))
          elseif (variables(i) == 0) then
            cycle
          else
            variables(i) = - variables(i)
            label = 'Vector plot (' // trim(descr(variables(i)))
     .              // ',' // trim(descr(variables(i)+1))//')'
            write (mfile,'(i2,a,i2,4a)') variables(i),' ',variables(i)+1
     .         ,'       ',options,'         ',trim(label)
          endif
        enddo
        write (mfile,*)

        close (mfile)

c     End program

      end subroutine createDrawInMfile

c     createDrawInCfile
c     ################################################################
      subroutine createDrawInCfile(nvar,binf,graphlabel,ivart,ivarx
     .                       ,ivary,descr,options,drawfile)

c     ----------------------------------------------------------------
c     Creates drawfile for "C-type" xdraw graphics.
c     ----------------------------------------------------------------

        implicit none

c     Call variables

        integer     :: nvar
        character(*):: binf,drawfile,graphlabel,options,descr(nvar)
        character(*):: ivart,ivarx,ivary

c     Local variables

        integer      :: i,j,cfile
        character(50):: label
        logical      :: flag

c     Begin program

        cfile = 111

        open (unit=cfile,file=drawfile,status='unknown')

        write (cfile,'(1a,/)') 
     .     'Type (G=Graph, C=Contour, M=TS contour):      C'
        write (cfile,'(a)') 'filename(s)'
        write (cfile,'(a,/)') trim(binf)

        write (cfile,'(2a,/)') 'comment:      ',graphlabel
        write (cfile,'(a)') 'variable names'
        write (cfile,'(a)')    '0        H'
        write (cfile,'(2a)')   'i0       ',trim(ivart)
        write (cfile,'(2a)')   'i1       ',trim(ivarx)
        write (cfile,'(2a,/)') 'i2       ',trim(ivary)
        write (cfile,'(a)') 'ix      iy      other           title'

        do i=1,nvar
          write (cfile,'(4a,i1,2a)') 'i2      i1'
     .         ,'     ',options,' -i0=',nvar-1
     .         ,'     ',trim(descr(i))
        enddo
        write (cfile,*)

        close (cfile)

c     End program

      end subroutine createDrawInCfile

c     openXDRAWGraphicsFiles
c     ##################################################################
      subroutine openXDRAWGraphicsFiles

        implicit none

c     Call variables

c     Begin program

        do igroup =1,ngroups
          open(unit=ugraph(igroup),file=graphfile(igroup)
     .            ,form='unformatted',status='replace')
          open(unit=uprofile(igroup),file=profilefile(igroup)
     .            ,form='unformatted',status='replace')
        enddo

        if (debug) open(unit=udebug,file=debugfile,form='unformatted')

      end subroutine openXDRAWGraphicsFiles

c     closeXDRAWGraphicsFiles
c     ##################################################################
      subroutine closeXDRAWGraphicsFiles

        do igroup=1,ngroups
          close(unit=ugraph(igroup))
          close(unit=uprofile(igroup))
        enddo

        if (debug) close(unit=udebug)

      end subroutine closeXDRAWGraphicsFiles

      end module XDRAW

c module HDF5_graphics
c ##################################################################
      module HDF5_graphics

#if defined(hdf5)
        use HDF5
#endif

        use graphics_variables

        use variables

        use timeStepping

        use generalPurposeFunctions

        character*(20) :: hdf5_fname = 'record.h5'

#if defined(hdf5)
        integer(hid_t) :: group_id,subgroup_id
        integer(hid_t) :: file_id        
        integer(hid_t) :: dataset_id     
        integer(hid_t) :: dataspace_id,node_dataspace_id
        integer(hid_t) :: attribute_id
        integer(hid_t) :: dataspace_attr_id
        integer(hid_t) :: attr_type_id

c       variables for parallel HDF5
        integer(hid_t) :: plist_id
        integer(hid_t) :: hs_dataspace_id
        integer(hid_t) :: p_dataspace_id,p_node_dataspace_id

        integer(4)     :: error,igrp

        integer(4), parameter :: ARRAY_RANK = 3
        integer(4), parameter :: ATTR_RANK = 1
        integer(4), parameter :: strlength = 20

        integer(hsize_t), dimension(ARRAY_RANK) ::
     .                         array_dim,node_array_dim
        integer(hsize_t), dimension(ATTR_RANK ) ::
     .                         attr_dim = (/ARRAY_RANK/)
        integer(4), dimension(7)                :: array_dim_i 

        character(20), dimension(ARRAY_RANK) :: coords_descr
     .                                            = (/'X','Y','Z'/)
        character(20)                        :: cell_coords_gname
     .                                            = 'cell_coords'
        character(20)                        :: node_coords_gname
     .                                            = 'node_coords'
        character(20)                        :: attr_name = 'coords'

        character(strlength), dimension(ARRAY_RANK) :: attr_data

      contains 

c     initializeHDF5file
c     ##############################################################
      subroutine initializeHDF5file

      implicit none

c     Call variables
      
c     Local variables

      integer(4) :: iattr

c     Begin program

c     Setup plot dump variables

      array_dim      = (/nxd,nyd,nzd/)
      node_array_dim = (/nxd+1,nyd+1,nzd+1/)

      !Attach location of node positions as attribute
      do iattr=1,attr_dim(1)
        attr_data(iattr) = trim('/'//trim(node_coords_gname)//
     $                     '/'//trim(coords_descr(iattr)))
      enddo

c     Exit if not in rank=0

#if defined(petsc)
      call MPI_Comm_rank(MPI_COMM_WORLD,my_rank,mpierr)
      if (my_rank /= 0) return
#endif

c     Initialize FORTRAN interface
      call h5open_f(error)

c     Create a new file using default properties.
      call h5fcreate_f(hdf5_fname, H5F_ACC_TRUNC_F, file_id, error)
      
c     Terminate access to the file.
      call h5fclose_f(file_id, error)

      call h5close_f(error)

c     End program

      end subroutine initializeHDF5file

c     createGroups
c     #################################################################
      subroutine createGroups(step)

      implicit none

c     Call variables
      
      integer(4)         :: step

c     Local variables

      integer(4)         :: ieq
      character*(20)     :: groupname

c     Begin program

#if defined(petsc)
      call MPI_Comm_rank(MPI_COMM_WORLD,my_rank,mpierr)
      if (my_rank /= 0) return
#endif

      groupname = trim('Timestep_'//int2char(step))

c     Initialize FORTRAN interface

      call h5open_f(error)

c     Create a new file using default properties.

      call h5fopen_f(hdf5_fname,H5F_ACC_RDWR_F,file_id,error,
     $               H5P_DEFAULT_F)

c     Create a group for coordinates data

      if (step == 0) then
        !Cell center positions
        call h5gcreate_f(file_id,cell_coords_gname,group_id,error,0)
        call h5gclose_f(group_id,error)
        !Node positions
        call h5gcreate_f(file_id,node_coords_gname,group_id,error,0)
        call h5gclose_f(group_id,error)
      endif

c     Create groups for computation data

      call h5gcreate_f(file_id,groupname,group_id,error, 0)

      do igrp=1,ngroups
        call h5gcreate_f(group_id,graph(igrp)%descr,subgroup_id,
     $                  error,0)
        call h5gclose_f(subgroup_id, error)
      enddo

      call h5gclose_f(group_id, error)

c     Terminate access to the file.
      call h5fclose_f(file_id, error)

      call h5close_f(error)

c     End program

      end subroutine createGroups

c     writeHDF5file
c     ##################################################################
      subroutine writeHDF5file
c     ------------------------------------------------------------------
c     Write the output date using HDF5 HyperSlab Column
c     ------------------------------------------------------------------

      implicit none

c     Call variables

c     Local variables

      integer(4)     :: i,j,k,ieq,icoords

      character*(20) :: timestep

      integer(hsize_t), dimension(ARRAY_RANK) :: count,node_count
      integer(hsize_t), dimension(ARRAY_RANK) :: offset

      real(8),allocatable,dimension(:,:,:,:) :: cell,node

c     Begin program

c     Create the output file and data groups

      timestep = trim('Timestep_'//int2char(itime))

      call createGroups(itime)

c     Initialize FORTRAN interface
      call h5open_f(error)

#if defined(petsc)
c PARALLEL
c     Set up file access property list with parallel I/O access
      call h5pcreate_f(H5P_FILE_ACCESS_F,plist_id,error)
      call h5pset_fapl_mpio_f(plist_id,MPI_COMM_WORLD,MPI_INFO_NULL,
     $                        error)

c     Open an existing file.
      call h5fopen_f(hdf5_fname,H5F_ACC_RDWR_F,file_id,error,
     $               access_prp = plist_id)
      call h5pclose_f(plist_id,error)
#else
c SERIAL
c     Open an existing file.
      call h5fopen_f(hdf5_fname,H5F_ACC_RDWR_F,file_id,error,
     $               H5P_DEFAULT_F)
#endif

c     Create the dataspace of actual data and its attributes

      !Cell-centered data
      call h5screate_simple_f(ARRAY_RANK,array_dim,dataspace_id,error)

      !Node-centered data
      call h5screate_simple_f(ARRAY_RANK,node_array_dim
     .                       ,node_dataspace_id,error)

      !Attribute data
      call h5screate_simple_f(ATTR_RANK,attr_dim,dataspace_attr_id
     $                       ,error)

#if defined(petsc)
c PARALLEL
c     Each process defines dataset in memory and writes it to the hyperslab
c     in the file.

      offset(1) = ilomg
      offset(2) = jlomg
      offset(3) = klomg

      !Cell-centered data
      count(1) = nxl
      count(2) = nyl
      count(3) = nzl

      call h5screate_simple_f(ARRAY_RANK,count,p_dataspace_id,error)

      !Node-centered data
      node_count(1) = nxl+1
      node_count(2) = nyl+1
      node_count(3) = nzl+1

      call h5screate_simple_f(ARRAY_RANK,node_count,p_node_dataspace_id
     .                       ,error)

c     Create property list for collective dataset write.
      call h5pcreate_f(H5P_DATASET_XFER_F,plist_id,error)
      call h5pset_dxpl_mpio_f(plist_id,H5FD_MPIO_COLLECTIVE_F,error)
#endif

c     Write the coordinate data

      allocate (cell(ilom:ihip,jlom:jhip,klom:khip,3)
     .         ,node(ilom:ihip,jlom:jhip,klom:khip,3))

      if (itime == 0) then

        !Cell center positions
        call h5gopen_f(file_id,cell_coords_gname,group_id,error) 

        do k=klom,khip
          do j=jlom,jhip
            do i=ilom,ihip
              call getCartesianCoordinates(i,j,k,iggx,iggy,iggz
     .          ,iig,jjg,kkg,cell(i,j,k,1),cell(i,j,k,2),cell(i,j,k,3))
            enddo
          enddo
        enddo

        do icoords=1,ARRAY_RANK

#if defined(petsc)
c PARALLEL
          call h5dcreate_f(group_id,coords_descr(icoords),
     $                     H5T_IEEE_F64BE,
     $                     dataspace_id,dataset_id,error)

c         Select hyperslab in the file.
          call h5dget_space_f(dataset_id,hs_dataspace_id,error)
          call h5sselect_hyperslab_f(hs_dataspace_id,H5S_SELECT_SET_F,
     $                               offset,count,error)    

c         Write data to the dataset
          call h5dwrite_f(dataset_id,H5T_NATIVE_DOUBLE,
     $                 cell(ilo:ihi,jlo:jhi,klo:khi,icoords),
     $                 array_dim_i,
     $                 error,file_space_id = hs_dataspace_id,
     $                 mem_space_id = p_dataspace_id,
     $                 xfer_prp = plist_id)

          call h5sclose_f(hs_dataspace_id,error)
#else
c SERIAL
          call h5dcreate_f(group_id,coords_descr(icoords),
     $                     H5T_IEEE_F64BE,
     $                     dataspace_id,dataset_id,error)

c         Write data to the dataset
          call h5dwrite_f(dataset_id,H5T_NATIVE_DOUBLE,
     $                    cell(ilo:ihi,jlo:jhi,klo:khi,icoords), !must specify index ranges
     $                    array_dim,error)
          
#endif
          call h5dclose_f(dataset_id,error)
        enddo

        call h5gclose_f(group_id, error)

        !Node positions
        call h5gopen_f(file_id,node_coords_gname,group_id,error) 

        call interpolateNodalData

        do icoords=1,ARRAY_RANK
          call h5dcreate_f(group_id,coords_descr(icoords),
     $                     H5T_IEEE_F64BE,
     $                     node_dataspace_id,dataset_id,error)

#if defined(petsc)
c PARALLEL
c         Select hyperslab in the file.
          call h5dget_space_f(dataset_id,hs_dataspace_id,error)
          call h5sselect_hyperslab_f(hs_dataspace_id,H5S_SELECT_SET_F,
     $                               offset,node_count,error)    

c         Write data to the dataset
          call h5dwrite_f(dataset_id,H5T_NATIVE_DOUBLE,
     $                 node(ilo:ihip,jlo:jhip,klo:khip,icoords),
     $                 array_dim_i,
     $                 error,
     .                 file_space_id = hs_dataspace_id,
     $                 mem_space_id = p_node_dataspace_id,
     $                 xfer_prp = plist_id)

          call h5sclose_f(hs_dataspace_id,error)
#else
c SERIAL
c         Write data to the dataset
          call h5dwrite_f(dataset_id,H5T_NATIVE_DOUBLE,
     $                    node(ilo:ihip,jlo:jhip,klo:khip,icoords),
     $                    node_array_dim,error)
#endif

          call h5dclose_f(dataset_id,error)
        enddo

        call h5gclose_f(group_id, error)

      endif

      deallocate(cell,node)

c     Write the actual data

      call h5gopen_f(file_id,timestep,group_id,error)

      do igrp=1,ngroups

c       Open an existing group of the specified file.
        call h5gopen_f(group_id,graph(igrp)%descr,subgroup_id,error)

        if (graph(igrp)%cartesian) then

c         Write "Attribute data"

          call h5tcopy_f(H5T_NATIVE_CHARACTER,attr_type_id,error)
          call h5tset_size_f(attr_type_id,strlength,error)

          call h5acreate_f(subgroup_id,attr_name,attr_type_id,
     $                     dataspace_attr_id,attribute_id,
     $                     error)

          call h5awrite_f(attribute_id,attr_type_id,
     $                    attr_data,
     $                    attr_dim,error)
          call h5aclose_f(attribute_id,error)

c         Write nodal data

          allocate (cell(ilom:ihip,jlom:jhip,klom:khip,1)
     .             ,node(ilom:ihip,jlom:jhip,klom:khip,1))

          do ieq=1,nqty(igrp)

c           Create the dataset
            call h5dcreate_f(subgroup_id
     $                      ,graph(igrp)%array_graph(ieq)%descr
     $                      ,H5T_IEEE_F64BE
     $                      ,node_dataspace_id,dataset_id,error)

c           Interpolate nodal values

            cell(:,:,:,1) = graph(igrp)%array_graph(ieq)%array

cc            write (*,*) cell(:,:,:,1)
cc            pause
            call interpolateNodalData

cc            write (*,*) node(:,:,:,1)
cc            pause

#if defined(petsc)
c PARALLEL
c           Select hyperslab in the file.
            call h5dget_space_f(dataset_id,hs_dataspace_id,error)
            call h5sselect_hyperslab_f(hs_dataspace_id,H5S_SELECT_SET_F
     $                                ,offset,node_count,error)    

c           Write data to the dataset
            call h5dwrite_f(dataset_id,H5T_NATIVE_DOUBLE
     $                     ,node(ilo:ihip,jlo:jhip,klo:khip,1)
     $                     ,array_dim_i
     $                     ,error
     $                     ,file_space_id = hs_dataspace_id
     $                     ,mem_space_id = p_node_dataspace_id
     $                     ,xfer_prp = plist_id)

c           Release resources
            call h5sclose_f(hs_dataspace_id,error)
#else
c SERIAL
c           Write data to the dataset
            call h5dwrite_f(dataset_id,H5T_NATIVE_DOUBLE
     $                     ,node(ilo:ihip,jlo:jhip,klo:khip,1)
     $                     ,node_array_dim,error)
#endif

c         Release resources

            call h5dclose_f(dataset_id,error)

          enddo

          deallocate(cell,node)

        else

c         Write cell-centered data

          do ieq=1,nqty(igrp)

c           Create the dataset
            call h5dcreate_f(subgroup_id
     $                      ,graph(igrp)%array_graph(ieq)%descr
     $                      ,H5T_IEEE_F64BE
     $                      ,dataspace_id,dataset_id,error)

#if defined(petsc)
c PARALLEL
c           Select hyperslab in the file.
            call h5dget_space_f(dataset_id,hs_dataspace_id,error)
            call h5sselect_hyperslab_f(hs_dataspace_id,H5S_SELECT_SET_F
     $                                ,offset,count,error)    

c           Write data to the dataset
            call h5dwrite_f(dataset_id,H5T_NATIVE_DOUBLE
     $                     ,graph(igrp)%array_graph(ieq)
     $                         %array(ilo:ihi,jlo:jhi,klo:khi)
     $                     ,array_dim_i
     $                     ,error,file_space_id = hs_dataspace_id
     $                     ,mem_space_id = p_dataspace_id
     $                     ,xfer_prp = plist_id)

c           Release resources
            call h5sclose_f(hs_dataspace_id,error)
#else
c SERIAL
c           Write data to the dataset
            call h5dwrite_f(dataset_id,H5T_NATIVE_DOUBLE
     $                     ,graph(igrp)%array_graph(ieq)
     $                                 %array(ilo:ihi,jlo:jhi,klo:khi) !must specify index ranges
     $                     ,array_dim,error)

#endif

c         Release resources

            call h5dclose_f(dataset_id,error)

          enddo

        endif

        call h5gclose_f(subgroup_id, error)

      enddo

c     Release resources

      call h5gclose_f(group_id, error)

#if defined(petsc)
c PARALLEL
      call h5sclose_f(p_dataspace_id,error)
      call h5sclose_f(p_node_dataspace_id,error)
      call h5pclose_f(plist_id,error)
#endif

      call h5sclose_f(dataspace_id,error)
      call h5sclose_f(node_dataspace_id,error)
      call h5fclose_f(file_id,error)
      call h5close_f (error)

c     End program

      contains

c     interpolateNodalData
c     ####################################################################
      subroutine interpolateNodalData

c     ---------------------------------------------------------------------
c     Interpolates cell-centered data to nodal data
c     ---------------------------------------------------------------------

        implicit none

        do k=klo,khip
          do j=jlo,jhip
            do i=ilo,ihip
cc              node(i,j,k,:) = 0.125
cc     .                       *(cell(i  ,j  ,k  ,:)+cell(i-1,j  ,k  ,:)
cc     .                        +cell(i  ,j-1,k  ,:)+cell(i-1,j-1,k  ,:)
cc     .                        +cell(i  ,j  ,k-1,:)+cell(i-1,j  ,k-1,:)
cc     .                        +cell(i  ,j-1,k-1,:)+cell(i-1,j-1,k-1,:))
cc              if (i == ilo .and. ilog == 1 .and. bcond(1) == SP) then
cc                node(i,j,k,1:2) = 0.25*(cell(i-1,j  ,k  ,1:2)
cc     .                                 +cell(i-1,j-1,k  ,1:2)
cc     .                                 +cell(i-1,j  ,k-1,1:2)
cc     .                                 +cell(i-1,j-1,k-1,1:2))
cc              endif

              if (i == ilo .and. ilog == 1 .and. bcond(1) == SP) then
                node(i,j,k,:) = 0.25*(cell(i-1,j  ,k  ,:)
     .                               +cell(i-1,j-1,k  ,:)
     .                               +cell(i-1,j  ,k-1,:)
     .                               +cell(i-1,j-1,k-1,:))
              else
                node(i,j,k,:) = 0.125
     .                       *(cell(i  ,j  ,k  ,:)+cell(i-1,j  ,k  ,:)
     .                        +cell(i  ,j-1,k  ,:)+cell(i-1,j-1,k  ,:)
     .                        +cell(i  ,j  ,k-1,:)+cell(i-1,j  ,k-1,:)
     .                        +cell(i  ,j-1,k-1,:)+cell(i-1,j-1,k-1,:))
              endif
            enddo
          enddo
        enddo

      end subroutine interpolateNodalData

      end subroutine writeHDF5file

#endif

      end module HDF5_graphics

c module graphics
c ######################################################################
      module graphics

        use grid

        use graphics_variables

        use XDRAW

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

        integer(4) :: igx,igy,igz,bcs(6)

c     Local variables

        integer(4) :: nx,ny,nz,isx,isy,isz

        integer(4) :: i,nplots
        
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

        call setGraphicsRange(nx,ny,nz,bcs)

c     Find grid in logical space

        allocate(xl(iming:imaxg),yl(jming:jmaxg),zl(kming:kmaxg))

        xl(:) = grid_params%xx(iming+isx:imaxg+isx)
        yl(:) = grid_params%yy(jming+isy:jmaxg+isy)
        zl(:) = grid_params%zz(kming+isz:kmaxg+isz)

c     Define graphics i/o and initialize graph arrays (external)

        call defineGraphics

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
        call initializeHDF5file()
#endif
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

#if defined(petsc)
        if (np == 1) call closeXDRAWGraphicsFiles
#else
        call closeXDRAWGraphicsFiles
#endif

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

        integer(4) :: nx,ny,nz,bcond(6)

c     Local variables
 
c     Begin program

        if (bcond(1) == PER) then
          iming = 1
          imaxg = nx+1
        elseif (bcond(1) == SP) then
cc          iming = 1
          iming = 0
          imaxg = nx+1
        else
          iming = 0
          imaxg = nx+1
        endif

        if (bcond(3) == PER) then
          jming = 1
          jmaxg = ny+1
        elseif (bcond(3) == SP) then
          jming = 1
          jmaxg = ny+1
        else
          jming = 0
          jmaxg = ny+1
        endif

        if (bcond(5) == PER) then
          kming = 1
          kmaxg = nz+1
        elseif (bcond(5) == SP) then
          kming = 1
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

#if defined(petsc)
        if (np == 1) call dumpXDRAWplots
#else
        call dumpXDRAWplots
#endif

#if defined(hdf5)
        call writeHDF5file
#endif

c     End program

      end subroutine dumpTimeStepPlots

c     allocateGraphicsVariables
c     ##################################################################
      subroutine allocateGraphicsVariables(ngroups)

        implicit none

c     Call variables

        integer(4) :: ngroups

c     Begin program

c     Allocate graph groups

        allocate(graph(ngroups),nqty(ngroups))

c     Allocate XDRAW quantities

        call allocateXDRAWvariables(ngroups)

      end subroutine allocateGraphicsVariables

c     deallocateGraphicsVariables
c     ##################################################################
      subroutine deallocateGraphicsVariables

        implicit none

c     Call variables

c     Begin program

        deallocate(graph,nqty)

#if !defined(petsc)
        call deallocateXDRAWvariables
#endif

      end subroutine deallocateGraphicsVariables

      end module graphics
