c module graphics
c ######################################################################
      module graphics

        use grid

        integer(4)     :: ngraph,ngroups
        parameter        (ngraph=30)

        character*(20) :: debugfile
        character*(20),allocatable,dimension(:) ::
     .                    graphfile,drawgraph,profilefile,drawprof
        character*(20),allocatable,dimension(:,:) :: prof_desc
        integer(4)    ,allocatable,dimension(:,:) :: prof_ivar,prof_log
        logical       ,allocatable,dimension(:)   :: prof_spline

        integer(4)     :: udebug
        integer(4),allocatable,dimension(:) :: ugraph,uprofile,nqty
        integer(4),allocatable,dimension(:,:) :: sel_gr

        logical        :: plot,debug

        real(8)        :: diagnostics(ngraph)

        character*(20) :: diag_desc(ngraph)

        integer(4)     :: diag_ivar(ngraph)
     .                   ,diag_log(ngraph)

        integer(4)     :: iming,imaxg,jming,jmaxg,kming,kmaxg,igroup
     .                   ,iggx,iggy,iggz,iig,jjg,kkg

        integer(4)     :: sel_diag(9),sel_graph(9),ndplot
        real(8)        :: dplot,tmplot

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

      contains

c     initializeGraphics
c     ###############################################################
      subroutine initializeGraphics (igx,igy,igz,bcond)
c     ---------------------------------------------------------------
c     Set graphics files and dumping intervals
c     ---------------------------------------------------------------

        implicit none

c     Call variables

        integer(4) :: igx,igy,igz,bcond(6)

c     Local variables

        integer(4) :: nx,ny,nz,isx,isy,isz

        integer(4) :: i,nplots

        character*(30) :: prof_title

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

        call setGraphicsRange(nx,ny,nz,bcond)

c     Find grid in logical space

        allocate(xl(iming:imaxg),yl(jming:jmaxg),zl(kming:kmaxg))

        xl(:) = grid_params%xx(iming+isx:imaxg+isx)
        yl(:) = grid_params%yy(jming+isy:jmaxg+isy)
        zl(:) = grid_params%zz(kming+isz:kmaxg+isz)

c     Define graphics i/o and initialize graph arrays (external)

        call defineGraphics

c     Initialize profile description array

        do igroup=1,ngroups
          where (graph(igroup)%array_graph%descr /= '')
            prof_desc(igroup,:) = 'Prof '
     .                            // graph(igroup)%array_graph%descr 
          elsewhere
            prof_desc(igroup,:) = ''
          end where
        enddo

c     Determine actual number of 3D plots

        do igroup=1,ngroups
          nqty(igroup) = 0
          do i = 1,size(graph(igroup)%array_graph)
            if (len(trim(graph(igroup)%array_graph(i)%descr)) == 0) exit
            nqty(igroup) = nqty(igroup) + 1
          enddo
        enddo

c     Open graphics files

        call openGraphicsFiles

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
     .                          ,'time','x','y',sel_gr(igroup,:)
     .                          ,graph(igroup)%array_graph%descr
     .                          ,'-L50',drawgraph(igroup))
        enddo

c       Profiles (selects same profiles as contour plots except vector plots)

        do igroup=1,ngroups
          nplots = count (sel_gr(igroup,:) /= 0)
          prof_title = 'Profiles '//graph(igroup)%descr
          call createDrawInGfile(nplots,profilefile(igroup),prof_title
     .          ,'x',abs(sel_gr(igroup,:)),prof_desc(igroup,:)
     .          ,prof_ivar(igroup,:),prof_log(igroup,:)
     .          ,drawprof(igroup),prof_spline(igroup))
        enddo

c     End program

      end subroutine initializeGraphics

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

        do igroup=1,ngroups
          call profplot(nqty(igroup),uprofile(igroup)
     .                 ,graph(igroup)%array_graph
     .                 ,graph(igroup)%cartesian,1,.false.)

          call contour_step(nqty(igroup),ugraph(igroup)
     .                     ,graph(igroup)%array_graph
     .                     ,graph(igroup)%cartesian,1)
        enddo

c     End program

      end subroutine dumpTimeStepPlots

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

        select case (comp)
        case(1)
          iimin = iming
          iimax = imaxg
          if (average) then
            jjmin = jming
            jjmax = jmaxg
          else
            jjmin = max(2*jmaxg/3,1)
            jjmax = max(2*jmaxg/3,1)
          endif
          kkmax = 1
          kkmin = 1

          llmin = iming
          llmax = imaxg
        case(2)

        case(3)

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
            profiles(:,igr) = profiles(:,igr)/(iimax-jjmin+1)
     .                                       /(kkmax-kkmin+1)
          case(3)
            profiles(:,igr) = profiles(:,igr)/(jjmax-jjmin+1)
     .                                       /(iimax-iimin+1)
          end select
        enddo

c     Dump plots

        do l=llmin,llmax
cc          write (*,*) yl(jmaxg/3)
          if (cartesian) then
            select case(comp)
            case(1)
              call getCartesianCoordinates(l,1,1,iggx,iggy,iggz
     .                                ,iig,jjg,kkg,car(1),car(2),car(3))
cc              car(:) = x_xi(xl(l),yl(1),zl(1))
            case(2)
              call getCartesianCoordinates(1,l,1,iggx,iggy,iggz
     .                                ,iig,jjg,kkg,car(1),car(2),car(3))
cc              car(:) = x_xi(xl(1),yl(l),zl(1))
            case(3)
              call getCartesianCoordinates(1,1,l,iggx,iggy,iggz
     .                                ,iig,jjg,kkg,car(1),car(2),car(3))
cc              car(:) = x_xi(xl(1),yl(1),zl(l))
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
c     Contours arr in xdraw format
c     Notes:
c      put the next 2 lines in main
c      open(unit=nunit,file='contour.bin',form='unformatted') before
c      close(unit=nunit) after
c     ---------------------------------------------------------------------

c     Call variables

        integer*4      nx,ny,iopt,nunit
        real*8         arr(0:nx+1,0:ny+1),xmin,xmax,ymin,ymax

c     Local variables

        integer*4      i,j

c     Begin program

        if(iopt.eq.0) then
          write(nunit) nx-1,ny-1,0
cc          write(nunit) nx+1,ny+1,0
          write(nunit) real(xmin,4),real(xmax,4)
     .                ,real(ymin,4),real(ymax,4) 
        endif
        write(nunit) ((real(arr(i,j),4),i=1,nx),j=1,ny)
cc        write(nunit) ((real(arr(i,j),4),i=0,nx+1),j=0,ny+1)

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

        k = (kmaxg-kming+1)/2
cc        k = 1

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

        call closeGraphicsFiles

        call deallocateGraphicsVariables

c     End program

      end subroutine finalizeGraphics

c     openGraphicsFiles
c     ##################################################################
      subroutine openGraphicsFiles

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

      end subroutine openGraphicsFiles

c     allocateGraphicsVariables
c     ##################################################################
      subroutine allocateGraphicsVariables(ngroups)

        implicit none

c     Call variables

        integer(4) :: ngroups

c     Begin program

        allocate(graphfile(ngroups)
     .          ,drawgraph(ngroups)
     .          ,profilefile(ngroups)
     .          ,drawprof(ngroups))

        allocate(ugraph(ngroups),uprofile(ngroups),nqty(ngroups))

        allocate(prof_desc(ngroups,ngraph),sel_gr(ngroups,9)
     .          ,prof_ivar(ngroups,ngraph),prof_log(ngroups,ngraph)
     .          ,prof_spline(ngroups))

        allocate(graph(ngroups))

      end subroutine allocateGraphicsVariables

c     deallocateGraphicsVariables
c     ##################################################################
      subroutine deallocateGraphicsVariables

        implicit none

c     Call variables

c     Begin program

        deallocate(graphfile,drawgraph,profilefile,drawprof)

        deallocate(ugraph,uprofile,nqty)

        deallocate(prof_desc,sel_gr,prof_ivar,prof_log,prof_spline)

        deallocate(graph)

      end subroutine deallocateGraphicsVariables

c     closeGraphicsFiles
c     ##################################################################
      subroutine closeGraphicsFiles

        do igroup=1,ngroups
          close(unit=ugraph(igroup))
          close(unit=uprofile(igroup))
        enddo

        if (debug) close(unit=udebug)

      end subroutine closeGraphicsFiles

c     createDrawInGfile
c     ################################################################
      subroutine createDrawInGfile(nvar,binf,graphtitle,ivar
     .                       ,dvariables,descr,ivariables,logvar
     .                       ,drawfile,spline)

c     ----------------------------------------------------------------
c     Creates the G-type (x-y curves) draw.in file for xdraw graphics.
c     On call sequence:
c       * nvar : number of variables in consideration for plotting
c       * binf : name of binary file to be read
c       * graphtitle : character array of arbitrary length
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
     .                ,logvar(nvar)
        character(*):: binf,drawfile,graphtitle,descr(ngraph),ivar
        logical     :: spline

c     Local variables

        integer      :: i,j,mfile
        character(50):: title
        logical      :: flag

c     Begin program

        mfile = 111

        open (unit=mfile,file=drawfile,status='unknown')

        write (mfile,'(1a,/)') 
     .     'Type (G=Graph, C=Contour, M=TS contour):      G'
        write (mfile,'(a)') 'filename(s)'
        write (mfile,'(a,/)') trim(binf)

        write (mfile,'(2a,/)') 'graph title: ',trim(graphtitle)
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
      subroutine createDrawInMfile(nvar,binf,graphtitle,ivart,ivarx
     .                       ,ivary,variables,descr,options,drawfile)

c     ----------------------------------------------------------------
c     Creates drawfile for xdraw graphics.
c     ----------------------------------------------------------------

        implicit none

c     Call variables

        integer     :: nvar,variables(nvar)
        character(*):: binf,drawfile,graphtitle,options,descr(ngraph)
        character(*):: ivart,ivarx,ivary

c     Local variables

        integer      :: i,j,mfile
        character(50):: title
        logical      :: flag

c     Begin program

        mfile = 111

        open (unit=mfile,file=drawfile,status='unknown')

        write (mfile,'(1a,/)') 
     .     'Type (G=Graph, C=Contour, M=TS contour):      M'
        write (mfile,'(a)') 'filename(s)'
        write (mfile,'(a,/)') trim(binf)

        write (mfile,'(2a,/)') 'comment:      ',graphtitle
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
            title = 'Vector plot (' // trim(descr(variables(i)))
     .              // ',' // trim(descr(variables(i)+1))//')'
            write (mfile,'(i2,a,i2,4a)') variables(i),' ',variables(i)+1
     .         ,'       ',options,'         ',title
          endif
        enddo
        write (mfile,*)

        close (mfile)

c     End program

      end subroutine createDrawInMfile

      end module graphics
