c module xdraw_graphics
c #####################################################################
      module xdraw_graphics

        use iosetup

        use xdraw_io

        use grid

        use variables

        character*(20) :: debugfile
        character*(20),allocatable,dimension(:) ::
     .                    graphfile,drawgraph,profilefile,drawprof
        character*(20),allocatable,dimension(:,:) :: prof_desc
        integer       ,allocatable,dimension(:,:) :: prof_ivar,prof_log
        logical       ,allocatable,dimension(:)   :: prof_spline

        integer        :: udebug

        integer   ,allocatable,dimension(:) :: ugraph,uprofile

        integer   ,allocatable,dimension(:,:) :: sel_gr

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
 
        integer    :: nplots,igraph

        character*(30) :: prof_label

c     Begin program

c     Open XDRAW graphics files

        call openXDRAWGraphicsFiles

c     Initialize contour files

        do igroup=1,ngroups
          if (graph(igroup)%cartesian) then
            call M_contour(ugraph(igroup),nqty(igroup)
     .                     ,graph(igroup)%array_graph(1:nqty(igroup))
     .                     ,0,car=gmetric%grid(iggx)%car)
          else
            call M_contour(ugraph(igroup),nqty(igroup)
     .                     ,graph(igroup)%array_graph(1:nqty(igroup))
     .                     ,0)
          endif
        enddo

c     Create draw*.in files

c       Contours

        do igroup=1,ngroups
          nplots = count (sel_gr(igroup,:) /= 0)

          call createDrawInMfile(nplots,graphfile(igroup)
     .                        ,graph(igroup)%cartesian
     .                        ,graph(igroup)%descr
     .                        ,'time',trim(cont_conf%label(1))
     .                        ,trim(cont_conf%label(2)),sel_gr(igroup,:)
     .                        ,graph(igroup)%array_graph%descr
     .                        ,'-L50',drawgraph(igroup))
        enddo

c       Profiles (selects same profiles as contour plots except vector plots)

        prof_desc = ''

        do igroup=1,ngroups

          do igraph=1,ngraph
            if (trim(graph(igroup)%array_graph(igraph)%descr) /='') then
              prof_desc(igroup,igraph) = 'Prof '
     .                  // trim(graph(igroup)%array_graph(igraph)%descr)
            else
              prof_desc(igroup,igraph) = ''
            endif
          enddo

cc          where (graph(igroup)%array_graph(:)%descr) /= '')
cc            prof_desc(igroup,:) = 'Prof '
cc     .                      // graph(igroup)%array_graph(:)%descr
cc          elsewhere
cc            prof_desc(igroup,:) = ''
cc          end where
cc        enddo
cc
cc        do igroup=1,ngroups

          nplots = count (sel_gr(igroup,:) /= 0)
          prof_label = 'Profiles '//graph(igroup)%descr

          call createDrawInGfile
     .          (nplots,ngraph,profilefile(igroup),prof_label
     .          ,trim(prof_conf%label),abs(sel_gr(igroup,:))
     .          ,prof_desc(igroup,:),prof_ivar(igroup,:)
     .          ,prof_log(igroup,:),drawprof(igroup)
     .          ,.true.,prof_spline(igroup))
        enddo

c     End program

      end subroutine initializeXDRAW

c     finalizeXDRAW
c     #################################################################
      subroutine finalizeXDRAW

c     ---------------------------------------------------------------
c     Finalizes XDRAW graphics
c     ---------------------------------------------------------------

        implicit none

c     Call variables

c     Local variables

c     Begin program

c     Close XDRAW graphics files

        call closeXDRAWGraphicsFiles

c     Deallocate XDRAW variables

cc        call deallocateXDRAWvariables

c     End program

      end subroutine finalizeXDRAW

c     allocateXDRAWvariables
c     ##################################################################
      subroutine allocateXDRAWvariables

        implicit none

c     Call variables

c     Begin program

        allocate(graphfile(ngroups)
     .          ,drawgraph(ngroups)
     .          ,profilefile(ngroups)
     .          ,drawprof(ngroups))

        allocate(ugraph(ngroups),uprofile(ngroups))

        allocate(prof_desc(ngroups,ngraph)
     .          ,prof_ivar(ngroups,ngraph),prof_log(ngroups,ngraph)
     .          ,prof_spline(ngroups),sel_gr(ngroups,xdraw_cont_lim))

      end subroutine allocateXDRAWvariables

c     deallocateXDRAWvariables
c     ##################################################################
      subroutine deallocateXDRAWvariables

        implicit none

c     Call variables

c     Begin program

        deallocate(graphfile,drawgraph,profilefile,drawprof)

        deallocate(ugraph,uprofile)

        deallocate(prof_desc,prof_ivar,prof_log,prof_spline,sel_gr)

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
     .                 ,graph(igroup)%cartesian,prof_conf%line,.false.)

cc          call M_contour(igroup,ugraph(igroup),nqty(igroup)
cc          call M_contour(ugraph(igroup),nqty(igroup)
cc     .                     ,graph(igroup)%array_graph(1:nqty(igroup))
cc     .                     ,graph(igroup)%cartesian,1)
          if (graph(igroup)%cartesian) then
            call M_contour(ugraph(igroup),nqty(igroup)
     .                     ,graph(igroup)%array_graph(1:nqty(igroup))
     .                     ,1,car=gmetric%grid(iggx)%car)
          else
            call M_contour(ugraph(igroup),nqty(igroup)
     .                     ,graph(igroup)%array_graph(1:nqty(igroup))
     .                     ,1)
          endif
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

        integer    :: nqty,unit,comp

        logical    :: cartesian,average

        type (graph_var_def),dimension(ngraph) :: array_graph

c     Local variables

        integer    :: i,j,k,iimin,iimax,jjmin,jjmax,kkmin,kkmax,igr
        integer    :: l,llmin,llmax,iplt,jplt,kplt
        real(8 )   :: car(3),coord
        real(8),allocatable,dimension(:,:) :: profiles

c     Begin program

        iplt = abs(prof_conf%coords(1))
        jplt = abs(prof_conf%coords(2))
        kplt = abs(prof_conf%coords(3))

        iimin = iming
        iimax = imaxg

        jjmin = jming
        jjmax = jmaxg

        kkmin = kming
        kkmax = kmaxg

        if (    (iplt < iimin .or. iplt > iimax)
     $      .or.(jplt < jjmin .or. jplt > jjmax)
     $      .or.(kplt < kkmin .or. kplt > kkmax)) then
           call pstop('profplot'
     .               ,'Input plotting coordinates out of range')
        endif

        select case (comp)
        case(1)
          if (.not.average) then
            jjmin = jplt
            jjmax = jplt

            kkmin = kplt
            kkmax = kplt
          endif

          llmin = iming
          llmax = imaxg

        case(2)
          if (.not.average) then
            iimin = iplt
            iimax = iplt

            kkmin = kplt
            kkmax = kplt
          endif

          llmin = jming
          llmax = jmaxg
        case(3)
          if (.not.average) then
            iimin = iplt
            iimax = iplt

            jjmin = jplt
            jjmax = jplt
          endif

          llmin = kming
          llmax = kmaxg

        case default

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
              call getCartesianCoordinates(l,jplt,kplt,iggx,iggy,iggz
     .                                ,iig,jjg,kkg,car(1),car(2),car(3))
            case(2)
              call getCartesianCoordinates(iplt,l,kplt,iggx,iggy,iggz
     .                                ,iig,jjg,kkg,car(1),car(2),car(3))
            case(3)
              call getCartesianCoordinates(iplt,jplt,l,iggx,iggy,iggz
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

          !Finalize M2 type contour plots
          if (graph(igroup)%cartesian)
     .         call M_contour(ugraph(igroup),nqty(igroup)
     .                    ,graph(igroup)%array_graph(1:nqty(igroup))
     .                    ,2,car=gmetric%grid(iggx)%car)

          close(unit=ugraph(igroup))
          close(unit=uprofile(igroup))
        enddo

        if (debug) close(unit=udebug)

      end subroutine closeXDRAWGraphicsFiles

      end module xdraw_graphics
