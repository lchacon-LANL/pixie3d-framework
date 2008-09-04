c module xdraw_graphics
c #####################################################################
      module xdraw_graphics

        use iosetup

        use graphics_variables

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
          call contour_step(igroup,ugraph(igroup),nqty(igroup)
     .                     ,graph(igroup)%array_graph
     .                     ,graph(igroup)%cartesian,0)
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

          call createDrawInGfile(nplots,profilefile(igroup),prof_label
     .          ,trim(prof_conf%label),abs(sel_gr(igroup,:))
     .          ,prof_desc(igroup,:),prof_ivar(igroup,:)
     .          ,prof_log(igroup,:),drawprof(igroup)
     .          ,prof_spline(igroup))
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

          call contour_step(igroup,ugraph(igroup),nqty(igroup)
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

c     contour_step
c     #####################################################################
      subroutine contour_step(igroup,nunit,nqty,array_graph,cartsn,iopt)

c     ---------------------------------------------------------------------
c     Dumps graph arrays in array_graph in M xdraw format.
c     ---------------------------------------------------------------------

        implicit none

c     Call variables

        integer    :: nqty,iopt,nunit,igroup

        logical    :: cartsn

        type (graph_var_def),dimension(ngraph) :: array_graph

c     Local variables

        integer    :: nsblk,nusblk
        integer    :: i,j,k,igr,plt(3)
        logical    :: cut(3)

c     Begin program

        nsblk  = 1              ! Number of structured blocks
        nusblk = 0              ! Number of unestructured blocks

c     Select 2D cuts

        plt = cont_conf%coords

        cut = .false.
        cut(cont_conf%plane) = .true.

c     Perform 2d cuts

        if (cut(1)) then
          call twodcut(plt(1),plt(1),plt(1),jming,jmaxg,kming,kmaxg
     .                ,'X-slice contour at i=')
        elseif (cut(2)) then
          call twodcut(plt(2),iming,imaxg,plt(2),plt(2),kming,kmaxg
     .                ,'Y-slice contour at j=')
        else
          call twodcut(plt(3),iming,imaxg,jming,jmaxg,plt(3),plt(3)
     .                ,'Z-slice contour at k=')
        endif

c     End program

      contains

c     twodcut
c     ################################################################
      subroutine twodcut(iplt,imng,imxg,jmng,jmxg,kmng,kmxg,messg)

        implicit none

c     Call variables

        integer      :: iplt,imng,imxg,jmng,jmxg,kmng,kmxg,i1,i2
        character(*) :: messg

c     Local variables

        integer    :: nxgr,nygr,nxgr_old=0
     .               ,nygr_old=0,zero=0,one=1
     .               ,nt=1,nt_prev=0

        logical    :: flag(ngroups)

c     Begin program

c     Write messg to standard output

        write (*,*) messg,iplt

        if (cartsn) then   !Type M2 contours (physical space)

          if (iopt == 0) then

c         Write initialization headings

            write(nunit) nsblk,nusblk,nqty

            flag(igroup) = .true.

          elseif (iopt == 1) then

            if (imng == imxg) then
              nxgr = jmxg-jmng
              nygr = kmxg-kmng
              call find_cart_plane(1,i1,i2)
cc              i1 = 2
cc              i2 = 3
            elseif (jmng == jmxg) then
              nxgr = imxg-imng
              nygr = kmxg-kmng
              call find_cart_plane(2,i1,i2)
cc              i1 = 1
cc              i2 = 3
            else
              nxgr = imxg-imng
              nygr = jmxg-jmng
              call find_cart_plane(3,i1,i2)
cc              i1 = 1
cc              i2 = 2
            endif

c         Write data

            if (flag(igroup)) then
              write(nunit) zero,nxgr,nygr,nt,zero,zero,zero
              flag(igroup) = .false.
            else
              write(nunit) zero,nxgr,nygr,nt,nxgr_old,nygr_old,nt_prev
            endif

            call write_coords(imng,imxg,jmng,jmxg,kmng,kmxg,i1,i2)

            write(nunit) one
            do igr=1,nqty
              write(nunit) (((real(array_graph(igr)%array(i,j,k),4)
     .                    ,i=imng,imxg),j=jmng,jmxg),k=kmng,kmxg)
            enddo

            nxgr_old = nxgr
            nygr_old = nygr
            nt_prev  = nt

          else

            write(nunit) zero,zero,zero,zero,nxgr_old,nygr_old,nt_prev

          endif

        else   !Type M contours (logical space)

          if(iopt == 0) then

c         Write initialization headings

            write(nunit) nsblk,nusblk,nqty

            if (imng == imxg) then
              write(nunit) (jmxg-jmng),(kmxg-kmng)
              i1 = 2
              i2 = 3
            elseif (jmng == jmxg) then
              write(nunit) (imxg-imng),(kmxg-kmng)
              i1 = 1
              i2 = 3
            else
              write(nunit) (imxg-imng),(jmxg-jmng)
              i1 = 1
              i2 = 2
            endif

            call write_coords(imng,imxg,jmng,jmxg,kmng,kmxg,i1,i2)

          else

c         Write data

            do igr=1,nqty
              write(nunit) (((real(array_graph(igr)%array(i,j,k),4)
     .                    ,i=imng,imxg),j=jmng,jmxg),k=kmng,kmxg)
            enddo

          endif

        endif

      end subroutine twodcut

c     write_coords
c     ################################################################
      subroutine write_coords(imng,imxg,jmng,jmxg,kmng,kmxg,i1,i2)

        implicit none

c     Call variables

        integer    :: imng,imxg,jmng,jmxg,kmng,kmxg,i1,i2

c     Local variables

        integer    :: ii

cc        real(8),allocatable,dimension(:,:,:,:) :: car
        real(8),pointer,dimension(:,:,:,:) :: car

c     Begin program

        allocate (car(imng:imxg,jmng:jmxg,kmng:kmxg,3))

        if (cartsn) then      !Write cartesian coordinates

          do k=kmng,kmxg
            do j=jmng,jmxg
              do i=imng,imxg
                call getCartesianCoordinates(i,j,k,iggx,iggy,iggz
     .            ,iig,jjg,kkg,car(i,j,k,1),car(i,j,k,2),car(i,j,k,3))
              enddo
            enddo
          enddo

        else             !Write logical coordinates

          do k=kmng,kmxg
            do j=jmng,jmxg
              do i=imng,imxg
                car(i,j,k,1) = xl(i)
                car(i,j,k,2) = yl(j)
                car(i,j,k,3) = zl(k)
              enddo
            enddo
          enddo

        endif

        write(nunit) (((real(car(i,j,k,i1),4),i=imng,imxg)
     .                                       ,j=jmng,jmxg)
     .                                       ,k=kmng,kmxg)
     .              ,(((real(car(i,j,k,i2),4),i=imng,imxg)
     .                                       ,j=jmng,jmxg)
     .                                       ,k=kmng,kmxg)

        deallocate(car)

      end subroutine write_coords

c     find_cart_plane
c     ################################################################
      subroutine find_cart_plane(iplane,i1,i2)

c     ----------------------------------------------------------------
c     Determines cartesian components that correspond to a logical
c     plane of interest (defined by iplane)
c     ----------------------------------------------------------------

        implicit none

c     Call variables

        integer    :: iplane,i1,i2

c     Local variables

        integer    :: ii
        real(8)    :: max_x,max_y,max_z,min_x,min_y,min_z

        real(8),pointer,dimension(:,:,:,:) :: car

c     Begin program

        car => gmetric%grid(iggx)%car

        if (iplane == 1) then   !We take i=1 as a sample; the result is valid for arbitrary i.
          max_x = maxval(car(1,:,:,1))
          max_y = maxval(car(1,:,:,2))
          max_z = maxval(car(1,:,:,3))

          min_x = minval(car(1,:,:,1))
          min_y = minval(car(1,:,:,2))
          min_z = minval(car(1,:,:,3))
        elseif (iplane == 2) then   !Same for j
          max_x = maxval(car(:,1,:,1))
          max_y = maxval(car(:,1,:,2))
          max_z = maxval(car(:,1,:,3))

          min_x = minval(car(:,1,:,1))
          min_y = minval(car(:,1,:,2))
          min_z = minval(car(:,1,:,3))
        elseif (iplane == 3) then   !Same for k
          max_x = maxval(car(:,:,1,1))
          max_y = maxval(car(:,:,1,2))
          max_z = maxval(car(:,:,1,3))
                                  
          min_x = minval(car(:,:,1,1))
          min_y = minval(car(:,:,1,2))
          min_z = minval(car(:,:,1,3))
        endif

        if (max_x == min_x) then
          i1 = 2
          i2 = 3
        endif

        if (max_y == min_y) then
          i1 = 1
          i2 = 3
        endif

        if (max_z == min_z) then
          i1 = 1
          i2 = 2
        endif

      end subroutine find_cart_plane

ccc     twodcut
ccc     ################################################################
cc      subroutine twodcut(iplt,imng,imxg,jmng,jmxg,kmng,kmxg,messg)
cc
cc        implicit none
cc
cc        integer      :: normal,iplt,imng,imxg,jmng,jmxg,kmng,kmxg,i1,i2
cc        character(*) :: messg
cc
cc        write (*,*) messg,iplt
cc
cc        if(iopt.eq.0) then
cc
ccc         Write initialization headings
cc
cc          write(nunit) nsblk,nusblk,nqty
cc          if (imng == imxg) then
cc            write(nunit) (jmxg-jmng),(kmxg-kmng)
cc            i1 = 2
cc            i2 = 3
cc          elseif (jmng == jmxg) then
cc            write(nunit) (imxg-imng),(kmxg-kmng)
cc            i1 = 1
cc            i2 = 3
cc          else
cc            write(nunit) (imxg-imng),(jmxg-jmng)
cc            i1 = 1
cc            i2 = 2
cc          endif
cc
cc          allocate (car(imng:imxg,jmng:jmxg,kmng:kmxg,3))
cc
cc          if (cartsn) then      !Write cartesian coordinates
cc
cc            do k=kmng,kmxg
cc              do j=jmng,jmxg
cc                do i=imng,imxg
cc                  call getCartesianCoordinates(i,j,k,iggx,iggy,iggz
cc     .              ,iig,jjg,kkg,car(i,j,k,1),car(i,j,k,2),car(i,j,k,3))
cc                enddo
cc              enddo
cc            enddo
cc
cc          else             !Write logical coordinates
cc
cc            do k=kmng,kmxg
cc              do j=jmng,jmxg
cc                do i=imng,imxg
cc                  car(i,j,k,1) = xl(i)
cc                  car(i,j,k,2) = yl(j)
cc                  car(i,j,k,3) = zl(k)
cc                enddo
cc              enddo
cc            enddo
cc
cc          endif
cc
cc          write(nunit) (((real(car(i,j,k,i1),4),i=imng,imxg)
cc     .                                         ,j=jmng,jmxg)
cc     .                                         ,k=kmng,kmxg)
cc     .                ,(((real(car(i,j,k,i2),4),i=imng,imxg)
cc     .                                         ,j=jmng,jmxg)
cc     .                                         ,k=kmng,kmxg)
cc
cc          deallocate(car)
cc
cc        else
cc
cc          do igr=1,nqty
cc            write(nunit) (((real(array_graph(igr)%array(i,j,k),4)
cc     .                    ,i=imng,imxg),j=jmng,jmxg),k=kmng,kmxg)
cc          enddo
cc
cc        endif
cc
cc      end subroutine twodcut

      end subroutine contour_step

ccc     contour
ccc     #####################################################################
cc      subroutine contour(arr,nx,ny,xmin,xmax,ymin,ymax,iopt,nunit)
cc      implicit none               !For safe fortran
ccc     ---------------------------------------------------------------------
ccc     Contours 2D array in xdraw format. In call:
ccc       * arr: 2D array to be plotted
ccc       * nx,ny: dimensions of array
ccc       * xmin,xmax,ymin,ymax: 2D domain limits
ccc       * iopt: whether to initialize xdraw plot (iopt=0) or not.
ccc       * nunit: integer file identifier.
ccc     ---------------------------------------------------------------------
cc
ccc     Call variables
cc
cc      integer    :: nx,ny,iopt,nunit
cc      real(8)    :: arr(nx,ny),xmin,xmax,ymin,ymax
cc
ccc     Local variables
cc
cc      integer    :: i,j
cc
ccc     Begin program
cc
cc      if(iopt == 0) then
cc        write(nunit) nx-1,ny-1,0
cc        write(nunit) real(xmin,4),real(xmax,4)
cc     .              ,real(ymin,4),real(ymax,4) 
cc      endif
cc      write(nunit) ((real(arr(i,j),4),i=1,nx),j=1,ny)
cc
ccc     End program
cc
cc      end subroutine contour

ccc     createDrawInGfile
ccc     ################################################################
cc      subroutine createDrawInGfile(nvar,binf,graphlabel,ivar
cc     .                       ,dvariables,descr,ivariables,logvar
cc     .                       ,drawfile,spline)
cc
ccc     ----------------------------------------------------------------
ccc     Creates the G-type (x-y curves) draw.in file for xdraw graphics.
ccc     On call sequence:
ccc       * nvar : number of variables in consideration for plotting
ccc       * binf : name of binary file to be read
ccc       * graphlabel : character array of arbitrary length
ccc       * ivar : independent variable axis name
ccc       * dvariables : integer array of size (nvar) identifying the 
ccc                      dependent variables (although they can also be 
ccc                      independent variables)
ccc       * descr      : character array containing dependent variables
ccc                      descriptions.
ccc       * ivariables : integer array specifying independent variables 
ccc                      for all the dependent variables.
ccc                      These include any of the dependent
ccc                      variable indeces plus 0 (indicating the true
ccc                      independent variable).
ccc       * logvar     : integer array specifying which dependent
ccc                      variables will be plotted in log scale.
ccc       * drawfile   : character variable specifying the name of the
ccc                      draw.in file.
ccc       * spline     : logical variable that specifies whether 
ccc                      xdraw should do spline interpolations or not.
ccc     ----------------------------------------------------------------
cc
cc        implicit none
cc
ccc     Call variables
cc
cc        integer      :: nvar,ivariables(ngraph),dvariables(nvar)
cc     .                 ,logvar(ngraph)
cc        character(*) :: binf,drawfile,graphlabel,ivar
cc        character(20):: descr(ngraph)
cc        logical      :: spline
cc
ccc     Local variables
cc
cc        integer      :: i,j,gfile,lngt
cc        character(50):: label
cc        logical      :: flag
cc
ccc     Begin program
cc
cc        gfile = 111
cc
cc        open (unit=gfile,file=drawfile,status='unknown')
cc
cc        write (gfile,'(1a,/)') 
cc     .     'Type (G=Graph, C=Contour, M=TS contour):      G'
cc        write (gfile,'(a)') 'filename(s)'
cc        write (gfile,'(a,/)') trim(binf)
cc
cc        write (gfile,'(2a,/)') 'graph title: ',trim(graphlabel)
cc        write (gfile,'(a)') 'variable names'
cc        write (gfile,'(2a)') ' 0       ',trim(ivar)
cc
cc        do i=1,size(descr)
cc          lngt = len(trim(descr(i)))
cc          if (lngt == 0) cycle
cc          write (gfile,'(i2,2a)') i  ,'       ',trim(descr(i))
cc        enddo
cc        write (gfile,'(/,a)') 'ix       iy        window title'
cc
cc        do i=1,size(descr)
cc          lngt = len(trim(descr(i)))
cc          if (lngt == 0) exit
cc          flag = .false.
cc          do j = 1,nvar
cc            if (dvariables(j).eq.i) then
cc              flag = .true.
cc              exit
cc            endif
cc          enddo
cc          if (flag) then
cc            if (spline) then
cc              if(logvar(i) == 1) then
cc                write (gfile,'(i2,a,i2,2a)') ivariables(i),'      '
cc     .                                ,i,'L&        ',trim(descr(i))
cc              else
cc                write (gfile,'(i2,a,i2,2a)') ivariables(i),'      '
cc     .                                ,i,'&        ',trim(descr(i))
cc              endif
cc            elseif(logvar(i) == 1) then
cc              write (gfile,'(i2,a,i2,2a)')  ivariables(i),'      '
cc     .                                ,i,'L        ',trim(descr(i))
cc            else
cc              write (gfile,'(i2,a,i2,2a)')  ivariables(i),'      '
cc     .                                ,i,'         ',trim(descr(i))
cc            endif
cc          else
cc            if(logvar(i) == 1) then
cc              write (gfile,'(a,i2,a,i2,2a)') ';',ivariables(i),'      '
cc     .                                ,i,'L        ',trim(descr(i))
cc
cc            else
cc              write (gfile,'(a,i2,a,i2,2a)') ';',ivariables(i),'      '
cc     .                                ,i,'         ',trim(descr(i))
cc            endif
cc          endif
cc        enddo
cc
cc        write (gfile,*)
cc
cc        close (gfile)
cc
ccc     End program
cc
cc      end subroutine createDrawInGfile
cc
ccc     createDrawInMfile
ccc     ################################################################
cc      subroutine createDrawInMfile(nvar,binf,cartesian,graphlabel
cc     .            ,ivart,ivarx,ivary,variables,descr,options,drawfile)
cc
ccc     ----------------------------------------------------------------
ccc     Creates drawfile for "M-type" xdraw graphics.
ccc     ----------------------------------------------------------------
cc
cc        implicit none
cc
ccc     Call variables
cc
cc        integer     :: nvar,variables(nvar)
cc        character(*):: binf,drawfile,graphlabel,options,descr(ngraph)
cc        character(*):: ivart,ivarx,ivary
cc        logical     :: cartesian
cc
ccc     Local variables
cc
cc        integer      :: i,j,mfile,lngt
cc        character(50):: label
cc        logical      :: flag
cc
ccc     Begin program
cc
cc        mfile = 111
cc
cc        open (unit=mfile,file=drawfile,status='unknown')
cc
cc        if (cartesian) then
cc          write (mfile,'(1a,/)') 
cc     .     'Type (G=Graph, C=Contour, M=TS contour):      M2'
cc        else
cc          write (mfile,'(1a,/)') 
cc     .     'Type (G=Graph, C=Contour, M=TS contour):      M'
cc        endif
cc        write (mfile,'(a)') 'filename(s)'
cc        write (mfile,'(a,/)') trim(binf)
cc
cc        write (mfile,'(2a,/)') 'comment:      ',graphlabel
cc        write (mfile,'(a)') 'independent variable names'
cc        write (mfile,'(2a)')   't       ',trim(ivart)
cc        write (mfile,'(2a)')   'x       ',trim(ivarx)
cc        write (mfile,'(2a,/)') 'y       ',trim(ivary)
cc        write (mfile,'(a)') 'dependent variable names'
cc        do i=1,size(descr)
cc          lngt = len(trim(descr(i)))
cc          if (lngt == 0) exit
cc          write (mfile,'(i2,2a)') i  ,'       ',descr(i)
cc        enddo
cc        write (mfile,'(/,a)') 'iqty       options     window title'
cc
cc        do i=1,size(variables)
cc          if (variables(i) > 0) then
cc            write (mfile,'(i2,4a)') variables(i)
cc     .         ,'         ',options
cc     .         ,'         ',trim(descr(variables(i)))
cc          elseif (variables(i) == 0) then
cc            cycle
cc          else
cc            variables(i) = - variables(i)
cc            label = 'Vector plot (' // trim(descr(variables(i)))
cc     .              // ',' // trim(descr(variables(i)+1))//')'
cc            write (mfile,'(i2,a,i2,4a)') variables(i),' ',variables(i)+1
cc     .         ,'       ',options,'         ',trim(label)
cc          endif
cc        enddo
cc        write (mfile,*)
cc
cc        close (mfile)
cc
ccc     End program
cc
cc      end subroutine createDrawInMfile
cc
ccc     createDrawInCfile
ccc     ################################################################
cc      subroutine createDrawInCfile(nvar,binf,graphlabel,ivart,ivarx
cc     .                       ,ivary,descr,options,drawfile)
cc
ccc     ----------------------------------------------------------------
ccc     Creates drawfile for "C-type" xdraw graphics.
ccc     ----------------------------------------------------------------
cc
cc        implicit none
cc
ccc     Call variables
cc
cc        integer     :: nvar
cc        character(*):: binf,drawfile,graphlabel,options,descr(nvar)
cc        character(*):: ivart,ivarx,ivary
cc
ccc     Local variables
cc
cc        integer      :: i,j,cfile
cc        character(50):: label
cc        logical      :: flag
cc
ccc     Begin program
cc
cc        cfile = 111
cc
cc        open (unit=cfile,file=drawfile,status='unknown')
cc
cc        write (cfile,'(1a,/)') 
cc     .     'Type (G=Graph, C=Contour, M=TS contour):      C'
cc        write (cfile,'(a)') 'filename(s)'
cc        write (cfile,'(a,/)') trim(binf)
cc
cc        write (cfile,'(2a,/)') 'comment:      ',graphlabel
cc        write (cfile,'(a)') 'variable names'
cc        write (cfile,'(a)')    '0        H'
cc        write (cfile,'(2a)')   'i0       ',trim(ivart)
cc        write (cfile,'(2a)')   'i1       ',trim(ivarx)
cc        write (cfile,'(2a,/)') 'i2       ',trim(ivary)
cc        write (cfile,'(a)') 'ix      iy      other           title'
cc
cc        do i=1,nvar
cc          write (cfile,'(4a,i2,2a)') 'i2      i1'
cc     .         ,'     ',options,' -i0=',i-1
cc     .         ,'     ',trim(descr(i))
cc        enddo
cc        write (cfile,*)
cc
cc        close (cfile)
cc
ccc     End program
cc
cc      end subroutine createDrawInCfile

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
     .         call contour_step(igroup,ugraph(igroup),nqty(igroup)
     .                          ,graph(igroup)%array_graph
     .                          ,graph(igroup)%cartesian,2)

          close(unit=ugraph(igroup))
          close(unit=uprofile(igroup))
        enddo

        if (debug) close(unit=udebug)

      end subroutine closeXDRAWGraphicsFiles

      end module xdraw_graphics
