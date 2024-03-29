
c ######################################################################
      module xdraw_io

      use graphics_io

      use io

      contains

c     contour
c     #####################################################################
      subroutine contour(arr,nx,ny,xmin,xmax,ymin,ymax,iopt,nunit)

c     ---------------------------------------------------------------------
c     Contours 2D array in xdraw format "C". In call:
c       * arr: 2D array to be plotted
c       * nx,ny: dimensions of array
c       * xmin,xmax,ymin,ymax: 2D domain limits
c       * iopt: whether to initialize xdraw plot (iopt=0) or not.
c       * nunit: integer file identifier.
c     ---------------------------------------------------------------------

      implicit none               !For safe fortran

c     Call variables

      integer    :: nx,ny,iopt,nunit
      real(8)    :: arr(nx,ny),xmin,xmax,ymin,ymax

c     Local variables

      integer    :: i,j

c     Begin program

      if(iopt == 0) then
        write(nunit) nx-1,ny-1,0
        write(nunit) real(xmin,4),real(xmax,4)
     .              ,real(ymin,4),real(ymax,4) 
      endif
      write(nunit) ((real(arr(i,j),4),i=1,nx),j=1,ny)

c     End program

      end subroutine contour

c     C_contour
c     #####################################################################
      subroutine C_contour(array,xmn,xmx,ymn,ymx,zmn,zmx,iopt,nunit)

c     ---------------------------------------------------------------------
c     Dumps graph arrays in array_graph in C format, but with arbitrary
c     cuts specified in cont_conf. On input:
c       * array: array to be plotted
c       * xmn,xmx,ymn,ymx,zmn,zmx: domain dimensions
c       * iopt:   = 0 --> Initialize file
c                 > 0 --> Append data
c       * nunit: file unit identification
c     ---------------------------------------------------------------------

        implicit none

c     Call variables

        integer :: iopt,nunit

        real(8) :: xmn,xmx,ymn,ymx,zmn,zmx
        real(8),dimension(:,:,:) :: array

c     Local variables

        integer :: i,j,k,imn,imx,jmn,jmx,kmn,kmx,plt(3)
        logical :: cut(3)

c     Begin program

        if (cont_conf%plane == 0) then
          write (*,*) '--> Plotting plane undefined in C_contour'
          write (*,*) '--> Aborting plot...'
          return
        endif

c     Select 2D cuts

        plt = cont_conf%coords

        cut = .false.
        cut(cont_conf%plane) = .true.

        imn = 1
        imx = size(array,1)

        jmn = 1
        jmx = size(array,2)

        kmn = 1
        kmx = size(array,3)

c     Perform 2d cuts

        if (cut(1)) then
          call twodcut(plt(1),plt(1),plt(1),jmn,jmx,kmn,kmx
     .                ,'XDRAW X-slice contour at i=')
        elseif (cut(2)) then
          call twodcut(plt(2),imn,imx,plt(2),plt(2),kmn,kmx
     .                ,'XDRAW Y-slice contour at j=')
        elseif (cut(3)) then
          call twodcut(plt(3),imn,imx,jmn,jmx,plt(3),plt(3)
     .                ,'XDRAW Z-slice contour at k=')
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

        logical    :: flag

c     Begin program

cc        write (*,*) messg,iplt

        if (iopt == 0) then

cc          write(nunit) (imxg-imng),(jmxg-jmng),(kmxg-kmng)
          if (imng == imxg) then
            write(nunit) (jmxg-jmng),(kmxg-kmng),0
            write(nunit) real(ymn,4),real(ymx,4)
     .                  ,real(zmn,4),real(zmx,4) 
          elseif (jmng == jmxg) then
            write(nunit) (imxg-imng),(kmxg-kmng),0
            write(nunit) real(xmn,4),real(xmx,4)
     .                  ,real(zmn,4),real(zmx,4) 
          else
            write(nunit) (imxg-imng),(jmxg-jmng),0
            write(nunit) real(xmn,4),real(xmx,4)
     .                  ,real(ymn,4),real(ymx,4) 
          endif

        endif

        write(nunit) (((real(array(i,j,k),4),i=imng,imxg)
     .                                      ,j=jmng,jmxg)
     .                                      ,k=kmng,kmxg)

      end subroutine twodcut

      end subroutine C_contour

c     M_contour
c     #####################################################################
      subroutine M_contour(nunit,nqty,array_graph,iopt,conf,car,xi
     $                    ,imax_g,imin_g,jmax_g,jmin_g,kmax_g,kmin_g)

c     ---------------------------------------------------------------------
c     Dumps graph arrays in array_graph in M and M2 xdraw format. It does
c     M2 plots (i.e, with coordinates) if car is present, M otherwise.
c     On input:
c       * nunit: file unit identification
c       * nqty: number of quantities to be plotted (contained in array_graph)
c       * array_graph: array of graphics variables of dimension nqty
c       * iopt:   = 0 --> Initialize file
c             (M) > 0 --> Append time stepping data
c             (M2)= 1 --> Append data
c             (M2)= 2 --> Finalize file
c       * car (optional): cartesian coordinates (for M2)
c     ---------------------------------------------------------------------

        implicit none

c     Call variables

        integer :: nqty,iopt,nunit!,igroup

        type (graph_var_def),dimension(nqty) :: array_graph

        real(8),optional,target,dimension(:,:,:,:) :: car,xi
        integer,optional :: imax_g,imin_g,jmax_g,jmin_g,kmax_g,kmin_g

        type(cont_info) :: conf
        
c     Local variables

        integer :: nsblk,nusblk
     $            ,imax_sv,imin_sv,jmax_sv,jmin_sv,kmax_sv,kmin_sv
        integer :: i,j,k,igr,plt(3)
        logical :: cut(3),M1,M2

c     Begin program

        nsblk  = 1              ! Number of structured blocks
        nusblk = 0              ! Number of unstructured blocks

        M1 = PRESENT(xi)
        M2 = PRESENT(car)
        
        imax_sv = imaxg
        imin_sv = iming
        jmax_sv = jmaxg
        jmin_sv = jming
        kmax_sv = kmaxg
        kmin_sv = kming

        if (PRESENT(imax_g)) imax_sv = imax_g
        if (PRESENT(imin_g)) imin_sv = imin_g
        if (PRESENT(jmax_g)) jmax_sv = jmax_g
        if (PRESENT(jmin_g)) jmin_sv = jmin_g
        if (PRESENT(kmax_g)) kmax_sv = kmax_g
        if (PRESENT(kmin_g)) kmin_sv = kmin_g

c     Select 2D cuts

        plt = conf%coords

        cut = .false.
        cut(conf%plane) = .true.

c     Perform 2d cuts

        if (cut(1)) then
          call twodcut(plt(1)
     $                ,plt(1),plt(1),jmin_sv,jmax_sv,kmin_sv,kmax_sv
     .                ,'XDRAW X-slice contour at i=')
        elseif (cut(2)) then
          call twodcut(plt(2)
     $                ,imin_sv,imax_sv,plt(2),plt(2),kmin_sv,kmax_sv
     .                ,'XDRAW Y-slice contour at j=')
        else
          call twodcut(plt(3)
     $                ,imin_sv,imax_sv,jmin_sv,jmax_sv,plt(3),plt(3)
     .                ,'XDRAW Z-slice contour at k=')
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

        integer :: nxgr,nygr,nxgr_old=0,nygr_old=0,zero=0,one=1
     .            ,nt=1,nt_prev=0

        logical,save :: flag

c     Begin program

c     Write messg to standard output

        write (*,*) messg,iplt

        if (M2) then   !Type M2 contours (physical space)

          select case(iopt)
          case(0)   !Write initialization headings

            write(nunit) nsblk,nusblk,nqty

            flag = .true.

          case(1)   !Write data

            if (imng == imxg) then
              nxgr = jmxg-jmng
              nygr = kmxg-kmng
              call find_cart_plane(1,i1,i2)
            elseif (jmng == jmxg) then
              nxgr = imxg-imng
              nygr = kmxg-kmng
              call find_cart_plane(2,i1,i2)
            else
              nxgr = imxg-imng
              nygr = jmxg-jmng
              call find_cart_plane(3,i1,i2)
            endif

            if (flag) then
              write(nunit) zero,nxgr,nygr,nt,zero,zero,zero
              flag = .false.
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

          case(2)  !Finalize

            write(nunit) zero,zero,zero,zero,nxgr_old,nygr_old,nt_prev

          case default

            write(*,*) 'Wrong value of iopt in XDRAW M_contour, iopt='
     .                ,iopt
            stop

          end select

        else   !Type M contours (logical space)

          select case(iopt)
          case(0)  !Write initialization headings

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

          case(1)   !Write data

            do igr=1,nqty
              write(nunit) (((real(array_graph(igr)%array(i,j,k),4)
     .                    ,i=imng,imxg),j=jmng,jmxg),k=kmng,kmxg)
            enddo

          case default

            write(*,*) 'Wrong value of iopt in XDRAW M_contour, iopt='
     .                 ,iopt
            stop

          end select

        endif

      end subroutine twodcut

c     write_coords
c     ################################################################
      subroutine write_coords(imng,imxg,jmng,jmxg,kmng,kmxg,i1,i2)

        implicit none

c     Call variables

        integer :: imng,imxg,jmng,jmxg,kmng,kmxg,i1,i2

c     Local variables

        integer :: ii,im,jm,km

        real(8),pointer,dimension(:,:,:,:) :: lcar,lxi

c     Begin program

        if (M2) then      !Write cartesian coordinates

          lcar => car     !lcar starts at i=j=k=1, not 0!!

          write(nunit) (((real(lcar(i+1,j+1,k+1,i1),4),i=imng,imxg)
     .                                                ,j=jmng,jmxg)
     .                                                ,k=kmng,kmxg)
     .                ,(((real(lcar(i+1,j+1,k+1,i2),4),i=imng,imxg)
     .                                                ,j=jmng,jmxg)
     .                                                ,k=kmng,kmxg)

          nullify(lcar)

        else             !Write logical coordinates

          lxi => xi      !lxi starts at i=j=k=1, not 0!!

          im = imng-1    !To make starting point @ i=1
          jm = jmng-1    !To make starting point @ j=1
          km = kmng-1    !To make starting point @ k=1

          write(nunit) (((real(lxi(i-im,j-jm,k-km,i1),4),i=imng,imxg)
     .                                                  ,j=jmng,jmxg)
     .                                                  ,k=kmng,kmxg)
     .                ,(((real(lxi(i-im,j-jm,k-km,i2),4),i=imng,imxg)
     .                                                  ,j=jmng,jmxg)
     .                                                  ,k=kmng,kmxg)
          nullify(lxi)
          
        endif

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

        integer :: iplane,i1,i2

c     Local variables

        integer :: ii
        real(8) :: max_x,max_y,max_z,min_x,min_y,min_z

c     Begin program

cc        select case(iplane)
cc        case(1)
cc          i1 = 2
cc          i2 = 3
cc        case(2)
cc          i1 = 1
cc          i2 = 3
cc        case(3)
cc          i1 = 1
cc          i2 = 2
cc        end select
cc
cc        return

        ii = 2 !This needs to be "2" because car starts at 1; BROKEN logic; needs fixing

        if (iplane == 1) then   
          max_x = maxval(car(ii,:,:,1))
          max_y = maxval(car(ii,:,:,2))
          max_z = maxval(car(ii,:,:,3))

          min_x = minval(car(ii,:,:,1))
          min_y = minval(car(ii,:,:,2))
          min_z = minval(car(ii,:,:,3))
        elseif (iplane == 2) then   !Same for j
          max_x = maxval(car(:,ii,:,1))
          max_y = maxval(car(:,ii,:,2))
          max_z = maxval(car(:,ii,:,3))

          min_x = minval(car(:,ii,:,1))
          min_y = minval(car(:,ii,:,2))
          min_z = minval(car(:,ii,:,3))
        elseif (iplane == 3) then   !Same for k
          max_x = maxval(car(:,:,ii,1))
          max_y = maxval(car(:,:,ii,2))
          max_z = maxval(car(:,:,ii,3))

          min_x = minval(car(:,:,ii,1))
          min_y = minval(car(:,:,ii,2))
          min_z = minval(car(:,:,ii,3))
        endif

        i1 = 0
        i2 = 0

cc        if (max_x-min_x<1d-14) then
cc          i1 = 2
cc          i2 = 3
cc        endif
cc
cc        if (max_y-min_y<1d-14) then
cc          i1 = 1
cc          i2 = 3
cc        endif
cc
cc        if (max_z-min_z<1d-14) then
cc          i1 = 1
cc          i2 = 2
cc        endif
        if (max_x-min_x<max_y-min_y .and.max_x-min_x<max_z-min_z) then
          i1 = 2
          i2 = 3
        endif

        if (max_y-min_y<max_x-min_x .and.max_y-min_y<max_z-min_z) then
          i1 = 1
          i2 = 3
        endif

        if (max_z-min_z<max_x-min_x .and.max_z-min_z<max_y-min_y) then
          i1 = 1
          i2 = 2
        endif

        if (i1 == 0 .or. i2 == 0) then
          write(*,*) 'Could not find cartesian plane in find_cart_plane'
          write(*,*) 'Plane=',iplane
          write(*,*) 'Min X, Max X=',min_x,max_x
          write(*,*) 'Min Y, Max Y=',min_y,max_y
          write(*,*) 'Min Z, Max Z=',min_z,max_z
          write(*,*) 'Aborting XDRAW ...'
          stop
        endif

      end subroutine find_cart_plane

      end subroutine M_contour

ccc     createDrawInGfile
ccc     ################################################################
cc      subroutine createDrawInGfile(nvar,ntot,binfile,graphlabel
cc     .                       ,ivar_name,dvar,dvar_name,ivar,logplot
cc     .                       ,drawfile,connect,spline,title)
cc
ccc     ----------------------------------------------------------------
ccc     Creates the G-type (x-y curves) draw.in file for xdraw graphics.
ccc     On call sequence:
ccc       * nvar       : number of variables in consideration for plotting
ccc       * ntot       : total number of variables
ccc       * binfile    : name of binary file to be read (string)
ccc       * graphlabel : character array of arbitrary length specifying
ccc                      graph title
ccc       * ivar_name  : name of independent variable (string)
ccc       * dvar       : integer array of size (nvar) identifying the 
ccc                      dependent variables
ccc       * dvar_name  : character array containing dependent variables
ccc                      descriptions.
ccc       * ivar       : integer array specifying independent variables 
ccc                      for all the dependent variables.
ccc                      These include any of the dependent
ccc                      variable indeces plus 0 (indicating the true
ccc                      independent variable).
ccc       * logplot    : integer array specifying which (in)dependent
ccc                      variables will be plotted in log scale.
ccc       * drawfile   : character variable specifying the name of the
ccc                      draw.in file.
ccc       * connect    : logical variable that specifies whether 
ccc                      xdraw should connect points or not.
ccc       * spline     : logical variable that specifies whether 
ccc                      xdraw should do spline interpolations or not.
ccc       * title (optional): character array containing title of
ccc                      xdraw windows (if not present, dvar_name is used)
ccc     ----------------------------------------------------------------
cc
cc        implicit none
cc
ccc     Call variables
cc
cc        integer,intent(in) :: nvar,ntot,ivar(ntot),dvar(nvar)
cc     .                       ,logplot(0:ntot)
cc        character(*),intent(in) :: binfile,drawfile,graphlabel,ivar_name
cc        character(*),intent(in) :: dvar_name(ntot)
cc        logical      :: spline,connect
cc        character(*),optional :: title(ntot)
cc
ccc     Local variables
cc
cc        integer      :: i,j,gfile,lngt
cc        character(50):: label,dscr(ntot)
cc        character(2) :: options
cc        character(6) :: str
cc        logical      :: flag
cc
ccc     Begin program
cc
cc        if (nvar > ntot) then
cc          write (*,*) 'Error in createDrawInGfile'
cc          write (*,*) 'More plotting variables than available'
cc          write (*,*) 'Aborting...'
cc          stop
cc        endif
cc
cc        gfile = 111
cc
cc        open (unit=gfile,file=drawfile,status='unknown')
cc
cc        write (gfile,'(1a,/)') 
cc     .     'Type (G=Graph, C=Contour, M=TS contour):      G'
cc        write (gfile,'(a)') 'filename(s)'
cc        write (gfile,'(a,/)') trim(binfile)
cc
cc        write (gfile,'(2a,/)') 'graph title: ',trim(graphlabel)
cc        write (gfile,'(a)') 'variable names'
cc        write (gfile,'(2a)') ' 0       ',trim(ivar_name)
cc
cc        do i=1,size(dvar_name)
cc          lngt = len(trim(dvar_name(i)))
cc          if (lngt == 0) cycle
cc          write (gfile,'(i2,2a)') i  ,'       ',trim(dvar_name(i))
cc        enddo
cc        write (gfile,'(/,a)') 'ix       iy        window title'
cc
cc        options = ''
cc        if (spline) options = trim(options) // '&'
cc        if (.not.connect) options = trim(options) // '.'
cc
cc        if (PRESENT(title)) then
cc          dscr = title
cc        else
cc          dscr = dvar_name
cc        endif
cc
cc        if (logplot(0) == 1) then
cc          str = 'L     '
cc        else
cc          str = '      '
cc        endif
cc
cc        do i=1,size(dscr)
cc          lngt = len(trim(dscr(i)))
cc          if (lngt == 0) exit
cc          flag = .false.
cc          do j = 1,nvar
cc            if (dvar(j).eq.i) then
cc              flag = .true.
cc              exit
cc            endif
cc          enddo
cc          if (flag) then
cc            if(logplot(i) == 1) then
cc              write (gfile,'(i2,a,i2,2a)')  ivar(i),str
cc     .                                ,i,'L'//trim(options)//'   '
cc     .                                ,trim(dscr(i))
cc            else
cc              write (gfile,'(i2,a,i2,2a)')  ivar(i),str
cc     .                                ,i,trim(options)//'      '
cc     .                                ,trim(dscr(i))
cc            endif
cc          else
cc            if(logplot(i) == 1) then
cc              write (gfile,'(a,i2,a,i2,2a)') ';',ivar(i),str
cc     .                                ,i,'L'//trim(options)//'   '
cc     .                                ,trim(dscr(i))
cc            else
cc              write (gfile,'(a,i2,a,i2,2a)') ';',ivar(i),str
cc     .                                ,i,trim(options)//'      '
cc     .                                ,trim(dscr(i))
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

c     createDrawInGfile
c     ################################################################
      subroutine createDrawInGfile(nvar,ntot,binfile,graphlabel
     .                       ,ivar_name,dvar,dvar_name,ivar,logplot
     .                       ,drawfile,connect,spline,title)

c     ----------------------------------------------------------------
c     Creates the G-type (x-y curves) draw.in file for xdraw graphics.
c     On call sequence:
c       * nvar       : number of variables in consideration for plotting
c       * ntot       : total number of variables
c       * binfile    : name of binary file to be read (string)
c       * graphlabel : character array of arbitrary length specifying
c                      graph title
c       * ivar_name  : name of independent variable (string)
c       * dvar       : integer array of size (nvar) identifying the 
c                      dependent variables
c       * dvar_name  : character array containing dependent variables
c                      descriptions.
c       * ivar       : integer array specifying independent variables 
c                      for all the dependent variables.
c                      These include any of the dependent
c                      variable indeces plus 0 (indicating the true
c                      independent variable).
c       * logplot    : integer array specifying which (in)dependent
c                      variables will be plotted in log scale.
c       * drawfile   : character variable specifying the name of the
c                      draw.in file.
c       * connect    : logical variable that specifies whether 
c                      xdraw should connect points or not.
c       * spline     : logical variable that specifies whether 
c                      xdraw should do spline interpolations or not.
c       * title (optional): character array containing title of
c                      xdraw windows (if not present, dvar_name is used)
c     ----------------------------------------------------------------

        implicit none

c     Call variables

        integer,intent(in) :: nvar,ntot,ivar(ntot),dvar(nvar)
     .                       ,logplot(0:ntot)
        character(*),intent(in) :: binfile(:)
     .                            ,drawfile,graphlabel,ivar_name
        character(*),intent(in) :: dvar_name(ntot)
        logical      :: spline,connect
        character(*),optional :: title(ntot)

c     Local variables

        integer      :: i,j,gfile,lngt,nttt
        character(50):: label,dscr(ntot)
        character(2) :: options
        character(6) :: str,ev,nnttt
        logical      :: flag

c     Begin program

        if (nvar > ntot) then
          write (*,*) 'Error in createDrawInGfile'
          write (*,*) 'More plotting variables than available'
          write (*,*) 'Aborting...'
          stop
        endif

        gfile = 111

        open (unit=gfile,file=drawfile,status='unknown')

        write (gfile,'(1a,/)') 
     .     'Type (G=Graph, C=Contour, M=TS contour):      G'
        write (gfile,'(a)') 'filename(s)'

        do i=1,size(binfile)
	 write(gfile,'(a)') trim(binfile(i))
	enddo

        write (gfile,'(/)') 

        write (gfile,'(2a,/)') 'graph title: ',trim(graphlabel)
        write (gfile,'(a)') 'variable names'
        write (gfile,'(2a)') ' 0       ',trim(ivar_name)

        do i=1,size(dvar_name)
          lngt = len(trim(dvar_name(i)))
          if (lngt == 0) cycle
          write (gfile,'(i2,2a)') i  ,'       ',trim(dvar_name(i))
        enddo
        write (gfile,'(/,a)') 'ix       iy        window title'

        options = ''
        if (spline) options = trim(options) // '&'
        if (.not.connect) options = trim(options) // '.'

        if (PRESENT(title)) then
          dscr = title
        else
          dscr = dvar_name
        endif

        if (logplot(0) == 1) then
          str = 'L     '
        else
          str = '      '
        endif

        do i=1,size(dscr)
          lngt = len(trim(dscr(i)))
          if (lngt == 0) exit
          flag = .false.
          do j = 1,nvar
            if (dvar(j).eq.i) then
              flag = .true.
              exit
            endif
          enddo
          if (flag) then
            if(logplot(i) == 1) then
              write (gfile,'(i2,a,i2,2a)')  ivar(i),str
     .                                ,i,'L'//trim(options)//'   '
     .                                ,trim(dscr(i))
            else
              write (gfile,'(i2,a,i2,2a)')  ivar(i),str
     .                                ,i,trim(options)//'      '
     .                                ,trim(dscr(i))
            endif
          else
            if(logplot(i) == 1) then
              write (gfile,'(a,i2,a,i2,2a)') ';',ivar(i),str
     .                                ,i,'L'//trim(options)//'   '
     .                                ,trim(dscr(i))
            else
              write (gfile,'(a,i2,a,i2,2a)') ';',ivar(i),str
     .                                ,i,trim(options)//'      '
     .                                ,trim(dscr(i))
            endif
          endif
        enddo

        write (gfile,*)

        close (gfile)

c     End program

      end subroutine createDrawInGfile

c     createDrawInMfile
c     ################################################################
      subroutine createDrawInMfile(nvar,binf,cartesian,graphlabel
     .            ,ivart,ivarx,ivary,variables,descr,options,drawfile)

c     ----------------------------------------------------------------
c     Creates drawfile for "M-type" (dynamic contour) xdraw graphics.
c     On call sequence:
c       * nvar       : number of variables to plot
c       * binf       : name of binary file to be read
c       * cartesian  : whether we are using logical or physical coord.
c       * graphlabel : character array of arbitrary length specifying
c                      graph title
c       * ivart      : time variable axis name
c       * ivarx      : X variable axis name
c       * ivary      : Y variable axis name
c       * variables  : variable specification (integer; negative
c                        integer means vector plot, using consecutive
c                        variable)
c       * descr      : character array containing dependent variables
c                      descriptions.
c       * options    : string with contour options
c       * drawfile   : string specifying the name of the draw.in file.
c     ----------------------------------------------------------------

        implicit none

c     Call variables

        integer     :: nvar,variables(:)
        character(*):: binf,drawfile,graphlabel,options,descr(:)
        character(*):: ivart,ivarx,ivary
        logical     :: cartesian

c     Local variables

        integer      :: i,j,mfile,lngt,lim
        character(50):: label,opt
        logical      :: flag

c     Begin program

        mfile = 111

        open (unit=mfile,file=drawfile,status='unknown')

        if (cartesian) then
          write (mfile,'(1a,/)') 
     .     'Type (G=Graph, C=Contour, M=TS contour):      M2'
          opt = options // ' -a'
        else
          write (mfile,'(1a,/)') 
     .     'Type (G=Graph, C=Contour, M=TS contour):      M'
          opt = options
        endif
        write (mfile,'(a)') 'filename(s)'
        write (mfile,'(a,/)') trim(binf)

        write (mfile,'(2a,/)') 'comment:      ',graphlabel
        write (mfile,'(a)') 'independent variable names'
        write (mfile,'(2a)')   't       ',trim(ivart)
        write (mfile,'(2a)')   'x       ',trim(ivarx)
        write (mfile,'(2a,/)') 'y       ',trim(ivary)
        write (mfile,'(a)') 'dependent variable names'
        do i=1,size(descr)
          lngt = len(trim(descr(i)))
          if (lngt == 0) exit
          write (mfile,'(i2,2a)') i  ,'       ',descr(i)
        enddo
        write (mfile,'(/,a)') 'iqty       options     window title'

        lim = min(i-1,size(variables))
        do i=1,lim
          if (variables(i) > 0) then
            write (mfile,'(i2,4a)') variables(i)
     .         ,'         ',trim(opt)
     .         ,'         ',trim(descr(variables(i)))
          elseif (variables(i) == 0) then
            cycle
          else
            variables(i) = - variables(i)
            label = 'Vector plot (' // trim(descr(variables(i)))
     .              // ',' // trim(descr(variables(i)+1))//')'
            write (mfile,'(i2,a,i2,4a)') variables(i),' ',variables(i)+1
     .         ,'       ',trim(opt),'         ',trim(label)
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
c     Creates drawfile for "C-type" (fixed contour type) xdraw graphics.
c     On call sequence:
c       * nvar : number of variables in consideration for plotting
c       * binf : name of binary file to be read by XDRAW
c       * graphlabel : character array of arbitrary length specifying
c                      graph title
c       * ivart      : time variable axis name
c       * ivarx      : X variable axis name
c       * ivary      : Y variable axis name
c       * descr      : character array containing dependent variables
c                      descriptions.
c       * options    : string with contour options
c       * drawfile   : string specifying the name of the draw.in file.
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
        write (cfile,'(2a)')   'i1       ',trim(ivary)
        write (cfile,'(2a,/)') 'i2       ',trim(ivarx)
        write (cfile,'(a)') 'ix      iy      other           title'

        do i=1,nvar
          if (i<=10) then
            write (cfile,'(4a,i1,2a)') 'i2      i1'
     .           ,'     ',options,' -i0=',i-1
     .           ,'     ',trim(descr(i))
          else
            write (cfile,'(4a,i2,2a)') 'i2      i1'
     .           ,'     ',options,' -i0=',i-1
     .           ,'     ',trim(descr(i))
          endif
        enddo
        write (cfile,*)

        close (cfile)

c     End program

      end subroutine createDrawInCfile

      end module xdraw_io
