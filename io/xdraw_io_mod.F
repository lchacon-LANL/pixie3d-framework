c module xdraw_io
c ######################################################################
      module xdraw_io

      use graphics_io

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

        integer    :: i,j,k,imn,imx,jmn,jmx,kmn,kmx,plt(3)
        logical    :: cut(3)

c     Begin program

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
     .                ,'X-slice contour at i=')
        elseif (cut(2)) then
          call twodcut(plt(2),imn,imx,plt(2),plt(2),kmn,kmx
     .                ,'Y-slice contour at j=')
        else
          call twodcut(plt(3),imn,imx,jmn,jmx,plt(3),plt(3)
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
      subroutine M_contour(nunit,nqty,array_graph,iopt,car)

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

        integer    :: nqty,iopt,nunit!,igroup

        type (graph_var_def),dimension(nqty) :: array_graph

        real(8),optional,target,dimension(:,:,:,:) :: car

c     Local variables

        integer    :: nsblk,nusblk
        integer    :: i,j,k,igr,plt(3)
        logical    :: cut(3),M2

c     Begin program

        nsblk  = 1              ! Number of structured blocks
        nusblk = 0              ! Number of unstructured blocks

        M2 = PRESENT(car)
          
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

        logical    :: flag

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

            write(*,*) 'Wrong value of iopt in M_contour, iopt=',iopt
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

            write(*,*) 'Wrong value of iopt in M_contour, iopt=',iopt
            stop

          end select

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

        real(8),pointer,dimension(:,:,:,:) :: lcar

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

          write(nunit) (((real(xi(i,j,k,i1),4),i=imng,imxg)
     .                                        ,j=jmng,jmxg)
     .                                        ,k=kmng,kmxg)
     .                ,(((real(xi(i,j,k,i2),4),i=imng,imxg)
     .                                        ,j=jmng,jmxg)
     .                                        ,k=kmng,kmxg)

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

        integer    :: iplane,i1,i2

c     Local variables

        integer    :: ii
        real(8)    :: max_x,max_y,max_z,min_x,min_y,min_z

c     Begin program

        ii = 2 !This needs to be 2 because car starts at 1; BROKEN logic; needs fixing

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

        if (i1 == 0 .or. i2 == 0) then
          write(*,*) 'Could not find cartesian plane in find_cart_plane'
          write(*,*) 'Aborting...'
          stop
        endif

      end subroutine find_cart_plane

      end subroutine M_contour

c     createDrawInGfile
c     ################################################################
      subroutine createDrawInGfile(nvar,ntot,binf,graphlabel,ivar
     .                       ,dvariables,descr,ivariables,logvar
     .                       ,drawfile,connect,spline,title)

c     ----------------------------------------------------------------
c     Creates the G-type (x-y curves) draw.in file for xdraw graphics.
c     On call sequence:
c       * nvar : number of variables in consideration for plotting
c       * ntot : total number of variables
c       * binf : name of binary file to be read
c       * graphlabel : character array of arbitrary length specifying
c                      graph title
c       * ivar       : independent variable axis name
c       * dvariables : integer array of size (nvar) identifying the 
c                      dependent variables
c       * descr      : character array containing dependent variables
c                      descriptions.
c       * ivariables : integer array specifying independent variables 
c                      for all the dependent variables.
c                      These include any of the dependent
c                      variable indeces plus 0 (indicating the true
c                      independent variable).
c       * logvar     : integer array specifying which (in)dependent
c                      variables will be plotted in log scale.
c       * drawfile   : character variable specifying the name of the
c                      draw.in file.
c       * connect    : logical variable that specifies whether 
c                      xdraw should connect points or not.
c       * spline     : logical variable that specifies whether 
c                      xdraw should do spline interpolations or not.
c       * title (optional): character array containing title of
c                      xdraw windows (if not present, descr is used)
c     ----------------------------------------------------------------

        implicit none

c     Call variables

        integer,intent(in) ::nvar,ntot,ivariables(ntot),dvariables(nvar)
     .                 ,logvar(0:ntot)
        character(*) :: binf,drawfile,graphlabel,ivar
        character(*),intent(in) :: descr(ntot)
        logical      :: spline,connect
        character(*),optional :: title(ntot)

c     Local variables

        integer      :: i,j,gfile,lngt
        character(50):: label,dscr(ntot)
        character(2) :: options
        character(6) :: str
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
        write (gfile,'(a,/)') trim(binf)

        write (gfile,'(2a,/)') 'graph title: ',trim(graphlabel)
        write (gfile,'(a)') 'variable names'
        write (gfile,'(2a)') ' 0       ',trim(ivar)

        do i=1,size(descr)
          lngt = len(trim(descr(i)))
          if (lngt == 0) cycle
          write (gfile,'(i2,2a)') i  ,'       ',trim(descr(i))
        enddo
        write (gfile,'(/,a)') 'ix       iy        window title'

        options = ''
        if (spline) options = trim(options) // '&'
        if (.not.connect) options = trim(options) // '.'

        if (PRESENT(title)) then
          dscr = title
        else
          dscr = descr
        endif

        if (logvar(0) == 1) then
          str = 'L     '
        else
          str = '      '
        endif

        do i=1,size(dscr)
          lngt = len(trim(dscr(i)))
          if (lngt == 0) exit
          flag = .false.
          do j = 1,nvar
            if (dvariables(j).eq.i) then
              flag = .true.
              exit
            endif
          enddo
          if (flag) then
            if(logvar(i) == 1) then
              write (gfile,'(i2,a,i2,2a)')  ivariables(i),str
     .                                ,i,'L'//trim(options)//'   '
     .                                ,trim(dscr(i))
            else
              write (gfile,'(i2,a,i2,2a)')  ivariables(i),str
     .                                ,i,trim(options)//'      '
     .                                ,trim(dscr(i))
            endif
          else
            if(logvar(i) == 1) then
              write (gfile,'(a,i2,a,i2,2a)') ';',ivariables(i),str
     .                                ,i,'L'//trim(options)//'   '
     .                                ,trim(dscr(i))
            else
              write (gfile,'(a,i2,a,i2,2a)') ';',ivariables(i),str
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

        integer     :: nvar,variables(xdraw_cont_lim)
        character(*):: binf,drawfile,graphlabel,options,descr(ngraph)
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

        lim = min(i-1,xdraw_cont_lim)
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
          write (cfile,'(4a,i2,2a)') 'i2      i1'
     .         ,'     ',options,' -i0=',i-1
     .         ,'     ',trim(descr(i))
        enddo
        write (cfile,*)

        close (cfile)

c     End program

      end subroutine createDrawInCfile

      end module xdraw_io
