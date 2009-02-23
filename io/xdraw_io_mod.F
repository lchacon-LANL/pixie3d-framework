c module xdraw_io
c ######################################################################
      module xdraw_io

      integer   :: ngraph,xdraw_cont_lim
      parameter   (ngraph=30,xdraw_cont_lim=16)

      contains

c     contour
c     #####################################################################
      subroutine contour(arr,nx,ny,xmin,xmax,ymin,ymax,iopt,nunit)

c     ---------------------------------------------------------------------
c     Contours 2D array in xdraw format. In call:
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
c       * logvar     : integer array specifying which dependent
c                      variables will be plotted in log scale.
c       * drawfile   : character variable specifying the name of the
c                      draw.in file.
c       * connect    : logical variable that specifies whether 
c                      xdraw should connect points or not.
c       * spline     : logical variable that specifies whether 
c                      xdraw should do spline interpolations or not.
c     ----------------------------------------------------------------

        implicit none

c     Call variables

        integer      :: nvar,ntot,ivariables(ntot),dvariables(nvar)
     .                 ,logvar(ntot)
        character(*) :: binf,drawfile,graphlabel,ivar
        character(*) :: descr(ntot)
        logical      :: spline,connect
        character(*),optional :: title(ntot)

c     Local variables

        integer      :: i,j,gfile,lngt
        character(50):: label
        character(2) :: options
        logical      :: flag

c     Begin program

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
          descr = title
        endif

        do i=1,size(descr)
          lngt = len(trim(descr(i)))
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
              write (gfile,'(i2,a,i2,2a)')  ivariables(i),'      '
     .                                ,i,'L'//trim(options)//'   '
     .                                ,trim(descr(i))
            else
              write (gfile,'(i2,a,i2,2a)')  ivariables(i),'      '
     .                                ,i,trim(options)//'      '
     .                                ,trim(descr(i))
            endif
          else
            if(logvar(i) == 1) then
              write (gfile,'(a,i2,a,i2,2a)') ';',ivariables(i),'      '
     .                                ,i,'L'//trim(options)//'   '
     .                                ,trim(descr(i))
            else
              write (gfile,'(a,i2,a,i2,2a)') ';',ivariables(i),'      '
     .                                ,i,trim(options)//'      '
     .                                ,trim(descr(i))
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
        write (cfile,'(2a)')   'i1       ',trim(ivarx)
        write (cfile,'(2a,/)') 'i2       ',trim(ivary)
        write (cfile,'(a)') 'ix      iy      other           title'

        do i=1,nvar
          write (cfile,'(4a,i1,2a)') 'i2      i1'
     .         ,'     ',options,' -i0=',i-1
     .         ,'     ',trim(descr(i))
        enddo
        write (cfile,*)

        close (cfile)

c     End program

      end subroutine createDrawInCfile

      end module xdraw_io
