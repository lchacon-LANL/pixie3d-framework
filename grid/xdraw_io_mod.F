
c module graphics_variables
c ######################################################################
      module graphics_variables

        integer        :: ngraph,xdraw_cont_lim,ngroups
        parameter        (ngraph=30,xdraw_cont_lim=16)

        real(8), allocatable, dimension(:) :: xl,yl,zl

        type :: graph_var_def
          real(8),pointer,dimension(:,:,:) :: array
          character(20)                    :: descr
          character(20)                    :: vector_name=''
        end type graph_var_def

        type :: graph_group
          type (graph_var_def),dimension(ngraph) :: array_graph
          logical                                :: cartesian
          character(20)                          :: descr
        end type graph_group

        type (graph_group),pointer,dimension(:) :: graph

        integer   ,allocatable,dimension(:) :: nqty

        logical    :: plot

#if defined(petsc)
        logical    :: hdf_plot=.true.
#else
        logical    :: hdf_plot=.false.
#endif

        integer    :: iming,imaxg,jming,jmaxg,kming,kmaxg,igroup
     .               ,iggx,iggy,iggz,iig,jjg,kkg

        integer    :: sel_diag (xdraw_cont_lim)
     .               ,sel_graph(xdraw_cont_lim)

        integer    :: ndplot
        real(8)    :: dplot,tmplot

        integer    :: iplot=1,jplot=1,kplot=1

        !Profile conf
        type :: prof_info
          integer        :: line = 1
          character*(3)  :: label = 'x'
          integer        :: coords(3) = (/1,1,1/)
        end type prof_info

        type(prof_info),save :: prof_conf

        !Contour conf
        type :: cont_info
          integer        :: plane = 3
          character*(3)  :: label(2) = (/'x','y'/)
          integer        :: coords(3)= (/1,1,1/)
        end type cont_info

        type(cont_info),save :: cont_conf

      end module graphics_variables

c module xdraw_io
c ######################################################################
      module xdraw_io

      use graphics_variables

cc      contains

      end module xdraw_io

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

cc      use graphics_variables

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

        use graphics_variables

        implicit none

c     Call variables

        integer      :: nvar,ivariables(ngraph),dvariables(nvar)
     .                 ,logvar(ngraph)
        character(*) :: binf,drawfile,graphlabel,ivar
        character(20):: descr(ngraph)
        logical      :: spline

c     Local variables

        integer      :: i,j,gfile,lngt
        character(50):: label
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
            if (spline) then
              if(logvar(i) == 1) then
                write (gfile,'(i2,a,i2,2a)') ivariables(i),'      '
     .                                ,i,'L&        ',trim(descr(i))
              else
                write (gfile,'(i2,a,i2,2a)') ivariables(i),'      '
     .                                ,i,'&        ',trim(descr(i))
              endif
            elseif(logvar(i) == 1) then
              write (gfile,'(i2,a,i2,2a)')  ivariables(i),'      '
     .                                ,i,'L        ',trim(descr(i))
            else
              write (gfile,'(i2,a,i2,2a)')  ivariables(i),'      '
     .                                ,i,'         ',trim(descr(i))
            endif
          else
            if(logvar(i) == 1) then
              write (gfile,'(a,i2,a,i2,2a)') ';',ivariables(i),'      '
     .                                ,i,'L        ',trim(descr(i))

            else
              write (gfile,'(a,i2,a,i2,2a)') ';',ivariables(i),'      '
     .                                ,i,'         ',trim(descr(i))
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
c     Creates drawfile for "M-type" xdraw graphics.
c     ----------------------------------------------------------------

        use graphics_variables

        implicit none

c     Call variables

        integer     :: nvar,variables(nvar)
        character(*):: binf,drawfile,graphlabel,options,descr(ngraph)
        character(*):: ivart,ivarx,ivary
        logical     :: cartesian

c     Local variables

        integer      :: i,j,mfile,lngt,lim
        character(50):: label
        logical      :: flag

c     Begin program

        mfile = 111

        open (unit=mfile,file=drawfile,status='unknown')

        if (cartesian) then
          write (mfile,'(1a,/)') 
     .     'Type (G=Graph, C=Contour, M=TS contour):      M2'
        else
          write (mfile,'(1a,/)') 
     .     'Type (G=Graph, C=Contour, M=TS contour):      M'
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

        use graphics_variables

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
          write (cfile,'(4a,i2,2a)') 'i2      i1'
     .         ,'     ',options,' -i0=',i-1
     .         ,'     ',trim(descr(i))
        enddo
        write (cfile,*)

        close (cfile)

c     End program

      end subroutine createDrawInCfile
