c module graphics_io
c ######################################################################
      module graphics_io

        use io, ONLY: find_unit
      
        integer   :: ngroups

        integer   :: ngraph   ,xdraw_cont_lim
        parameter   (ngraph=40,xdraw_cont_lim=16)

        real(8), allocatable, dimension(:) :: xl,yl,zl
        real(8), allocatable, dimension(:,:,:,:) :: xi

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

        logical    :: skip_diag=.false.

#if defined(petsc)
        logical    :: hdf_plot=.true.,hdf_debug=.false.
     .               ,xdraw_plot=.true.
#else
        logical    :: hdf_plot=.false.,hdf_debug=.false.
     .               ,xdraw_plot=.true.
#endif

#if defined(adios)
        logical    :: adios_plot=.true.
#else
        logical    :: adios_plot=.false.
#endif

        logical    :: no_graphics_io=.false.

        integer    :: iming,imaxg,jming,jmaxg,kming,kmaxg,igroup
     .               ,iggx,iggy,iggz,iig,jjg,kkg

        integer    :: sel_diag (xdraw_cont_lim)
     .               ,sel_graph(xdraw_cont_lim)

        integer    :: ndplot,nplot
        real(8)    :: dplot,tplot

        integer    :: iplot,jplot,kplot

        !Profile conf
        type :: prof_info
          integer        :: line = 1
          character(3)   :: label = 'x  '
          integer        :: coords(3) = (/1,1,1/)
        end type prof_info

        type(prof_info),save :: prof_conf

        !Contour conf
        type :: cont_info
          integer        :: plane = 0
          character(3)   :: label(2) = (/'x  ','y  '/)
          integer        :: coords(3)= (/1,1,1/)
        end type cont_info

        type(cont_info),save :: cont_conf

      contains

c     readGraphicsConfig
c     ######################################################################
      function readGraphicsConfig(file) result(ioerr)

c     -----------------------------------------------------------------
c     Initializes grid parameters
c     -----------------------------------------------------------------

      implicit none

c     Call variables

      integer :: ioerr

      character(*) :: file

c     Local variables

      integer :: i,ufile

      namelist /graphics/ sel_diag,sel_graph
     .                   ,ndplot,dplot
     .                   ,hdf_plot,hdf_debug,xdraw_plot,adios_plot
     .                   ,prof_conf,cont_conf
     .                   ,iplot,jplot,kplot,no_graphics_io
     .                   ,skip_diag

c     Graphics defaults

      ndplot = 0     !Plotting interval (integer units)
      dplot  = 0d0   !Plotting interval (time units)

      iplot = 1
      jplot = 1
      kplot = 1

cc      sel_diag  =(/ (i,i=1,xdraw_cont_lim) /)
cc      sel_graph =(/ (i,i=1,xdraw_cont_lim) /)
      sel_diag  = 0
      sel_graph = 0

c     Read graphics initialization parameters

      ufile = find_unit(25)
      open(unit=ufile,file=trim(file),status='old')
      read(ufile,graphics,iostat=ioerr)
      close(unit=ufile)

c     Initialize if not set

      if (maxval(sel_diag) == 0) sel_diag =(/ (i,i=1,xdraw_cont_lim) /)
      if (maxval(sel_graph)== 0) sel_graph=(/ (i,i=1,xdraw_cont_lim) /)

c     End program

      end function readGraphicsConfig

      end module graphics_io
