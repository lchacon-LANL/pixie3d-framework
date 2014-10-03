c module graphics_io
c ######################################################################
      module graphics_io

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

        logical    :: plot,skip_diag=.false.

#if defined(petsc)
        logical    :: hdf_plot=.true.,xdraw_plot=.true.
#else
        logical    :: hdf_plot=.false.,xdraw_plot=.true.
#endif

#if defined(adios)
        logical    :: adios_plot=.true.
#else
        logical    :: adios_plot=.false.
#endif

        integer    :: iming,imaxg,jming,jmaxg,kming,kmaxg,igroup
     .               ,iggx,iggy,iggz,iig,jjg,kkg

        integer    :: sel_diag (xdraw_cont_lim)
     .               ,sel_graph(xdraw_cont_lim)

        integer    :: ndplot
        real(8)    :: dplot,tmplot

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

      integer :: i

      namelist /graphics/ sel_diag,sel_graph
     .                   ,ndplot,dplot
     .                   ,hdf_plot,xdraw_plot,adios_plot
     .                   ,prof_conf,cont_conf
     .                   ,iplot,jplot,kplot

c     Graphics defaults

      ndplot = 0     !Plotting interval (integer units)
      dplot  = 0d0   !Plotting interval (time units)

      iplot = 1
      jplot = 1
      kplot = 1

      sel_diag  =(/ (i,i=1,xdraw_cont_lim) /)
      sel_graph =(/ (i,i=1,xdraw_cont_lim) /)

c     Read graphics initialization parameters

      open(unit=25,file=trim(file),status='old')
      read(25,graphics,iostat=ioerr)
      close(unit=25)

c     End program

      end function readGraphicsConfig

      end module graphics_io
