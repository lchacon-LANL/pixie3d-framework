c module graphics_variables
c ######################################################################
      module graphics_variables

        use xdraw_io

        integer        :: ngroups

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

        integer    :: iplot,jplot,kplot

        !Profile conf
        type :: prof_info
          integer        :: line = 1
          character(3)   :: label = 'x'
          integer        :: coords(3) = (/1,1,1/)
        end type prof_info

        type(prof_info),save :: prof_conf

        !Contour conf
        type :: cont_info
          integer        :: plane = 3
          character(3)   :: label(2) = (/'x','y'/)
          integer        :: coords(3)= (/1,1,1/)
        end type cont_info

        type(cont_info),save :: cont_conf

      end module graphics_variables
