      module F_call_C
#include "cprecision.h"
#if defined(gfortran)
      use iso_c_binding

      interface
         subroutine DEVICE_START_TIMING() 
     .        bind(c, name='device_start_timing')
         use iso_c_binding
         implicit none    
         end subroutine DEVICE_START_TIMING
      end interface

      interface
         subroutine DEVICE_STOP_TIMING() 
     .        bind(c, name='device_stop_timing')
         use iso_c_binding
         implicit none    
         end subroutine DEVICE_STOP_TIMING
      end interface

      interface
         subroutine DEVICE_PRINT_MSG() bind(c, name='device_print_msg')
         use iso_c_binding
         implicit none    
         end subroutine DEVICE_PRINT_MSG
      end interface

      interface
         subroutine DEVICE_FINALIZE() bind(c, name='device_finalize')
         use iso_c_binding
         implicit none    
         end subroutine DEVICE_FINALIZE
      end interface

      interface
         subroutine LCG_RAND(rx,seed)
     .        bind(c, name='LCG_rand')
         use iso_c_binding
         implicit none
         real (c_double) :: rx
         integer (c_int) :: seed
         end subroutine LCG_RAND
      end interface

      interface
         subroutine RANDN(rx,seed)
     .        bind(c, name='randn')
         use iso_c_binding
         implicit none
         real (c_double) :: rx
         integer (c_int) :: seed
         end subroutine RANDN
      end interface

      interface
         subroutine DEVICE_SETUP_PIC(nxin,xin,explicit_pic,quiet_start
     .        ,q_m,v0,vth,npcles,n_sp) 
     .        bind(c, name='device_setup_pic')
         use iso_c_binding
         implicit none
         integer (c_int),value :: nxin
         real (c_float_kind)  :: xin(*)
         logical (c_bool) :: explicit_pic,quiet_start         
         real (c_float_kind)  :: q_m(*),v0(*),vth(*)
         integer (c_int) :: npcles(*)
         integer (c_int) :: n_sp
         end subroutine DEVICE_SETUP_PIC
      end interface

      interface
         subroutine DEVICE_ALLOCATE_FIELD_MEMORY(nx) 
     .        bind(c, name='device_allocate_field_memory')
         use iso_c_binding
         implicit none
         integer (c_int), value :: nx
         end subroutine DEVICE_ALLOCATE_FIELD_MEMORY
      end interface

      interface
         subroutine DEVICE_INIT_VP_EXP(dt) 
     .        bind(c, name='device_init_vp_exp')
         use iso_c_binding
         implicit none
         real (c_float_kind), value :: dt
         end subroutine DEVICE_INIT_VP_EXP
      end interface

      interface
         subroutine DEVICE_FIND_RHO(rho,nx) 
     .        bind(c, name='device_find_rho')
         use iso_c_binding
         implicit none
         real (c_float_kind) :: rho(*)
         integer (c_int), value :: nx
         end subroutine DEVICE_FIND_RHO
      end interface

      interface
         subroutine DEVICE_FIND_CURRENT(j,nx) 
     .        bind(c, name='device_find_current')
         use iso_c_binding
         implicit none
         real (c_float_kind) :: j(*)
         integer (c_int), value :: nx
         end subroutine DEVICE_FIND_CURRENT
      end interface

      interface
         subroutine DEVICE_PERTURB_PCLES(nh1,eps,dt,rho,nx)
     .        bind(c, name='device_perturb_pcles')
         use iso_c_binding
         implicit none
         integer (c_short),value :: nh1
         real    (c_float_kind),value :: eps, dt
         real (c_float_kind) :: rho(*)
         integer (c_int), value :: nx
         end subroutine DEVICE_PERTURB_PCLES
      end interface

      interface
         subroutine DEVICE_UPDATE_TNP_PCLES()
     .        bind(c, name='device_update_tnp_pcles')
         use iso_c_binding
         implicit none
         end subroutine DEVICE_UPDATE_TNP_PCLES
      end interface

      interface
         subroutine DEVICE_SET_FIELD(sendbuf, count,itime) 
     .        bind(c, name='device_set_field')
         use iso_c_binding
         implicit none
         real (c_float_kind) :: sendbuf(*)
         integer (c_int), value :: count,itime
         end subroutine DEVICE_SET_FIELD
      end interface

      interface
         subroutine DEVICE_MOVE_PCLES(dt,rho,nx) 
     .        bind(c, name='device_move_pcles')
         use iso_c_binding
         implicit none        
         real (c_float_kind), value :: dt
         real (c_float_kind) :: rho(*)
         integer (c_int), value :: nx
         end subroutine DEVICE_MOVE_PCLES
      end interface

      interface
         subroutine DEVICE_MOVE_PCLES_ACC_CN(dt,j,nx,itime,dE,
     .        is_residual)
     .        bind(c, name='device_move_pcles_acc_cn')
         use iso_c_binding
         implicit none        
         real (c_float_kind), value :: dt
         real (c_float_kind) :: j(*),dE(*)
         integer (c_int), value :: nx,itime
         logical (c_bool), value ::is_residual
         end subroutine DEVICE_MOVE_PCLES_ACC_CN
      end interface

      interface
         subroutine DEVICE_MOMENTS(m,ke)
     .        bind(c, name='device_moments')
         use iso_c_binding
         implicit none             
         real (c_float_kind) :: m(*),ke(*)
         end subroutine DEVICE_MOMENTS
      end interface
#endif

      end module F_call_C



