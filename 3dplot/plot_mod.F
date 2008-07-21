c module icond
c ######################################################################
      module icond

        implicit none

        integer :: nh1,nh2,nh3

        logical :: odd,random

        real(8),dimension(10) :: pert

      end module icond

c module newtongm
c ######################################################################
      module newtongm

        integer    ::  maxitnwt,maxitgm,iguess,maxksp,method,global

        real(8)    ::  tolgm,rtol,atol,stol,mf_eps

      end module newtongm


c module precond_setup
c ######################################################################
      module precond_setup

        integer       :: nsweep,maxvcyc,precpass,debug_it
        real(8)       :: mgtol,omega

        character*(10):: precon

        logical       :: asm_PC

      end module precond_setup

c module counters
c ######################################################################
      module counters

        integer    :: itnewt,itgmres,itwhistler
        integer    :: gmres_tot,newt_tot,wh_tot

      end module counters

c module timeStepping
c ######################################################################
      module timeStepping

        integer    :: numtime,ndstep,nrstep,sm_pass,sm_flag,inewtime

        real(8)    :: dt,dtbase,cnfactor,cnf_d,alpha,gammat
        real(8)    :: tmax,dstep,rstep,time,dfreq,dtexp,damp,dt0

        real(8)    :: vx_max,vy_max,vz_max,bx_max,by_max,bz_max

        logical    :: timecorr,restart,source,predictor,vol_wgt
     .               ,postprocess,chk_src,eigen

        real(8)    :: tmrst

        integer    :: itime,nrst

      end module timeStepping
