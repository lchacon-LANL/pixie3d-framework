c processOldSolution
c####################################################################
      subroutine processOldSolution(x,imingc,imaxgc,jmingc,jmaxgc
     .                             ,kmingc,kmaxgc,gmits,nwits)
c--------------------------------------------------------------------
c     Calculates Fi(Uj) in equations of the form: dt Ui + Fi(Uj) = 0
c--------------------------------------------------------------------

      use variables

      use timeStepping

      use iosetup

      use counters

      implicit none

c Call variables

      integer(4)  ::  imingc,imaxgc,jmingc,jmaxgc,kmingc,kmaxgc
     .               ,gmits,nwits

      type(petsc_var) :: x(imingc:imaxgc,jmingc:jmaxgc,kmingc:kmaxgc)

c Local variables

      integer(4) :: i,j,k,il,jl,kl,ieq,ii
     .             ,imingcl,imaxgcl,jmingcl,jmaxgcl,kmingcl,kmaxgcl

      real(8)    :: mag

c Begin program

c Find local limits

      call fromGlobalToLocalLimits(imingc ,jmingc ,kmingc
     $                            ,imingcl,jmingcl,kmingcl,1,1,1)
      call fromGlobalToLocalLimits(imaxgc ,jmaxgc ,kmaxgc
     $                            ,imaxgcl,jmaxgcl,kmaxgcl,1,1,1)

c Unpack petsc array

      do ieq=1,neqd
        u_n%array_var(ieq)
     .       %array(imingcl:imaxgcl,jmingcl:jmaxgcl,kmingcl:kmaxgcl)
     .       = x(imingc:imaxgc,jmingc:jmaxgc,kmingc:kmaxgc)%var(ieq)
      enddo

c Evaluate nonlinear function Fi(Uj) at time level (n+1)

      call evaluateNonlinearFunction(u_n,fold)

c Time level plots (xdraw)

      if (nrst.eq.ndstep.or.tmrst.ge.dstep) then
        nrst  = 0
        if (itime.gt.0) tmrst = tmrst - dstep
        call writeRecordFile(urecord,itime,time,dt,u_n)
cc        if (my_rank == 0) write (*,*) itime,time
      endif

c Output per time step

      itgmres = gmits
      itnewt  = nwits

      if (my_rank == 0) call output

c End program

      end subroutine processOldSolution
