c writeOutputData
c ######################################################################
      subroutine writeOutputData(array,imingc,imaxgc,jmingc,jmaxgc
     .                                ,kmingc,kmaxgc,gmits,nwits)

      use variables

      use timeStepping

      use iosetup

      use counters

      implicit none

c Call variables

      integer(4)      :: imingc,imaxgc,jmingc,jmaxgc,kmingc,kmaxgc
     .                   ,gmits,nwits

      type(petsc_var) ::array(imingc:imaxgc,jmingc:jmaxgc,kmingc:kmaxgc)

c Local variables

      integer(4)      :: imingcl,imaxgcl,jmingcl,jmaxgcl,kmingcl,kmaxgcl

      type(petsc_array) :: petscarray

      integer(4)        :: ierr,i,j,k,ieq

      real(8) :: mag

c Begin program

      ierr = 0

      call allocatePetscType(petscarray)

      call fromGlobalToLocalLimits(imingc ,jmingc ,kmingc
     $                            ,imingcl,jmingcl,kmingcl)
      call fromGlobalToLocalLimits(imaxgc ,jmaxgc ,kmaxgc
     $                            ,imaxgcl,jmaxgcl,kmaxgcl)

      petscarray%array(imingcl:imaxgcl,jmingcl:jmaxgcl,kmingcl:kmaxgcl)
     .         = array(imingc :imaxgc ,jmingc :jmaxgc ,kmingc :kmaxgc )

c Map old solution

      u_n = petscarray

c Time level plots (xdraw)

      if (nrst.eq.ndstep.or.tmrst.ge.dstep) then
        nrst  = 0
        if (itime.gt.0) tmrst = tmrst - dstep
        call writeRecordFile(itime,time,dt,u_n)
cc        write (*,*) itime,time
      endif

c Output per time step

      itgmres = gmits
      itnewt  = nwits

      if (my_rank == 0) call output

c Deallocate memory

      call deallocatePetscType(petscarray)

      end subroutine writeOutputData
