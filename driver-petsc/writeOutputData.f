c writeOutputData
c ######################################################################
      subroutine writeOutputData(array,imin,imax,jmin,jmax,kmin,kmax
     .                             ,gmits,nwits)

      use variables

      use timeStepping

      use iosetup

      use counters

      implicit none

c Call variables

      integer(4)      :: imin,imax,jmin,jmax,kmin,kmax,nwits,gmits

      type(petsc_var) :: array(imin:imax,jmin:jmax,kmin:kmax)

c Local variables

      type(petsc_array) :: petscarray

      integer(4)        :: ierr,i,j,k,ieq

      real(8) :: mag

c Begin program

      ierr = 0

      call allocatePetscType(petscarray)

      petscarray%array(imin:imax,jmin:jmax,kmin:kmax)
     .         = array(imin:imax,jmin:jmax,kmin:kmax)

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

      call output


c Deallocate memory

      call deallocatePetscType(petscarray)

      end subroutine writeOutputData
