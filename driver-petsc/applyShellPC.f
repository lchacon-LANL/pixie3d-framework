c applyShellPC
c######################################################################
      subroutine applyShellPC(yarr,xarr,imin,imax,jmin,jmax,kmin,kmax)

c----------------------------------------------------------------------
c     Initializes MG and creates grid
c----------------------------------------------------------------------

      use parameters

      use grid

      use variables

      use timeStepping

      use newtongm

      use constants

      use iosetup

      use icond

      implicit none

c Call variables

      integer(4)      :: imin,imax,jmin,jmax,kmin,kmax

      type(petsc_var) :: yarr(imin:imax,jmin:jmax,kmin:kmax)
     .                  ,xarr(imin:imax,jmin:jmax,kmin:kmax)

c Local variables

      integer(4) :: ieq,il,jl,kl,i,j,k,ii,iout
      real(8)    :: x(ntotd),y(ntotd)

c Begin program

      x = 0d0

      do k = kmin,kmax
        do j = jmin,jmax
          do i = imin,imax
            call fromGlobalToLocalLimits(i,j,k,il,jl,kl,1,1,1)
            ii = vecPos(neqd,il,jl,kl,1,1,1)
            do ieq=1,neqd
              y(ii+ieq) = yarr(i,j,k)%var(ieq)
            enddo
          enddo
        enddo
      enddo

c Call PC fortran setup routine

      iout = ilevel - 3

      call applyPreconditioner(ntotd,y,x,iout)

c Return preconditioned vector

      do k = kmin,kmax
        do j = jmin,jmax
          do i = imin,imax
            call fromGlobalToLocalLimits(i,j,k,il,jl,kl,1,1,1)
            ii = vecPos(neqd,il,jl,kl,1,1,1)
            do ieq=1,neqd
              xarr(i,j,k)%var(ieq) = x(ii+ieq)
            enddo
          enddo
        enddo
      enddo

c End program

      end subroutine applyShellPC
