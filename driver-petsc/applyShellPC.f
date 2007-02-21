c applyShellPC
c######################################################################
      subroutine applyShellPC(yarr,xarr,imin,imax,jmin,jmax,kmin,kmax)

c----------------------------------------------------------------------
c     Applies PC on local domain. Input vector yarr, output vector xarr.
c----------------------------------------------------------------------

      use parameters

      use grid

      use variables

      use timeStepping

      use newtongm

      use constants

      use iosetup

      use icond

      use precond_setup

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

c Setup parallel BC flags for PC

      call setASMflag(asm_PC.and.(np>1))

      iout = ilevel - 3

      if (iout >= 0) then
        write (*,*) 'Proc ',my_rank,': asm    in applyPC',asm
        write (*,*) 'Proc ',my_rank,': par_bc in applyPC',par_bc
      endif

c Call fortran PC routine

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
