c evaluateNonlinearResidual
c####################################################################
      subroutine evaluateNonlinearResidual(ntot,x,f)
c--------------------------------------------------------------------
c     Calculates nonlinear residuals. 
c--------------------------------------------------------------------

      use mg_setup

      use variables

      use timeStepping

      implicit none

c Call variables

      integer  :: ntot
      real*8      x(ntot),f(ntot)

c Local variables

      integer          :: i,j,ieq,ii

      double precision :: dudt(neqd),cnf(neqd),one_over_dt(neqd)

      type (var_array) :: varray

c Begin program

c Unpack vector x, taking boundary conditions from t=n solution

      call mapVectorToStructure(x,varray,u_n)

c Evaluate nonlinear function

      call evaluateNonlinearFunction(varray,f)

c Calculate residuals

      call defineTSParameters(cnf,one_over_dt)
      
      do j = 1,nyd
        do i = 1,nxd

          ii = neqd*(i-1) + neqd*nxd*(j-1)

          do ieq=1,neqd
            dudt(ieq) = (varray%array_var(ieq)%array(i,j)
     .                -     u_n%array_var(ieq)%array(i,j) )
     .                *one_over_dt(ieq)
          enddo

          f(ii+1:ii+neqd) = dudt(:) - (1.-cnf)*f   (ii+1:ii+neqd)
     .                              -     cnf *fold(ii+1:ii+neqd)
     .                    + fsrc(ii+1:ii+neqd)

        enddo
      enddo

      f = dx(ngrd)*dy(ngrd)*f

c End program

      call deallocateDerivedType(varray)

      end subroutine

c evaluateNonlinearFunction
c####################################################################
      subroutine  evaluateNonlinearFunction(varray,fi)
c--------------------------------------------------------------------
c     Stores evaluation of Nonlinear function Fi(Uj) in vector fi.
c     The Uj's are given in varray.
c--------------------------------------------------------------------

      use parameters

      use variable_setup

      implicit none

c Call variables

      double precision :: fi(ntotd)
      type (var_array) :: varray

c Local variables

      integer          :: i,j,ii
      double precision :: tmp(neqd)

c Begin program

c Prepare auxiliar quantities (vx,vy,jz,lapbz)

      call setupNonlinearFunction(varray)

c Store function evaluation

      do j = 1,nyd
        do i = 1,nxd
          ii = neqd*(i-1) + neqd*nxd*(j-1)
          call nonlinearRHS(i,j,varray,tmp)
          fi(ii+1:ii+neqd) = tmp
        enddo
      enddo

c End program

      end subroutine
