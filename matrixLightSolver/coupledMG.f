c TODO list:
c
c 1) Implement direct block solver for neq > 2
c
c 2) Clean mg_setup module (including pertaining routines)
c
c 4) Eliminate hardwired boundary conditions (Dirichtlet, periodic) in
c    higher-order MG restriction and prolongation.

c cmgmeth
c####################################################################
      recursive subroutine cmgmeth(neq,ntot,y,x,matvec,smoother,igridn
     .                 ,igridmin,orderres,orderprol,vcyc,iter,omega
     .                 ,tol,guess,fdiag,out)
c-------------------------------------------------------------------
c     This subroutine solves A x = y for the vector x using a coupled
c     multigrid correction scheme. The coarse grid operators are 
c     implemented matrix-free using the 'matvec' routine, and 
c     smoothing is done in a coupled manner using 'smoother'.
c-------------------------------------------------------------------

      use mg_setup

      implicit none

c Call variables

      integer*4      neq,ntot,igridn,vcyc,iter,guess,out,igridmin
      real*8         x(ntot),y(ntot),omega,tol

      external       matvec,smoother

      integer*4      orderres,orderprol

      logical        fdiag

c Local variables

cc      real*8         diag(neq,2*ntot)
      double precision, allocatable, dimension(:,:):: diag

      real*8         xx(2*ntot),yy(2*ntot),wrk(2*ntot)
      integer*4      izero,i,j,ii,ivcyc,igrid,nn,isig

      real*8         omega_res,omega_prol,rr0,rr1,mag,mag1
      integer*4      iter_res,iter_prol,guess2,outc

      real*8         dummy(ntot),rr(ntot)

      save           diag

c Begin program

      omega_res = omega
      iter_res  = iter

      omega_prol = omega
      iter_prol  = iter

      if (out.ge.2.and.igridmin.lt.igridn) write (*,5)

      outc = out
      if (igridmin.lt.igridn) outc = out - 2

      allocate(diag(neq,2*ntot))

c Set pointers

      allocate(istart(igridn),ntotv(igridn))

      istart(igridmin-1) = 1
      do i = igridmin,igridn
        istart(i) = neq*nxvp(i-1)*nyvp(i-1) + istart (i-1)
        ntotv (i) = neq*nxvp(i)*nyvp(i)
      enddo

c Find diagonal for smoothers

      if (fdiag) then
        if (out.ge.2) write (*,*) 'Forming diagonal...'
        call find_mf_diag_neq(neq,ntot,matvec,diag,igridn,igridmin)
        if (out.ge.2) write (*,*) 'Finished!'
      endif

c Initialize local solution vector (xx) and local r.h.s. (yy)

      xx  = 0d0
      yy  = 0d0
      wrk = 0d0

      if (guess.eq.0) then
        x = 0d0
        yy(istart(igridn):istart(igridn) + ntot - 1) = y(:)
      else
        xx(istart(igridn):istart(igridn) + ntot - 1) = x(:)
        yy(istart(igridn):istart(igridn) + ntot - 1) = y(:)
      endif

c Compute initial residual and check convergence

      if (guess.eq.0) then
        rr0 = sqrt(sum(y*y))
      else
        call matvec(0,ntot,x,dummy,igridn)
        rr = y - dummy
        rr0 = sqrt(sum(rr*rr))
      endif

      if (rr0.lt.1d-16*ntot) then
        if (out.ge.1) write (*,*) 'Initial solution seems exact in MG'
        deallocate(istart,ntotv,diag)
        return
      endif

      rr1 = rr0

c Start V-cycle

      guess2 = guess

      do ivcyc = 1,vcyc

        if (outc.ge.1.and.igridmin.lt.igridn) write (*,7) ivcyc

c     Perform Vcycle recursively

        if (igridn.eq.igridmin) then
          call solve (igridn,tol,iter_res,omega_res)
          exit
        else
          call vcycle(igridn)
        endif

c     Check MG convergence

        call matvec(0,ntot,xx(istart(igridn)),dummy,igridn)

        rr   = y - dummy
        mag  = sqrt(sum(rr*rr))
        mag1 = mag/rr1

        if (out.ge.3) then
          write (*,10) mag,mag1
        elseif (out.ge.2) then
          write (*,20) mag,mag1,ivcyc
        endif

        rr1 = mag

        if (mag/rr0 < tol) exit

      enddo

c Map solution from local vector xx to external vector x

      x(:) = xx(istart(igridn):istart(igridn) + ntot - 1)

c MG convergence info

      mag1 = mag/rr0

      if (igridmin.lt.igridn) then
        if (out.eq.1) then
          write (*,20) mag,mag1,vcyc
        elseif (out.ge.2.and.vcyc.gt.1) then
          write (*,*) 
          write (*,*) 'Final MG convergence info:'
          write (*,20) mag,mag1,vcyc
        endif
      endif

c End program

      deallocate(istart,ntotv,diag)

      return

 5    format (/,' MG method output:')
 7    format (/,' MG V-cycle #:',i3)
 10   format (  ' MG residual:',1p1e12.4,'; Ratio:',1p1e12.4)
 20   format (  ' MG residual:',1p1e12.4,'; Ratio:',1p1e12.4,
     .          '; V-cycle #:'i3)

      contains

c     vcycle
c     ###################################################################
        recursive subroutine vcycle(igr)

          implicit none
          integer :: igr,isigm,isig,nn

c       Begin program

          nn    = ntotv(igr)
          isig  = istart(igr)
          isigm = istart(igr-1)

c       Relax error/solution on grid number igr/igridn (find new xx)

          call solve(igr,tol,iter_res,omega_res)

c       Evaluate residual (ie wrk = yy - A xx = yy - wrk )

          call matvec(0,nn,xx(isig),wrk(isig),igr)

          wrk(isig:isig+nn-1) = yy(isig:isig+nn-1) - wrk(isig:isig+nn-1)

c       Restrict residual( i.e. yy_c = R * yy_f = R * wrk ) down a grid

          call crestrict (neq
     .                   ,yy(isigm),ntotv(igr-1),nxvp(igr-1),nyvp(igr-1)
     .                   ,wrk(isig),ntotv(igr  ),nxvp(igr  ),nyvp(igr  )
     .                   ,orderres,igr,4d0)

          guess2 = 0

c       If on coarsest grid, solve for error, else descend a grid level

          if (igr-1.eq.igridmin) then
            call solve (igr-1,tol,iter_res,omega_res)
          else
            call vcycle(igr-1)
          endif

c       Cycle back up to grid igr updating errors (xx)
c       with fixed R.H.S. (yy)

          guess2 = 1

c       Prolong error (wrk = P * xx_1) on grid 1 to grid 2

          call cprolong (neq
     .                  ,wrk(isig),ntotv(igr  ),nxvp(igr  ),nyvp(igr  )
     .                  ,xx(isigm),ntotv(igr-1),nxvp(igr-1),nyvp(igr-1)
     .                  ,orderprol,igr-1)

c       Update existing error on grid 2 (i.e. xx_2): xx_2 = xx_2 + wrk 

          xx(isig:isig+nn-1) = xx(isig:isig+nn-1) + wrk(isig:isig+nn-1)

c       Relax updated error on igr (i.e. xx_igr)

          call solve(igr,tol,iter_prol,omega_prol)

        end subroutine vcycle

c       solve
c       ###################################################################
        subroutine solve(igr,tol,iter,omega)

          implicit none
          integer :: igr,nn,isig,iter
          double precision :: tol,omega


          nn   = ntotv(igr)
          isig = istart(igr)

          call smoother(neq,nn,yy(isig),xx(isig),matvec,diag(1,isig)
     .                 ,igr,iter,omega,tol,guess2,outc)

        end subroutine solve

      end subroutine

*deck cprolong
c######################################################################
      subroutine cprolong(neq,xf,ntotf,nxf,nyf,xc,ntotc,nxc,nyc,order
     .                   ,igc)
c----------------------------------------------------------------------
c     This is a prolongation routine for MG. 
c
c     If order = 0, it employs simple injection.
c     If order > 1, it employs spline interpolation of order "order".
c----------------------------------------------------------------------

      use mg_setup

      implicit none            ! For safe Fortran

c Call variables

      integer      neq,ntotc,nxc,nyc,ntotf,nxf,nyf,order,igc
      real*8       xc(ntotc),xf(ntotf)

c Local variables
 
      real*8       xxc(ntotc/neq),xxf(ntotf/neq)

      integer*4    ic,jc,if,jf,iic,iif,i,ieq,nntotc,nntotf

c Begin program

c Injection

      if (order.eq.0) then

      do jc = 1,nyc

        do ic = 1,nxc
          jf = 2*jc
          if = 2*ic

          do ieq = 1,neq
            iic = ieq + neq*(ic-1 + nxc*(jc-1))
            iif = ieq + neq*(if-1 + nxf*(jf-1))

            xf(iif                ) = xc(iic) 
            xf(iif - neq          ) = xc(iic) 
            xf(iif - neq*nxf      ) = xc(iic) 
            xf(iif - neq*nxf - neq) = xc(iic) 
          enddo
        enddo

      enddo

      else

c Interpolation

cc        order = 0

        nntotc = ntotc/neq
        nntotf = ntotf/neq

        do ieq = 1,neq

c     Unpack vector

          do i = 1,nntotc
            xxc(i) = xc(neq*(i-1)+ieq)
          enddo

c     Use scalar prolongation

          call prolong(xxf,nntotf,nxf,nyf
     .                ,xxc,nntotc,nxc,nyc,order,igc)

c     Repack prolonged vector

          do i=1,nntotf
            xf(neq*(i-1)+ieq)=xxf(i)
          enddo

        enddo

      endif

c End program

      return
      end subroutine

*deck crestrict
c######################################################################
      subroutine crestrict(neq,xc,ntotc,nxc,nyc,xf,ntotf,nxf,nyf,order
     .                   ,igf,volf)
c----------------------------------------------------------------------
c     This is a restriction routine for MG. 
c
c     If order = 0, it employs simple agglomeration.
c     If order > 1, it employs spline interpolation of order "order".
c----------------------------------------------------------------------

      use mg_setup

      implicit none        !For safe fortran

c Call variables

      integer*4    ntotc,nxc,nyc,ntotf,nxf,nyf,order,igf,neq

      real*8       xc(ntotc),xf(ntotf),volf

c Local variables
 
      real*8       xxc(ntotc/neq),xxf(ntotf/neq)

      integer*4    ic,jc,if,jf,iic,iif,i,ieq,nntotc,nntotf

c Begin program

c Agglomeration

      if (order.eq.0) then

        do jc = 1,nyc

          do ic = 1,nxc
            jf = 2*jc
            if = 2*ic

            do ieq = 1,neq
              iic = ieq + neq*(ic-1 + nxc*(jc-1))
              iif = ieq + neq*(if-1 + nxf*(jf-1))

              xc(iic) = (xf(iif)        +xf(iif        -neq)
     .                 + xf(iif-nxf*neq)+xf(iif-nxf*neq-neq))/4d0*volf
            enddo
          enddo

        enddo

      else

c Interpolation

cc        order = 0

        nntotc = ntotc/neq
        nntotf = ntotf/neq

        do ieq = 1,neq

c     Unpack vector

          do i = 1,nntotf
            xxf(i) = xf(neq*(i-1)+ieq)
          enddo

c     Use scalar restriction

          call restrict(xxc,nntotc,nxc,nyc
     .                 ,xxf,nntotf,nxf,nyf,order,igf,volf)

c     Repack restricted vector

          do i=1,nntotc
            xc(neq*(i-1)+ieq)=xxc(i)
          enddo

        enddo

      endif

c End program

      return
      end subroutine

*deck prolong
c######################################################################
      subroutine prolong(xf,ntotf,nxf,nyf,xc,ntotc,nxc,nyc,order,igc)
c----------------------------------------------------------------------
c     This is a prolongation routine for MG. 
c
c     If order = 0, it employs simple injection.
c     If order > 1, it employs spline interpolation of order "order".
c----------------------------------------------------------------------

      use mg_setup

      implicit none            ! For safe Fortran

c Call variables

      integer      ntotc,nxc,nyc,ntotf,nxf,nyf,order,igc
      real*8       xc(ntotc),xf(ntotf)

c Local variables
 
      real*8       uc(0:nxc+1,0:nyc+1)

      real*8       xxf,yyf,ff,xx(nxc+2),yy(nyc+2)

      integer*4    ic,if,jc,jf,ii,i,igf,iic,iif

c Interpolation

      integer*4    kx,ky,nx,ny,dim,flg
      real*8, dimension(:),allocatable:: tx,ty,work
      real*8, dimension(:,:),allocatable:: bcoef

      real*8       db2val
      external     db2val

c Begin program

c Injection

      if (order.eq.0) then

        do jc = 1,nyc
          do ic = 1,nxc
            jf = 2*jc
            if = 2*ic
            iic = ic + nxc*(jc-1)
            iif = if + nxf*(jf-1)

            xf(iif)           = xc(iic) 
            xf(iif-1)         = xc(iic) 
            xf(iif - nxf)     = xc(iic) 
            xf(iif - nxf - 1) = xc(iic) 
          enddo
        enddo

      else

c Interpolation

c     Determine current coarse grid

c     Map vectors into arrays

        xx(1:nxc+2) = xl(0:nxc+1,igc)
        yy(1:nyc+2) = yl(0:nyc+1,igc)

        do ic = 1,nxc
          do jc = 1,nyc
            ii = ic + nxc*(jc-1)
            uc(ic,jc) = xc(ii)
          enddo
        enddo

c     Impose boundary conditions on arrays

c       Top and bottom (residuals are zero at Dirichlet boundaries)

        do ic=0,nxc+1
          uc(ic,nyc+1) = 0d0
          uc(ic,0)     = 0d0
        enddo

c       Periodic BC's

        do jc=0,nyc+1
          uc(0    ,jc) = uc(nxc-1,jc)
          uc(nxc+1,jc) = uc(2  ,jc)
        enddo

c     Calculate interpolation

        flg = 0
        kx = order+1
        ky = order+1
        nx = nxc+2
        ny = nyc+2
        dim = nx*ny + max(2*kx*(nx+1),2*ky*(ny+1))

        allocate(tx(nx+kx))
        allocate(ty(nx+ky))
        allocate(work(dim))
        allocate(bcoef(nx,ny))

        call db2ink(xx,nx,yy,ny,uc,nx,kx,ky,tx,ty,bcoef,work,flg)

        do jc = 1,nyc
          do ic = 1,nxc
            jf = 2*jc
            if = 2*ic
            iif = if + nxf*(jf-1)

            igf = igc+1         !Consider next finer grid

            xxf = xl(if,igf)
            yyf = yl(jf,igf)

            xf(iif) =
     .             db2val(xxf,yyf,0,0,tx,ty,nx,ny,kx,ky,bcoef,work)

            xxf = xl(if-1,igf)
            yyf = yl(jf  ,igf)

            xf(iif-1) =
     .             db2val(xxf,yyf,0,0,tx,ty,nx,ny,kx,ky,bcoef,work)

            xxf = xl(if-1,igf)
            yyf = yl(jf-1,igf)

            xf(iif-1-nxf) =
     .             db2val(xxf,yyf,0,0,tx,ty,nx,ny,kx,ky,bcoef,work)

            xxf = xl(if  ,igf)
            yyf = yl(jf-1,igf)

            xf(iif-nxf) =
     .             db2val(xxf,yyf,0,0,tx,ty,nx,ny,kx,ky,bcoef,work)

          enddo
        enddo

        deallocate(tx)
        deallocate(ty)
        deallocate(work)
        deallocate(bcoef)

      endif

c End program

      return
      end subroutine

*deck restrict
c######################################################################
      subroutine restrict(xc,ntotc,nxc,nyc,xf,ntotf,nxf,nyf,order
     .                   ,igf,volf)
c----------------------------------------------------------------------
c     This is a restriction routine for MG. 
c
c     If order = 0, it employs simple agglomeration.
c     If order > 1, it employs spline interpolation of order "order".
c----------------------------------------------------------------------

      use mg_setup

      implicit none        !For safe fortran

c Call variables

      integer*4    ntotc,nxc,nyc,ntotf,nxf,nyf,order,igf

      real*8       xc(ntotc),xf(ntotf),volf

c Local variables
 
      real*8       uf(0:nxf+1,0:nyf+1)

      real*8       xxc,yyc,xx(nxf+2),yy(nyf+2)

      integer*4    ic,if,jc,jf,iif,iic,i,igc

c Interpolation

      integer*4    kx,ky,nx,ny,dim,flg
      real*8, dimension(:)  ,allocatable:: tx,ty,work
      real*8, dimension(:,:),allocatable:: bcoef

      real*8       db2val
      external     db2val

c Begin program

c Agglomeration

      if (order.eq.0) then

        do jc = 1,nyc
          do ic = 1,nxc
            jf = 2*jc
            if = 2*ic
            iic = ic + nxc*(jc-1)
            iif = if + nxf*(jf-1)

            xc(iic) = (xf(iif)     + xf(iif-1)
     .               + xf(iif-nxf) + xf(iif-nxf-1) )/4d0*volf
          enddo
        enddo

      else

c Interpolation

c     Map vectors into arrays

        xx(1:nxf+2) = xl(0:nxf+1,igf)
        yy(1:nyf+2) = yl(0:nyf+1,igf)

        do if = 1,nxf
          do jf = 1,nyf
            iif = if + nxf*(jf-1)
            uf(if,jf) = xf(iif)
          enddo
        enddo

c     Impose boundary conditions on arrays (hardwired --> fix)

c       Top and bottom (residuals are zero at Dirichlet boundaries)

        do if=0,nxf+1
          uf(if,nyf+1) = 0d0
          uf(if,0)     = 0d0
        enddo

c       Periodic BC's

        do jf=0,nyf+1
          uf(0    ,jf) = uf(nxf-1,jf)
          uf(nxf+1,jf) = uf(2    ,jf)
        enddo

c     Calculate interpolation

        flg = 0
        kx = order + 1
        ky = order + 1
        nx = nxf+2
        ny = nyf+2
        dim = nx*ny + max(2*kx*(nx+1),2*ky*(ny+1))

        allocate(tx(nx+kx))
        allocate(ty(nx+ky))
        allocate(work(dim))
        allocate(bcoef(nx,ny))

        call db2ink(xx,nx,yy,ny,uf,nx,kx,ky,tx,ty,bcoef,work,flg)

        do jc = 1,nyc
          do ic = 1,nxc

            iic = ic + nxc*(jc-1)

            igc = igf-1         !Consider next coarser grid

            xxc = xl(ic,igc)
            yyc = yl(jc,igc)
            
            xc(iic) = 
     .           volf*db2val(xxc,yyc,0,0,tx,ty,nx,ny,kx,ky,bcoef,work)

          enddo
        enddo

        deallocate(tx)
        deallocate(ty)
        deallocate(work)
        deallocate(bcoef)

      endif

c End program

      return
      end subroutine

c cjacobi
c####################################################################
      subroutine cjacobi (neq,ntot,rr,zz,matvec,diag,igrid,iter
     .                   ,omega,tol,guess,out)
c-------------------------------------------------------------------
c     Does matrix-free jacobi iteration on system A*zz = rr. 
c     In call, zz contains initial guess.
c-------------------------------------------------------------------

      use mg_setup

      implicit none        !For safe fortran

c Call variables

      integer*4   ntot,igrid,iter,guess,out,neq
      real*8      rr(ntot),zz(ntot),diag(neq,ntot)
      real*8      omega,tol

      external    matvec

c Local variables

      integer*4   i,j,itr,nn,ieq,ii,iii
      real*8      mag0,mag1,mag,yy(ntot)
      real*8      dummy(neq),rhs(neq),delta

c Begin program

      if (out.ge.2) write (*,*)

      nn = ntot

      if (guess.eq.0) then
        do ii=1,nn
          zz(ii) = 0d0
        enddo
      endif

c Jacobi iteration

      do itr=1,iter

        call matvec(0,ntot,zz,yy,igrid)

        mag1 = mag
        mag  = 0d0

        select case (neq)
        case (1)

          do ii = 1,nn
            rhs(neq) = rr(ii) - yy(ii)
            dummy(neq) = rhs(neq)/diag(neq,ii)
            zz(ii) = zz(ii) + omega*dummy(neq)
            mag = mag + rhs(neq)**2
          enddo

        case (2)

          do ii = 1,nn/neq

            do ieq = 1,neq
              iii = neq*(ii-1) + ieq 
              rhs(ieq) = rr(iii) - yy(iii)
            enddo

            iii = neq*(ii-1)
            delta = diag(1,iii+1)*diag(2,iii+2)
     .             -diag(2,iii+1)*diag(1,iii+2)
            dummy(1) = (rhs(1)*diag(2,iii+2)-rhs(2)*diag(2,iii+1))/delta
            dummy(2) = (rhs(2)*diag(1,iii+1)-rhs(1)*diag(1,iii+2))/delta

            do ieq=1,neq
              zz(iii+ieq) = zz(iii+ieq) + omega*dummy(ieq)
            enddo

            mag = mag + sum(rhs**2)
cc            mag = mag + rhs(2)**2
          enddo

        case default

          write (*,*) 'Jacobi smoothing for neq > 2 not implemented yet'
          stop

        end select

c     Check convergence

        mag = sqrt(mag)

        if (itr.eq.1) then
          mag0 = mag
          if (out.ge.2) write (*,10) itr,mag,mag/mag0
        else
          if (out.ge.2) write (*,20) itr,mag,mag/mag0,mag/mag1
cc          if (mag/mag0.lt.tol.or.mag.lt.1d-16) exit
          if (mag/mag0.lt.tol) exit
        endif

      enddo

      if (out.ge.2) write (*,*)
      if (out.eq.1) write (*,10) itr,mag,mag/mag0

cc      iter = itr 
cc      tol  = mag

c End program

      return

 10   format (' JB Iteration:',i4,'; Residual:',1p1e10.2,
     .        '; Ratio:',1p1e10.2)
 20   format (' JB Iteration:',i4,'; Residual:',1p1e10.2,
     .        '; Ratio:',1p1e10.2,'; Damping:',1p1e10.2)

      end subroutine

c csgs
c####################################################################
      subroutine csgs(neq,ntot,rr,zz,matvec,diag,igrid,iter
     .               ,omega,tol,guess,out)
      implicit none       !For safe fortran
c-------------------------------------------------------------------
c     Does matrix-free SGS iteration on system, A zz = rr. In call, 
c     zz contains initial guess.
c-------------------------------------------------------------------

c Call variables

      integer*4   ntot,igrid,iter,guess,out,neq
      real*8      rr(ntot),zz(ntot),diag(neq,ntot)
      real*8      omega,tol

      external    matvec

c Local variables

      integer*4   i,j,itr,nn,ieq,ii,iii
      real*8      mag0,mag1,mag,yy(ntot)
      real*8      dummy(neq),rhs(neq),delta

c Begin program

      if (out.ge.2) write (*,*)

      nn = ntot

      if (guess.eq.0) then
        do ii=1,nn
          zz(ii) = 0d0
        enddo
      endif

c SGS iteration

      do itr=1,iter

        mag1 = mag
        mag = 0d0

        select case (neq)
        case (1)

c       Forward pass

          do ii = 1,nn
            call matvec(ii,ntot,zz,yy,igrid)
            rhs(neq) = rr(ii) - yy(ii)
            dummy(neq) = rhs(neq)/diag(neq,ii)
            zz(ii) = zz(ii) + omega*dummy(neq)
          enddo

c       Backward pass

          do ii = nn,1,-1
            call matvec(ii,ntot,zz,yy,igrid)
            rhs(neq) = rr(ii) - yy(ii)
            dummy(neq) = rhs(neq)/diag(neq,ii)
            zz(ii) = zz(ii) + omega*dummy(neq)
            mag = mag + rhs(neq)**2
          enddo

        case (2)

c       Forward pass

          do ii = 1,nn/neq

            do ieq = 1,neq
              iii = neq*(ii-1) + ieq 
              call matvec(iii,ntot,zz,yy,igrid)
              rhs(ieq) = rr(iii) - yy(iii)
            enddo

            iii = neq*(ii-1)
            delta = diag(1,iii+1)*diag(2,iii+2)
     .             -diag(2,iii+1)*diag(1,iii+2)
            dummy(1) = (rhs(1)*diag(2,iii+2)-rhs(2)*diag(2,iii+1))/delta
            dummy(2) = (rhs(2)*diag(1,iii+1)-rhs(1)*diag(1,iii+2))/delta

            do ieq=1,neq
              zz(iii+ieq) = zz(iii+ieq) + omega*dummy(ieq)
            enddo

          enddo

c       Backward pass

          do ii = nn/neq,1,-1

            do ieq = 1,neq
              iii = neq*(ii-1) + ieq 
              call matvec(iii,ntot,zz,yy,igrid)
              rhs(ieq) = rr(iii) - yy(iii)
            enddo

            iii = neq*(ii-1)
            delta = diag(1,iii+1)*diag(2,iii+2)
     .             -diag(2,iii+1)*diag(1,iii+2)
            dummy(1) = (rhs(1)*diag(2,iii+2)-rhs(2)*diag(2,iii+1))/delta
            dummy(2) = (rhs(2)*diag(1,iii+1)-rhs(1)*diag(1,iii+2))/delta

            do ieq=1,neq
              zz(iii+ieq) = zz(iii+ieq) + omega*dummy(ieq)
            enddo

            mag = mag + sum(rhs**2)
cc            mag = mag + rhs(2)**2
          enddo

        case default

          write (*,*) 'SGS smoothing for neq > 2 not implemented yet'
          stop

        end select

c     Check convergence

        mag = sqrt(mag)

        if (itr.eq.1) then
          mag0 = mag
          if (out.ge.2) write (*,10) itr,mag,mag/mag0
        else
          if (out.ge.2) write (*,20) itr,mag,mag/mag0,mag/mag1
          if (mag/mag0.lt.tol) exit
        endif

      enddo

      if (out.ge.2) write (*,*)
      if (out.eq.1) write (*,10) itr,mag,mag/mag0

cc      iter = itr
cc      tol  = mag

c End program

      return

 10   format (' SGS Iteration:',i4,'; Residual:',1p1e10.2,
     .        '; Ratio:',1p1e10.2)
 20   format (' SGS Iteration:',i4,'; Residual:',1p1e10.2,
     .        '; Ratio:',1p1e10.2,'; Damping:',1p1e10.2)

      end subroutine

c find_mf_diag_neq
c####################################################################
      subroutine find_mf_diag_neq(neq,ntot,matvec,diag,ngrid,igridmin)
c--------------------------------------------------------------------
c     Finds diagonal elements matrix-free using subroutine matvec.
c--------------------------------------------------------------------

      use mg_setup

      implicit none      !For safe fortran

c Call variables

      integer*4      neq,ntot,ngrid,igridmin
      real*8         diag(neq,2*ntot)

      external       matvec

c Local variables

      real*8         x1(ntot),dummy(ntot)
      integer*4      i,j,ii,nn,jj,i1,j1,i2,j2
      integer*4      nx,ny,igrid,iii,ieq

c Begin program

      do igrid = ngrid,igridmin,-1

c     Form diagonal terms for smoother

        nx  = nxvp(igrid)
        ny  = nyvp(igrid)

        nn = neq*nx*ny

        do ii = 1,nn
          x1    (ii) = 0d0
          dummy (ii) = 0d0
        enddo

c Finds block diagonals for neq equations.

        do ii = 1,nn/neq

          jj  = (ii-1)*neq + istart(igrid) - 1

          do ieq = 1,neq

c       Find column vector ii

          if (mod(ii,nx).eq.1) then
            iii = neq*(ii-1) + ieq
            x1(iii) = 1d0
            iii = neq*(ii+nx-2) + ieq
            x1(iii) = 1d0
          elseif (mod(ii,nx).eq.0) then
            iii = neq*(ii-1) + ieq
            x1(iii) = 1d0
            iii = neq*(ii-nx) + ieq
            x1(iii) = 1d0
          else
            iii = neq*(ii-1) + ieq
            x1(iii) = 1d0
          endif

          call matvec(ii,nn,x1,dummy,igrid)

          if (mod(ii,nx).eq.1) then
            iii = neq*(ii-1) + ieq
            x1(iii) = 0d0
            iii = neq*(ii+nx-2) + ieq
            x1(iii) = 0d0
          elseif (mod(ii,nx).eq.0) then
            iii = neq*(ii-1) + ieq
            x1(iii) = 0d0
            iii = neq*(ii-nx) + ieq
            x1(iii) = 0d0
          else
            iii = neq*(ii-1) + ieq
            x1(iii) = 0d0
          endif

          diag(ieq,jj+1:jj+neq) = dummy((ii-1)*neq+1:ii*neq)

          enddo

        enddo

      enddo

c End program

      return
      end subroutine
