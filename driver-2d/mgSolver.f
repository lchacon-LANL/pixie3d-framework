*deck mgprep
c #####################################################################
      subroutine mgprep(pos,idiagp,x,y,igridn,guess)

c ---------------------------------------------------------------------
c   This subroutine solves, A x = y for the vector x
c   using a multigrid correction scheme. The coarse
c   grid operators are also stored in (a).
c
c   The integer array idiagp contains the information
c   that relates the actual column to the main diagonal
c   since only the non-zero diagonals of A are stored, i.e.,
c
c     i = row number
c     j = column number = i + idiagp(k)
c     k = non-zero diagonal index numbered from left to right
c ---------------------------------------------------------------------

      use parameters

      use mg_setup

      use precond_setup

      implicit none

c Call variables

      integer*4  igridn,idiagp(ndiagdp,igridn),guess
      real*8     pos(ndiagdp,ntotd2p),x(ntotdp),y(ntotdp)

c Local variables

      real*8     xx(ntotd2p),yy(ntotd2p),wrk(ntotd2p)
      real*8     rave,rold
      integer*4  izero,i,j,ii,ivcyc,igrid,igp

c External

      external matmul

c Begin program

c Initialize local solution vector (xx) and local r.h.s. (yy)

      do i = 1,ntotd2p
        xx(i) = 0d0
        yy(i) = 0d0
        wrk(i)= 0d0
      enddo

      if (guess.eq.0) then

        do i = 1,ntotvp(igridn)
          ii = istartp(igridn) + i - 1
          yy(ii) = y(i)
        enddo

      else

        do i = 1,ntotvp(igridn)
          ii = istartp(igridn) + i - 1
          yy(ii) = y(i)
          xx(ii) = x(i)
        enddo

      endif

c  Start V-cycle

      do ivcyc = 1,maxvcyc

c     Cycle down to coasest grid, igrid = 1 

        do igrid = igridn,1,-1

          if(igrid .ne. 1) then

c         Relax error on grid number igrid

            izero = 1
            if(ivcyc .gt. 1 .and. igrid .eq. igridn) izero = 0

            call sgs_mg(ntotvp(igrid),nxvp(igrid),nyvp(igrid),
     &                  pos(1,istartp(igrid)),
     &                  idiagp(1,igrid),ndiagdp,xx(istartp(igrid)),
     &                                yy(istartp(igrid)),izero,igrid)

c         Evaluate updated R.H.S. (ie yy)

c           wrk = A xx

            call matmul(ntotvp(igrid),pos(1,istartp(igrid)),
     &             xx(istartp(igrid)),wrk,ndiagdp,idiagp(1,igrid),igrid)

c           yy = yy - A xx = yy - wrk

            call vecad2(ntotvp(igrid),-1.d0,wrk,1.d0,yy(istartp(igrid)))

c         Restrict R.H.S ( i.e. yy ) down a grid

c           yy_c = R * yy_f

            call mgrestp(yy(istartp(igrid-1)),ntotvp(igrid-1),
     &              nxvp(igrid-1),nyvp(igrid-1),
     &           wrk,ntotvp(igrid),nxvp(igrid),nyvp(igrid))

          else

c         Solve for error on igrid = 1 (i.e. xx_1)

            call sgs_mg(ntotvp(igrid),nxvp(igrid),nyvp(igrid),
     &                  pos(1,istartp(igrid)),
     &                  idiagp(1,igrid),ndiagdp,xx(istartp(igrid)),
     &                                    yy(istartp(igrid)),1,igrid)

c         Prolongate error (i.e. xx) on grid 1 to grid 2

c           wrk = P * xx_1 

            if (igridn .gt. 1) then

              call mgprolp(wrk,ntotvp(igrid+1),
     &              nxvp(igrid+1),nyvp(igrid+1),
     &              xx(istartp(igrid)),ntotvp(igrid),
     &              nxvp(igrid),nyvp(igrid))

            endif

c         Update existing error on grid 2 (i.e. xx_2)

c           xx_2 = xx_2 + wrk 

            igp = igrid+1
            call vecad2(ntotvp(igp),1d0,xx(istartp(igp)),1d0,wrk)

          endif

        enddo

c     Cycle back up to grid number igridn updating errors (xx)
c     with fixed R.H.S. (yy)

        do igrid = 2,igridn

c       Relax updated error on igrid (i.e. xx_igrid)

          call sgs_mg(ntotvp(igrid),nxvp(igrid),nyvp(igrid),
     &                  pos(1,istartp(igrid)),
     &                  idiagp(1,igrid),ndiagdp,xx(istartp(igrid)),
     &                                    yy(istartp(igrid)),0,igrid)

          if( igrid .ne. igridn )then

c         Prolongate updated error  (i.e. xx) to igrid + 1 
c         wrk = P * xx_igrid

            call mgprolp(wrk,ntotvp(igrid+1),
     &              nxvp(igrid+1),nyvp(igrid+1),
     &              xx(istartp(igrid)),ntotvp(igrid),
     &              nxvp(igrid),nyvp(igrid) )

c         Ppdate existing error on grid igrid+1 (i.e. xx_igrid+1)
c         xx_igrid+1 = xx_igrid+1 + wrk 

            igp = igrid + 1
            call vecad2(ntotvp(igp),1d0,xx(istartp(igp)),1d0,wrk)

          endif

        enddo

      enddo               !End V-cycle loop

c Transfer local solution to global solution vector

      do i = 1,ntotvp(igridn)
        ii = istartp(igridn) + i - 1
        x(i) = xx(ii)
      enddo

c End program

      return
      end

c mgprolp
c##################################################################
      subroutine mgprolp(xf,ntotf,nxf,nyf,xc,ntotc,nxc,nyc)
cc      implicit none            ! For safe Fortran
c------------------------------------------------------------------
c     This is a 1-D prolongation routine for multigrid
c     preconditioning, using bilinear interpolation
c------------------------------------------------------------------

      implicit none

c Call variables

      integer      ntotc,nxc,nyc,ntotf,nxf,nyf
      real*8       xc(ntotc),xf(ntotf)

c Local variables
 
      real*8       uc(0:nxc+1,0:nyc+1),vc(0:nxc+1,0:nyc+1),
     .             uf(nxf,nyf),vf(nxf,nyf),c(4,5),ff

      integer*4    ic,if,jc,jf,ii,i,ig

c Begin program

cc      ig = floor(log(1d0*nxd/nxc)/log(2d0)+0.5)

cc      ff = 2.**(ig-1)/(1.+2.**ig)

cc      write (*,*) nxd,nxc,ig,ff
cc      stop

c Map vectors into arrays

      do ic = 1,nxc
        do jc = 1,nyc
          ii = ic + nxc*(jc-1)
          uc(ic,jc) = xc(ii)
        enddo
      enddo

c Impose boundary conditions on arrays

c     Top and bottom (Newton updates are zero at boundary)

      do ic=0,nxc+1
        uc(ic,nyc+1) = 0d0
        uc(ic,0)     = 0d0
      enddo

c     Periodic BC's

      do jc=0,nyc+1
        uc(0    ,jc) = uc(nxc,jc)
        uc(nxc+1,jc) = uc(1  ,jc)
      enddo

c Calculate interpolation

      do jc = 1,nyc
        do ic = 1,nxc
          jf = 2*jc
          if = 2*ic

c         Binomial interpolation coefficients

          do i = 1,4
            c(i,1) = 9.0d0
            c(i,2) = 3.0d0
            c(i,3) = 3.0d0
            c(i,4) = 1.0d0
            c(i,5) = 16.0d0
          enddo

c         Interpolation coefficients for dirichlet and periodic BC's

          if(ic .ne. 1 .and. ic .ne. nxc .and.
     .       jc .ne. 1 .and. jc .ne. nyc)then

            goto 3000

          elseif(ic .eq. 1 .and.
     .           jc .ne. 1 .and. jc .ne. nyc)then

            c(2,1) = 3d0
            c(2,2) = 3d0
            c(2,3) = 1d0
            c(2,4) = 1d0
            c(2,5) = 8d0

            c(3,1) = 3d0
            c(3,2) = 3d0
            c(3,3) = 1d0
            c(3,4) = 1d0
            c(3,5) = 8d0

          elseif(ic .eq. nxc.and.
     .           jc .ne. 1  .and. jc .ne. nyc)then

            c(1,1) = 3d0
            c(1,2) = 3d0
            c(1,3) = 1d0
            c(1,4) = 1d0
            c(1,5) = 8d0

            c(4,1) = 3d0
            c(4,2) = 3d0
            c(4,3) = 1d0
            c(4,4) = 1d0
            c(4,5) = 8d0

          elseif(jc .eq. 1 .and.
     .           ic .ne. 1 .and. ic .ne. nxc)then

cc            c(4,1) = 9d0*(1.-ff)
cc            c(4,2) = 3d0*(1.-ff)
cc            c(4,3) = 9.*ff
cc            c(4,4) = 3.*ff
cc            c(4,5) = 12d0
     
cc            c(3,1) = 9d0*(1.-ff)
cc            c(3,2) = 3d0*(1.-ff)
cc            c(3,3) = 9.*ff
cc            c(3,4) = 3.*ff
cc            c(3,5) = 12d0

            c(4,1) = 6d0
            c(4,2) = 2d0
            c(4,3) = 3d0
            c(4,4) = 1d0
            c(4,5) = 12d0
     
            c(3,1) = 6d0
            c(3,2) = 2d0
            c(3,3) = 3d0
            c(3,4) = 1d0
            c(3,5) = 12d0

          elseif(jc .eq. nyc.and.
     .           ic .ne. 1  .and. ic .ne. nxc)then

            c(1,1) = 6d0
            c(1,2) = 2d0
            c(1,3) = 3d0
            c(1,4) = 1d0
            c(1,5) = 12d0
     
            c(2,1) = 6d0
            c(2,2) = 2d0
            c(2,3) = 3d0
            c(2,4) = 1d0
            c(2,5) = 12d0

          elseif(ic .eq. 1 .and. jc .eq. 1) then

            c(2,1) = 2d0
            c(2,2) = 2d0
            c(2,3) = 1d0
            c(2,4) = 1d0
            c(2,5) = 6d0

            c(3,1) = 2d0
            c(3,2) = 2d0
            c(3,3) = 1d0
            c(3,4) = 1d0
            c(3,5) = 6d0

            c(4,1) = 2d0
            c(4,2) = 2d0
            c(4,3) = 1d0
            c(4,4) = 1d0
            c(4,5) = 6d0

          elseif(ic .eq. 1  .and. jc .eq. nyc) then

            c(1,1) = 2d0
            c(1,2) = 2d0
            c(1,3) = 1d0
            c(1,4) = 1d0
            c(1,5) = 6d0

            c(2,1) = 2d0
            c(2,2) = 2d0
            c(2,3) = 1d0
            c(2,4) = 1d0
            c(2,5) = 6d0

            c(3,1) = 2d0
            c(3,2) = 2d0
            c(3,3) = 1d0
            c(3,4) = 1d0
            c(3,5) = 6d0

          elseif(ic .eq. nxc  .and. jc .eq. 1) then

            c(1,1) = 2d0
            c(1,2) = 2d0
            c(1,3) = 1d0
            c(1,4) = 1d0
            c(1,5) = 6d0

            c(3,1) = 2d0
            c(3,2) = 2d0
            c(3,3) = 1d0
            c(3,4) = 1d0
            c(3,5) = 6d0

            c(4,1) = 2d0
            c(4,2) = 2d0
            c(4,3) = 1d0
            c(4,4) = 1d0
            c(4,5) = 6d0

          elseif(ic .eq. nxc  .and. jc .eq. nyc) then

            c(1,1) = 2d0
            c(1,2) = 2d0
            c(1,3) = 1d0
            c(1,4) = 1d0
            c(1,5) = 6d0

            c(2,1) = 2d0
            c(2,2) = 2d0
            c(2,3) = 1d0
            c(2,4) = 1d0
            c(2,5) = 6d0

            c(4,1) = 2d0
            c(4,2) = 2d0
            c(4,3) = 1d0
            c(4,4) = 1d0
            c(4,5) = 6d0

          endif

 3000     continue

          uf(if,jf)     = (c(1,1)*uc(ic,jc)+c(1,2)*uc(ic+1,jc)
     .          +c(1,3)*uc(ic,jc+1)+c(1,4)*uc(ic+1,jc+1))/c(1,5)
          uf(if-1,jf)   = (c(2,1)*uc(ic,jc)+c(2,2)*uc(ic-1,jc)
     .          +c(2,3)*uc(ic,jc+1)+c(2,4)*uc(ic-1,jc+1))/c(2,5)
          uf(if-1,jf-1) = (c(3,1)*uc(ic,jc)+c(3,2)*uc(ic-1,jc)
     .          +c(3,3)*uc(ic,jc-1)+c(3,4)*uc(ic-1,jc-1))/c(3,5)
          uf(if,jf-1)   = (c(4,1)*uc(ic,jc)+c(4,2)*uc(ic+1,jc)
     .          +c(4,3)*uc(ic,jc-1)+c(4,4)*uc(ic+1,jc-1))/c(4,5)

         enddo
       enddo

c Map arrays in vector

       do if = 1,nxf
         do jf = 1,nyf
           ii = if + nxf*(jf-1)
           xf(ii) = uf(if,jf) 
         enddo
       enddo

c End program

      return
      end

c mgrestp
c######################################################################
      subroutine mgrestp(xc,ntotc,nxc,nyc,xf,ntotf,nxf,nyf)
      implicit none        !For safe fortran
c----------------------------------------------------------------------
c     This is a piecewise constant restriction routine for MG
c----------------------------------------------------------------------

c Call variables

      integer*4    ntotc,nxc,nyc,ntotf,nxf,nyf

      real*8       xc(ntotc), xf(ntotf)

c Local variables

      integer*4    ic,if,jc,jf,iic,iif

c Begin program

c Restrict (piece-wise constant)

      do jc = 1,nyc
        do ic = 1,nxc
          jf = 2*jc
          if = 2*ic
          iic = ic + nxc*(jc-1)
          iif = if + nxf*(jf-1)

          xc(iic) = xf(iif)     + xf(iif-1)
     .            + xf(iif-nxf) + xf(iif-nxf-1)

        enddo
      enddo

c End program

      return
      end

c mgrestvar
c######################################################################
      subroutine mgrestvar(xc,ntotc,nxc,nyc,xf,ntotf,nxf,nyf)
      implicit none        !For safe fortran
c----------------------------------------------------------------------
c     This is a piecewise constant restriction routine for MG
c----------------------------------------------------------------------

c Call variables

      integer*4    ntotc,nxc,nyc,ntotf,nxf,nyf

      real*8       xc(ntotc), xf(ntotf)

c Local variables

      integer*4    ic,if,jc,jf,iic,iif

c Begin program

c Restrict (piece-wise constant)

      do jc = 1,nyc
        do ic = 1,nxc
          jf = 2*jc
          if = 2*ic
          iic = ic + nxc*(jc-1)
          iif = if + nxf*(jf-1)

          xc(iic) = ( xf(iif)     + xf(iif-1)
     .              + xf(iif-nxf) + xf(iif-nxf-1))/4d0

        enddo
      enddo

c End program

      return
      end

*deck matmul
c ####################################################################
      subroutine matmul(nn,a,x,y,ndiag,idiag,igrid)
c --------------------------------------------------------------------
c   This subroutine computes
c
c     y = A*x
c
c   where A is the Jacobian matrix of a 2-D
c   natural convection problem.
c
c   The integer array idiag contains the information
c   that relates the actual column to the main diagonal
c   since only the non-zero diagonals of A are stored, i.e.,
c
c     i = row number
c     j = column number = i + idiag(k)
c     k = non-zero diagonal index numbered from left to right
c --------------------------------------------------------------------

      use mg_setup

      implicit none

c Call variables

      integer*4   nn,ndiag,idiag(ndiag),igrid
      real*8      a(ndiag,nn),x(nn),y(nn)

c Local variables

      integer*4   i,ig,jg,iim,ip,nx,ny

c Begin program

      nx = nxvp(igrid)
      ny = nyvp(igrid)

      do jg = 1,ny
        do ig = 1,nx

          i = ig + nx*(jg-1)

          iim = i-1
          if (ig.eq.1) iim = nx-1 + nx*(jg-1) 

          ip = i+1
          if (ig.eq.nx) ip = 2 + nx*(jg-1)

          if (i .le. nx) then
            y(i) = a(2,i) * x(iim )
     &           + a(3,i) * x(i   )
     &           + a(4,i) * x(ip  )
     &           + a(5,i) * x(i+nx)
          elseif (i .gt. nn-nx) then
            y(i) = a(1,i) * x(i-nx)
     &           + a(2,i) * x(iim )
     &           + a(3,i) * x(i   )
     &           + a(4,i) * x(ip  )
          else
            y(i) = a(1,i) * x(i-nx)
     &           + a(2,i) * x(iim )
     &           + a(3,i) * x(i   )
     &           + a(4,i) * x(ip  )
     &           + a(5,i) * x(i+nx)
          endif

        enddo
      enddo

      return
      end

*deck sgs_mg
c ######################################################################
      subroutine sgs_mg(ntot,nx,ny,a,idiag,ndiag,x,y,izero,igrid)

c ----------------------------------------------------------------------
c   This subroutine solves, Ax = y for the vector x using a point 
c   SGS iteration.
c
c   The integer array idiag contains the information
c   that relates the actual column to the main diagonal
c   since only the non-zero diagonals of A are stored, i.e.,
c
c     i = row number
c     j = column number = i + idiag(k)
c     k = non-zero diagonal index numbered from left to right
c----------------------------------------------------------------------

      use precond_setup

      implicit none

c Input variables

      integer*4    ntot,nx,ny,ndiag,izero,igrid,idiag(ndiag)
      real*8       a(ndiag,ntot),x(ntot),y(ntot)

c Local variables

      integer*4    nn,i,j,k,mdiag,ii,ip,iim,ig,jg

c Begin program

      nn = ntot
      mdiag = 3

      do j = 1,nsweep

c     Forward pass

       do jg = 1,ny
       do ig = 1,nx

         i = ig + nx*(jg-1)

         iim = i-1
         if (ig.eq.1) iim = nx-1 + nx*(jg-1) 

         ip = i+1
         if (ig.eq.nx) ip = 2 + nx*(jg-1)

         if (i .le. nx) then
           x(i) = ( y(i) 
     &                   - (a(2,i) * x(iim ))
     &                   - (a(4,i) * x(ip  ))
     &                   - (a(5,i) * x(i+nx))
     &                                           ) / a(3,i)
         elseif (i .gt. nn-nx) then
           x(i) = ( y(i) - (a(1,i) * x(i-nx))
     &                   - (a(2,i) * x(iim ))
     &                   - (a(4,i) * x(ip  ))
     &                                           ) / a(3,i)
         else
           x(i) = ( y(i) - (a(1,i) * x(i-nx))
     &                   - (a(2,i) * x(iim ))
     &                   - (a(4,i) * x(ip  ))
     &                   - (a(5,i) * x(i+nx))
     &                                           ) / a(3,i)
         endif

       enddo
       enddo

c     Backward pass

       do jg = ny,1,-1
       do ig = nx,1,-1

         i = ig + nx*(jg-1)

         iim = i-1
         if (ig.eq.1) iim = nx-1 + nx*(jg-1) 

         ip = i+1
         if (ig.eq.nx) ip = 2 + nx*(jg-1)

         if (i .le. nx) then
           x(i) = ( y(i) 
     &                   - (a(2,i) * x(iim ))
     &                   - (a(4,i) * x(ip  ))
     &                   - (a(5,i) * x(i+nx))
     &                                           ) / a(3,i)
         elseif (i .gt. nn-nx) then
           x(i) = ( y(i) - (a(1,i) * x(i-nx))
     &                   - (a(2,i) * x(iim ))
     &                   - (a(4,i) * x(ip  ))
     &                                           ) / a(3,i)
         else
           x(i) = ( y(i) - (a(1,i) * x(i-nx))
     &                   - (a(2,i) * x(iim ))
     &                   - (a(4,i) * x(ip  ))
     &                   - (a(5,i) * x(i+nx))
     &                                           ) / a(3,i)
         endif
       enddo
       enddo

      enddo

      return
      end

c jacobi
c ######################################################################
      subroutine jacobi_mg(ntot,nx,ny,a,am,idiag,ndiag,x,y,izero,igrid)

c ----------------------------------------------------------------------
c   This subroutine solves, Ax = y for the vector x using a point 
c   Jacobi iteration.
c
c   The integer array idiag contains the information
c   that relates the actual column to the main diagonal
c   since only the non-zero diagonals of A are stored, i.e.,
c
c     i = row number
c     j = column number = i + idiag(k)
c     k = non-zero diagonal index numbered from left to right
c----------------------------------------------------------------------

      use precond_setup

      implicit none

c Input variables

      integer*4    ntot,nx,ny,ndiag,izero,igrid,idiag(ndiag)
      real*8       a(ndiag,ntot),am(ndiag,ntot),x(ntot),y(ntot)

c Local variables

      integer*4    nn,j,ii,ip,iim,ig,jg
      real*8       omega,xold(ntot)

c Begin program

      nn = ntot

      omega = .7

c Jacobi sweep

      do j = 1,nsweep

        do jg = 1,ny
          do ig = 1,nx
            ii = ig + nx*(jg-1)
            xold(ii) = x(ii)
          enddo
        enddo

        do jg = 1,ny
          do ig = 1,nx

            ii = ig + nx*(jg-1)

            iim = ii-1
            if (ig.eq.1) iim = nx-1 + nx*(jg-1) 

            ip = ii+1
            if (ig.eq.nx) ip = 2 + nx*(jg-1)

            if (ii.le.nx) then
              x(ii) = (1.-omega)*xold(ii)
     &             + omega*( y(ii) - (a(2,ii) * xold(iim  ))
     &                             - (a(4,ii) * xold(ip   ))
     &                             - (a(5,ii) * xold(ii+nx))
     &                                                    ) / a(3,ii)
            elseif (ii.gt.nn-nx) then
              x(ii) = (1.-omega)*xold(ii)
     &             + omega*( y(ii) - (a(1,ii) * xold(ii-nx))
     &                             - (a(2,ii) * xold(iim  ))
     &                             - (a(4,ii) * xold(ip   ))
     &                                                    ) / a(3,ii)
            else
              x(ii) = (1.-omega)*xold(ii)
     &             + omega*( y(ii) - (a(1,ii) * xold(ii-nx))
     &                             - (a(2,ii) * xold(iim  ))
     &                             - (a(4,ii) * xold(ip   ))
     &                             - (a(5,ii) * xold(ii+nx))
     &                                                    ) / a(3,ii)
            endif

          enddo
        enddo

      enddo

c End program

      return
      end
