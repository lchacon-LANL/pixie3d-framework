c module imposeBCinterface
c #####################################################################
      module imposeBCinterface

        use BCS_variables

        use singularBCinterface

        use dirichletBCinterface

        use neumannBCinterface

        INTERFACE setBC
          module procedure imposeBConScalar,imposeBConVector
        end INTERFACE

        integer(4),private :: bctype,dim,loc,ibc,i,j,k

        !Define order of interpolation for EDGE and CORNER treatment
        !(linear to preserve nearest-neighbors stencil; this is 
        ! important for colored diagonal formation)
        integer(4),private :: order=1

      contains

c     imposeBConScalar
c     #################################################################
      subroutine imposeBConScalar(ieq,nx,ny,nz,array,array0,bcond
     .                           ,iigx,iigy,iigz,i1,i2,j1,j2,k1,k2
     .                           ,iorder)
c     -----------------------------------------------------------------
c     Imposes BC on scalar quantities.
c     -----------------------------------------------------------------

      implicit none

c     Call variables

      integer(4) :: ieq,bcond(6),nx,ny,nz,iigx,iigy,iigz
      real(8)    :: array (0:nx+1,0:ny+1,0:nz+1)
     .             ,array0(0:nx+1,0:ny+1,0:nz+1)
      integer(4),optional,intent(IN) :: i1,i2,j1,j2,k1,k2,iorder

c     Local variables

c     Begin program

c     Initialize interpolation order

      if (PRESENT(iorder)) order = iorder

c     Initialize BC grid quantities

      igxbc = iigx
      igybc = iigy
      igzbc = iigz

      nxbc = nx
      nybc = ny
      nzbc = nz

c     Initialize global limits

      if (     PRESENT(i1).and.PRESENT(i2)
     .    .and.PRESENT(j1).and.PRESENT(j2)
     .    .and.PRESENT(k1).and.PRESENT(k2)) then
        iimin = i1
        iimax = i2
        jjmin = j1
        jjmax = j2
        kkmin = k1
        kkmax = k2
      else
        iimin = grid_params%ilo(iigx)
        iimax = grid_params%ihi(iigx)
        jjmin = grid_params%jlo(iigy)
        jjmax = grid_params%jhi(iigy)
        kkmin = grid_params%klo(iigz)
        kkmax = grid_params%khi(iigz)
      endif

c     Impose BCs

      do bctype=1,BCLIM             !Enforces a particular order in the BCs (see grid_mod.f)
        do dim=1,3
          do loc=0,1
            ibc = (1+loc)+2*(dim-1)
            if (bcond(ibc) == bctype) then
              call FillGhostNodes(ieq,1,1,dim,loc,bctype,array,array0)
            endif
          enddo
        enddo
      enddo

c     Singular point boundary condition

      if (bcond(1) == SP) call singularBC(array,2)

c     Fill edges

      do dim=1,3
        do loc=0,1
          call fillEdges(array,array0,dim,loc,bcond)
        enddo
      enddo

c     Fill corners

      do dim=1,3
        do loc=0,1
          call fillCorners(array,array0,dim,loc,bcond)
        enddo
      enddo

ccc     Synchronize periodic boundaries
cc
cc      bctype=PER
cc
cc      do dim=1,3
cc        do loc=0,1
cc          ibc = (1+loc)+2*(dim-1)
cc          if (bcond(ibc) == bctype) then
cc            call FillGhostNodes(ieq,1,1,dim,loc,bctype,array,array0)
cc          endif
cc        enddo
cc      enddo

c     End program

      end subroutine imposeBConScalar

c     imposeBConVector
c     #################################################################
      subroutine imposeBConVector(fcomp,nx,ny,nz,v_cnv,v_cov,var0,bcond
     .                           ,iigx,iigy,iigz,i1,i2,j1,j2,k1,k2
     .                           ,iorder)
c     -----------------------------------------------------------------
c     Imposes BC on vector field
c     -----------------------------------------------------------------

      implicit none

c     Call variables

      integer(4) :: bcond(6,3),fcomp,nx,ny,nz,iigx,iigy,iigz
      real(8)    :: v_cnv(0:nx+1,0:ny+1,0:nz+1,3)
     .             ,var0 (0:nx+1,0:ny+1,0:nz+1,3)
     .             ,v_cov(0:nx+1,0:ny+1,0:nz+1,3)
      integer(4),optional,intent(IN) :: i1,i2,j1,j2,k1,k2,iorder

c     Local variables

      integer(4) :: ivar,ieq,ibc,loc,dim,bctype
     .             ,ilmin,jlmin,klmin,ilmax,jlmax,klmax
      logical    :: cov_to_cnv

c     Begin program

c     Initialize interpolation order

      if (PRESENT(iorder)) order = iorder

c     Initialize BC grid quantities

      igxbc = iigx
      igybc = iigy
      igzbc = iigz

      nxbc = nx
      nybc = ny
      nzbc = nz

c     Initialize global limits

      if (     PRESENT(i1).and.PRESENT(i2)
     .    .and.PRESENT(j1).and.PRESENT(j2)
     .    .and.PRESENT(k1).and.PRESENT(k2)) then
        iimin = i1
        iimax = i2
        jjmin = j1
        jjmax = j2
        kkmin = k1
        kkmax = k2
      else
        iimin = grid_params%ilo(iigx)
        iimax = grid_params%ihi(iigx)
        jjmin = grid_params%jlo(iigy)
        jjmax = grid_params%jhi(iigy)
        kkmin = grid_params%klo(iigz)
        kkmax = grid_params%khi(iigz)
      endif

c     Find cov components (we assume cnv components are defined)

#if defined(petsc)
      call fromGlobalToLocalLimits(iimin,jjmin,kkmin,ilmin,jlmin,klmin)
      call fromGlobalToLocalLimits(iimax,jjmax,kkmax,ilmax,jlmax,klmax)

      do k=klmin-1,klmax+1
        do j=jlmin-1,jlmax+1
          do i=ilmin-1,ilmax+1
#else
      do k=kkmin-1,kkmax+1
        do j=jjmin-1,jjmax+1
          do i=iimin-1,iimax+1
#endif
            call transformFromCurvToCurv(i,j,k,igxbc,igybc,igzbc
     .            ,v_cov(i,j,k,1),v_cov(i,j,k,2),v_cov(i,j,k,3)
     .            ,v_cnv(i,j,k,1),v_cnv(i,j,k,2),v_cnv(i,j,k,3),.false.)
          enddo
        enddo
      enddo

c     Impose BCs

      cov_to_cnv = .false.

      do bctype=1,BCLIM            !Enforces a particular order in the BCs (see grid_mod.f)
        do dim=1,3
          do loc=0,1
            ibc = (1+loc)+2*(dim-1)

            do ivar = 1,3
              ieq = ivar + fcomp - 1

              if (abs(bcond(ibc,ivar)) == bctype) then

                if (bcond(ibc,ivar) < 0) then
                  call FillGhostNodes(ieq,ivar,3,dim,loc,bctype,v_cov
     .                               ,var0)
                  cov_to_cnv = .true.
                else
                  call FillGhostNodes(ieq,ivar,3,dim,loc,bctype,v_cnv
     .                               ,var0)
                endif

             endif

            enddo

          enddo
        enddo
      enddo

c     Synchronize covariant and contravariant components => contravariant

      if (cov_to_cnv) call synchronize_cnv(v_cnv,v_cov,bcond)

c     Impose vector singular point BCs

      if (bcond(1,1) == SP) call singularBC(v_cnv,.false.,2)

c     Fill edges

      do dim=1,3
        do loc=0,1
          do ivar=1,3
            call fillEdges(v_cnv(:,:,:,ivar),var0(:,:,:,ivar)
     .                    ,dim,loc,bcond)
          enddo
        enddo
      enddo

c     Fill corners

      do dim=1,3
        do loc=0,1
          do ivar=1,3
            call fillCorners(v_cnv(:,:,:,ivar),var0(:,:,:,ivar)
     .                      ,dim,loc,bcond)
          enddo
        enddo
      enddo

ccc     Synchronize periodic boundaries
cc
cc      bctype=PER
cc
cc      do dim=1,3
cc        do loc=0,1
cc          ibc = (1+loc)+2*(dim-1)
cc          do ivar = 1,3
cc            ieq = ivar + fcomp - 1
cc            if (bcond(ibc,ivar) == bctype) then
cc              call FillGhostNodes(ieq,ivar,3,dim,loc,bctype,v_cnv,var0)
cc            endif
cc          enddo
cc        enddo
cc      enddo

c     Find covariant components at ALL boundary ghost cells

      do dim = 1,3
        do loc = 0,1

          call findLoopLimits(dim,loc,iimin-1,iimax+1
     .                               ,jjmin-1,jjmax+1
     .                               ,kkmin-1,kkmax+1
     .                               ,imin,imax,jmin,jmax,kmin,kmax)

          do i=imin,imax
            do j=jmin,jmax
              do k=kmin,kmax
                call transformFromCurvToCurv(i,j,k,igxbc,igybc,igzbc
     .            ,v_cov(i,j,k,1),v_cov(i,j,k,2),v_cov(i,j,k,3)
     .            ,v_cnv(i,j,k,1),v_cnv(i,j,k,2),v_cnv(i,j,k,3),.false.)
              enddo
            enddo
          enddo

        enddo
      enddo

c     End program

      end subroutine imposeBConVector

c     synchronize_cnv
c     #################################################################
      subroutine synchronize_cnv(v_cnv,v_cov,bcond)
c     -----------------------------------------------------------------
c     Finds all contravariant components at boundaries with BC in terms
c     of covariant components.
c
c     On input, tangential covariant components and normal contravariant
c     components are known at ghost cells. On output, all contravariant
c     components are known at ghost cells.
c     -----------------------------------------------------------------

      implicit none

c     Call variables

      integer(4) :: bcond(6,3)
      real(8)    :: v_cnv(0:nxbc+1,0:nybc+1,0:nzbc+1,3)
     .             ,v_cov(0:nxbc+1,0:nybc+1,0:nzbc+1,3)

c     Local variables

      integer(4) :: i,j,k,dim,loc,ig,jg,kg,ivar
      real(8)    :: x1,x2,x3,gsuper(3,3),gsub(3,3)

c     Begin program

      do dim = 1,3
        do loc = 0,1
          ibc = (1+loc)+2*(dim-1)

          do ivar = 1,3

            if (ivar == dim) then  !Select tangential components
              cycle
            elseif (bcond(ibc,ivar) < 0) then

              call findLoopLimits(dim,loc,iimin,iimax
     .                                   ,jjmin,jjmax
     .                                   ,kkmin,kkmax
     .                                   ,imin,imax,jmin,jmax,kmin,kmax)

              if (do_nothing) exit

              select case (ibc)
              case (1)

                do i=imin,imax
                  do j=jmin,jmax
                    do k=kmin,kmax

                      gsuper = gmetric%grid(igxbc)%gsup(i-1,j,k,:,:)

                      v_cov(i-1,j,k,1) = -(gsuper(1,2)*v_cov(i-1,j,k,2)
     .                                    +gsuper(1,3)*v_cov(i-1,j,k,3)
     .                                    -v_cnv(i-1,j,k,1))/gsuper(1,1)

                      call transformFromCurvToCurv(i-1,j,k,igxbc,igybc,igzbc
     .               ,v_cov(i-1,j,k,1),v_cov(i-1,j,k,2),v_cov(i-1,j,k,3)
     .               ,v_cnv(i-1,j,k,1),v_cnv(i-1,j,k,2),v_cnv(i-1,j,k,3)
     .               ,.true.)

                    enddo
                  enddo
                enddo

              case (2)

                do i=imin,imax
                  do j=jmin,jmax
                    do k=kmin,kmax

                      gsuper = gmetric%grid(igxbc)%gsup(i+1,j,k,:,:)

                      v_cov(i+1,j,k,1) = -(gsuper(1,2)*v_cov(i+1,j,k,2)
     .                                    +gsuper(1,3)*v_cov(i+1,j,k,3)
     .                                    -v_cnv(i+1,j,k,1))/gsuper(1,1)

                      call transformFromCurvToCurv(i+1,j,k,igxbc,igybc,igzbc
     .               ,v_cov(i+1,j,k,1),v_cov(i+1,j,k,2),v_cov(i+1,j,k,3)
     .               ,v_cnv(i+1,j,k,1),v_cnv(i+1,j,k,2),v_cnv(i+1,j,k,3)
     .               ,.true.)

                    enddo
                  enddo
                enddo

              case (3)

                do i=imin,imax
                  do j=jmin,jmax
                    do k=kmin,kmax

                      gsuper = gmetric%grid(igxbc)%gsup(i,j-1,k,:,:)

                      v_cov(i,j-1,k,2) = -(gsuper(2,1)*v_cov(i,j-1,k,1)
     .                                    +gsuper(2,3)*v_cov(i,j-1,k,3)
     .                                    -v_cnv(i,j-1,k,2))/gsuper(2,2)

                      call transformFromCurvToCurv(i,j-1,k,igxbc,igybc,igzbc
     .               ,v_cov(i,j-1,k,1),v_cov(i,j-1,k,2),v_cov(i,j-1,k,3)
     .               ,v_cnv(i,j-1,k,1),v_cnv(i,j-1,k,2),v_cnv(i,j-1,k,3)
     .               ,.true.)

                    enddo
                  enddo
                enddo

              case (4)

                do i=imin,imax
                  do j=jmin,jmax
                    do k=kmin,kmax

                      gsuper = gmetric%grid(igxbc)%gsup(i,j+1,k,:,:)

                      v_cov(i,j+1,k,2) = -(gsuper(2,1)*v_cov(i,j+1,k,1)
     .                                    +gsuper(2,3)*v_cov(i,j+1,k,3)
     .                                    -v_cnv(i,j+1,k,2))/gsuper(2,2)

                      call transformFromCurvToCurv(i,j+1,k,igxbc,igybc,igzbc
     .               ,v_cov(i,j+1,k,1),v_cov(i,j+1,k,2),v_cov(i,j+1,k,3)
     .               ,v_cnv(i,j+1,k,1),v_cnv(i,j+1,k,2),v_cnv(i,j+1,k,3)
     .               ,.true.)

                    enddo
                  enddo
                enddo

              case (5)

                do i=imin,imax
                  do j=jmin,jmax
                    do k=kmin,kmax

                      gsuper = gmetric%grid(igxbc)%gsup(i,j,k-1,:,:)

                      v_cov(i,j,k-1,3) = -(gsuper(3,1)*v_cov(i,j,k-1,1)
     .                                    +gsuper(3,2)*v_cov(i,j,k-1,2)
     .                                    -v_cnv(i,j,k-1,3))/gsuper(3,3)

                      call transformFromCurvToCurv(i,j,k-1,igxbc,igybc,igzbc
     .               ,v_cov(i,j,k-1,1),v_cov(i,j,k-1,2),v_cov(i,j,k-1,3)
     .               ,v_cnv(i,j,k-1,1),v_cnv(i,j,k-1,2),v_cnv(i,j,k-1,3)
     .               ,.true.)
                    enddo
                  enddo
                enddo

              case (6)

                do i=imin,imax
                  do j=jmin,jmax
                    do k=kmin,kmax

                      gsuper = gmetric%grid(igxbc)%gsup(i,j,k+1,:,:)

                      v_cov(i,j,k+1,3) = -(gsuper(3,1)*v_cov(i,j,k+1,1)
     .                                    +gsuper(3,2)*v_cov(i,j,k+1,2)
     .                                    -v_cnv(i,j,k+1,3))/gsuper(3,3)

                      call transformFromCurvToCurv(i,j,k+1,igxbc,igybc,igzbc
     .               ,v_cov(i,j,k+1,1),v_cov(i,j,k+1,2),v_cov(i,j,k+1,3)
     .               ,v_cnv(i,j,k+1,1),v_cnv(i,j,k+1,2),v_cnv(i,j,k+1,3)
     .               ,.true.)
                    enddo
                  enddo
                enddo

              end select

            endif
          enddo

        enddo
      enddo

c     End program

      end subroutine synchronize_cnv

c     fillGhostNodes
c     ###############################################################
      subroutine fillGhostNodes(ieq,ivar,nvar,dim,loc,bctype
     .                         ,array,array0)

c     ---------------------------------------------------------------
c     Sets adequate boundary conditions on array.
c
c     On input:
c       * ieq    -> equation identifier
c       * ivar   -> vector component
c       * nvar   -> vector dimension
c       * dim    -> dimension (1 -> X, 2 -> Y, 3 -> Z)
c       * loc    -> location in dimension (0 -> right, 1 -> left)
c       * bctype -> type of BC (dirichlet, neumann, periodic, etc.)
c       * array  -> real array with ghost-nodes
c       * array0 -> auxiliary real array
c     ---------------------------------------------------------------

      implicit none       !For safe fortran

c     Call variables

      integer(4) :: ieq,dim,loc,bctype,nvar,ivar
      real(8)    :: array (0:nxbc+1,0:nybc+1,0:nzbc+1,nvar)
     .             ,array0(0:nxbc+1,0:nybc+1,0:nzbc+1,nvar)

c     Local variables

      integer(4) :: neq,ibc,i,j,k,ig,jg,kg

c     Begin program

      nnvar = nvar

c     Determine boundary limits

      call findLoopLimits(dim,loc,iimin,iimax,jjmin,jjmax,kkmin,kkmax
     .                           ,imin ,imax ,jmin ,jmax ,kmin ,kmax)

      if (do_nothing) return

c     Allocate rhs

      select case(dim)
      case(1)
        allocate(rhs(0:nybc+1,0:nzbc+1))
      case(2)
        allocate(rhs(0:nxbc+1,0:nzbc+1))
      case(3)
        allocate(rhs(0:nxbc+1,0:nybc+1))
      end select

c     Find BC update

      ibc = (1+loc)+2*(dim-1)

      if (nvar == 1) then
        select case(bctype)
        case(PER)
          call periodicBC(array(:,:,:,nvar),ibc)
        case(EQU)
          call dirichletBC(array(:,:,:,nvar),array0(:,:,:,nvar)
     .                    ,ieq,dim,loc,order)
        case(DIR)
          call dirichletBC(array(:,:,:,nvar),zeros,ieq,dim,loc,order)
        case(NEU)
          call neumannBC(array(:,:,:,nvar),ieq,dim,loc)
        case(EXT)
          call fillFaces(array(:,:,:,nvar),dim,loc)
        case default
          write (*,*) 'BC',bctype,' not implemented'
          stop
        end select
      else
        select case(bctype)
        case(PER)
          call periodicBC(array(:,:,:,ivar),ibc)
        case(EQU)
          call dirichletBC(ivar,array,array0,ieq,dim,loc,order)
        case(DIR)
          call dirichletBC(ivar,array,vzeros,ieq,dim,loc,order)
        case(NEU)
          call neumannBC(ivar,array,ieq,dim,loc)
        case(EXT)
          call fillFaces(array(:,:,:,ivar),dim,loc)
        case default
          write (*,*) 'BC',bctype,' not implemented'
          stop
        end select
      endif

c     Update BC ghost nodes

cc      select case (ibc)
cc      case (1)                  !x0
cc        array(0     ,jmin:jmax,kmin:kmax,ivar)= rhs(jmin:jmax,kmin:kmax)
cc      case (2)                  !x1
cc        array(nxbc+1,jmin:jmax,kmin:kmax,ivar)= rhs(jmin:jmax,kmin:kmax)
cc      case (3)                  !y0
cc        array(imin:imax,0     ,kmin:kmax,ivar)= rhs(imin:imax,kmin:kmax)
cc      case (4)                  !y1
cc        array(imin:imax,nybc+1,kmin:kmax,ivar)= rhs(imin:imax,kmin:kmax)
cc      case (5)                  !z0
cc        array(imin:imax,jmin:jmax,0     ,ivar)= rhs(imin:imax,jmin:jmax)
cc      case (6)                  !z1
cc        array(imin:imax,jmin:jmax,nzbc+1,ivar)= rhs(imin:imax,jmin:jmax)
cc      case default
cc        write (*,*) 'Boundary',ibc,' non existent'
cc        stop
cc      end select

      deallocate(rhs)

c     End

      end subroutine fillGhostNodes

c     periodicBC
c     #################################################################
      subroutine periodicBC(array,ibc)
c     -----------------------------------------------------------------
c     Imposes periodic BC. On input:
c        * array -> 3D array to impose BC's on.
c        * ibc   -> integer boundary identifier
c                   (0 -> X1, 1->X2, 2->y1, 3->y2, 4->z1, 5->z2)
c     -----------------------------------------------------------------

      implicit none

c     Call variables

      integer(4) :: ibc
      real(8)    :: array (0:nxbc+1,0:nybc+1,0:nzbc+1)

c     Local variables

c     Begin program

#if !defined(petsc)

      select case (ibc)
      case (1)                  !x0
       array(0     ,jmin:jmax,kmin:kmax)=array(nxbc,jmin:jmax,kmin:kmax)
      case (2)                  !x1                                       
       array(nxbc+1,jmin:jmax,kmin:kmax)=array(1   ,jmin:jmax,kmin:kmax)
      case (3)                  !y0                                       
       array(imin:imax,0     ,kmin:kmax)=array(imin:imax,nybc,kmin:kmax)
      case (4)                  !y1                                       
       array(imin:imax,nybc+1,kmin:kmax)=array(imin:imax,1   ,kmin:kmax)
      case (5)                  !z0                                       
       array(imin:imax,jmin:jmax,0     )=array(imin:imax,jmin:jmax,nzbc)
      case (6)                  !z1                                       
       array(imin:imax,jmin:jmax,nzbc+1)=array(imin:imax,jmin:jmax,1   )
      case default
       write (*,*) 'Boundary',ibc,' non existent'
       stop
      end select

cc      select case (ibc)
cc      case (1)
cc        rhs(jmin:jmax,kmin:kmax) = array(nxbc,jmin:jmax,kmin:kmax)
cc      case (2)
cc        rhs(jmin:jmax,kmin:kmax) = array(1   ,jmin:jmax,kmin:kmax)
cc      case (3)
cc        rhs(imin:imax,kmin:kmax) = array(imin:imax,nybc,kmin:kmax)
cc      case (4)
cc        rhs(imin:imax,kmin:kmax) = array(imin:imax,1   ,kmin:kmax)
cc      case (5)
cc        rhs(imin:imax,jmin:jmax) = array(imin:imax,jmin:jmax,nzbc)
cc      case (6)
cc        rhs(imin:imax,jmin:jmax) = array(imin:imax,jmin:jmax,1   )
cc      end select
cc
cc#else !Do nothing
cc
cc      select case (ibc)
cc      case (1)
cc        rhs(jmin:jmax,kmin:kmax) = array(0     ,jmin:jmax,kmin:kmax)
cc      case (2)
cc        rhs(jmin:jmax,kmin:kmax) = array(nxbc+1,jmin:jmax,kmin:kmax)
cc      case (3)
cc        rhs(imin:imax,kmin:kmax) = array(imin:imax,0     ,kmin:kmax)
cc      case (4)
cc        rhs(imin:imax,kmin:kmax) = array(imin:imax,nybc+1,kmin:kmax)
cc      case (5)
cc        rhs(imin:imax,jmin:jmax) = array(imin:imax,jmin:jmax,0     )
cc      case (6)
cc        rhs(imin:imax,jmin:jmax) = array(imin:imax,jmin:jmax,nzbc+1)
cc      end select

#endif

c     End program

      end subroutine periodicBC

c     fillFaces
c     #################################################################
      subroutine fillFaces(array,dim,loc)
c     -----------------------------------------------------------------
c     Fills ghost cells in faces by extrapolation (6 in 3D). On input:
c        * array  -> 3D array to impose BC's on.
c        * dim    -> dimension (1 -> X, 2 -> Y, 3 -> Z)
c        * loc    -> location in dimension (0 -> right, 1 -> left)
c     This routine is only used when no good information is available
c     to fill the ghost cells of a particular quantity (for instance,
c     during MG restriction of arrays that contain EQU BCs).
c     -----------------------------------------------------------------

      implicit none

c     Call variables

      integer(4) :: dim,loc
      real(8)    :: array (0:nxbc+1,0:nybc+1,0:nzbc+1)

c     Local variables

      integer(4) :: ito,jto,kto,ifrom,jfrom,kfrom
     .             ,imn,imx,jmn,jmx,kmn,kmx,ibc,bcond(6)

c     Begin program

      ibc = (1+loc)+2*(dim-1)

      bcond = 0
      bcond(ibc) = EXT

c     Determine boundary limits

      call findLoopLimits(dim,loc,iimin,iimax,jjmin,jjmax,kkmin,kkmax
     .                           ,imin ,imax ,jmin ,jmax ,kmin ,kmax)

      if (do_nothing) return

c     Select edge

      select case (ibc)
      case (1) !On face #1, fill face i=0

        if (imin == 1) then

          ito   = imin-1
          ifrom = min(imin+order,nxbc)

          imn   = min(ito,ifrom)
          imx   = max(ito,ifrom)

          do j=jmin,jmax
            do k=kmin,kmax
              call extrapolateToCorner(ibc,ito,j,k,imn,j,k
     .                  ,imx,j,k,array(imn:imx,j,k),bcond)
            enddo
          enddo

        endif

      case (2) !On face #2, fill face i=nx+1

        if (imax == nxbc) then

          ito   = imax+1
          ifrom = max(imax-order,1)

          imn   = min(ito,ifrom)
          imx   = max(ito,ifrom)

          do j=jmin,jmax
            do k=kmin,kmax
              call extrapolateToCorner(ibc,ito,j,k,imn,j,k
     .                  ,imx,j,k,array(imn:imx,j,k),bcond)
            enddo
          enddo

        endif

      case (3) !On face #3, fill face j=0

        if (jmin == 1) then

          jto   = jmin-1
          jfrom = min(jmin+order,nybc)

          jmn   = min(jto,jfrom)
          jmx   = max(jto,jfrom)

          do k=kmin,kmax
            do i=imin,imax
                call extrapolateToCorner(ibc,i,jto,k,i,jmn,k
     .                  ,i,jmx,k,array(i,jmn:jmx,k),bcond)
            enddo
          enddo

        endif

      case (4) !On face #4, fill face j=ny+1

        if (jmax == nybc) then

          jto   = jmax+1
          jfrom = max(jmax-order,1)

          jmn   = min(jto,jfrom)
          jmx   = max(jto,jfrom)

          do k=kmin,kmax
            do i=imin,imax
              call extrapolateToCorner(ibc,i,jto,k,i,jmn,k
     .                  ,i,jmx,k,array(i,jmn:jmx,k),bcond)
            enddo
          enddo

        endif

      case (5) !On face #5, fill face k=0

        if (kmin == 1) then

          kto   = kmin-1
          kfrom = min(kmin+order,nzbc)

          kmn   = min(kto,kfrom)
          kmx   = max(kto,kfrom)

          do i=imin,imax
            do j=jmin,jmax
              call extrapolateToCorner(ibc,i,j,kto,i,j,kmn
     .                  ,i,j,kmx,array(i,j,kmn:kmx),bcond)
            enddo
          enddo

        endif

      case (6) !On face #6, fill face k=nz+1

        if (kmax == nzbc) then

          kto   = kmax+1
          kfrom = max(kmax-order,1)

          kmn   = min(kto,kfrom)
          kmx   = max(kto,kfrom)

          do i=imin,imax
            do j=jmin,jmax
              call extrapolateToCorner(ibc,i,j,kto,i,j,kmn
     .                  ,i,j,kmx,array(i,j,kmn:kmx),bcond)
            enddo
          enddo

        endif

      end select

c     End program

      end subroutine fillFaces

c     fillEdges
c     #################################################################
      subroutine fillEdges(array,array0,dim,loc,bcond)
c     -----------------------------------------------------------------
c     Fills ghost cells in edges (12 in 3D). On input:
c        * array  -> 3D array to impose BC's on.
c        * array0 -> 3D array containing equilibrium BCs
c        * dim    -> dimension (1 -> X, 2 -> Y, 3 -> Z)
c        * loc    -> location in dimension (0 -> right, 1 -> left)
c     -----------------------------------------------------------------

      implicit none

c     Call variables

      integer(4) :: dim,loc,bcond(6)
      real(8)    :: array (0:nxbc+1,0:nybc+1,0:nzbc+1)
     .             ,array0(0:nxbc+1,0:nybc+1,0:nzbc+1)

c     Local variables

      integer(4) :: ito,jto,kto,ifrom,jfrom,kfrom
     .             ,imn,imx,jmn,jmx,kmn,kmx,ibc

c     Begin program

      ibc = (1+loc)+2*(dim-1)

c     Determine boundary limits

      call findLoopLimits(dim,loc,iimin,iimax,jjmin,jjmax,kkmin,kkmax
     .                           ,imin ,imax ,jmin ,jmax ,kmin ,kmax)

      if (do_nothing) return

c     Select edge

      select case (ibc)
      case (1) !On face #1, fill edges y=jmin-1,y=jmax+1, vary z

        if (imin == 1) then
          ito   = imin-1
          ifrom = min(imin+order,nxbc)

          imn   = min(ito,ifrom)
          imx   = max(ito,ifrom)

          if (jmin == 1) then
            jto   = jmin-1
            jfrom = min(jmin+order,nybc)

            jmn   = min(jto,jfrom)
            jmx   = max(jto,jfrom)

            do k=kmin,kmax
              if ((bcond(1)==EQU.or.bcond(3)==EQU) .and. order == 0)then
                array(ito,jto,k) = array0(ito,jto,k)
              else
                call extrapolateToCorner(ibc,ito,jto,k,imn,jmn,k
     .                  ,imx,jmx,k,array(imn:imx,jmn:jmx,k),bcond)
              endif
            enddo
          endif

          if (jmax == nybc) then
            jto   = jmax+1
            jfrom = max(jmax-order,1)

            jmn   = min(jto,jfrom)
            jmx   = max(jto,jfrom)

            do k=kmin,kmax
              if ((bcond(1)==EQU.or.bcond(4)==EQU).and. order == 0)then
                array(ito,jto,k) = array0(ito,jto,k)
              else
                call extrapolateToCorner(ibc,ito,jto,k,imn,jmn,k
     .                  ,imx,jmx,k,array(imn:imx,jmn:jmx,k),bcond)
              endif
            enddo
          endif

        endif

      case (2) !On face #2, fill edges y=jmin-1,y=jmax+1, vary z

        if (imax == nxbc) then
          ito   = imax+1
          ifrom = max(imax-order,1)

          imn   = min(ito,ifrom)
          imx   = max(ito,ifrom)

          if (jmin == 1) then
            jto   = jmin-1
            jfrom = min(jmin+order,nybc)

            jmn   = min(jto,jfrom)
            jmx   = max(jto,jfrom)

            do k=kmin,kmax
              if ((bcond(2)==EQU.or.bcond(3)==EQU).and. order == 0)then
                array(ito,jto,k) = array0(ito,jto,k)
              else
                call extrapolateToCorner(ibc,ito,jto,k,imn,jmn,k
     .                  ,imx,jmx,k,array(imn:imx,jmn:jmx,k),bcond)
              endif
            enddo
          endif

          if (jmax == nybc) then
            jto   = jmax+1
            jfrom = max(jmax-order,1)

            jmn   = min(jto,jfrom)
            jmx   = max(jto,jfrom)

            do k=kmin,kmax
              if ((bcond(2)==EQU.or.bcond(4)==EQU).and. order == 0)then
                array(ito,jto,k) = array0(ito,jto,k)
              else
                call extrapolateToCorner(ibc,ito,jto,k,imn,jmn,k
     .                  ,imx,jmx,k,array(imn:imx,jmn:jmx,k),bcond)
              endif
            enddo
          endif

        endif

      case (3) !On face #3, fill edges z=kmin-1,z=kmax+1, vary x

        if (jmin == 1) then

          jto   = jmin-1
          jfrom = min(jmin+order,nybc)

          jmn   = min(jto,jfrom)
          jmx   = max(jto,jfrom)

          if (kmin == 1) then

            kto   = kmin-1
            kfrom = min(kmin+order,nzbc)

            kmn   = min(kto,kfrom)
            kmx   = max(kto,kfrom)

            do i=imin,imax
              if ((bcond(3)==EQU.or.bcond(5)==EQU).and. order == 0)then
                array(i,jto,kto) = array0(i,jto,kto)
              else
                call extrapolateToCorner(ibc,i,jto,kto,i,jmn,kmn
     .                  ,i,jmx,kmx,array(i,jmn:jmx,kmn:kmx),bcond)
              endif
            enddo
          endif

          if (kmax == nzbc) then

            kto   = kmax+1
            kfrom = max(kmax-order,1)

            kmn   = min(kto,kfrom)
            kmx   = max(kto,kfrom)

            do i=imin,imax
              if ((bcond(3)==EQU.or.bcond(6)==EQU).and. order == 0)then
                array(i,jto,kto) = array0(i,jto,kto)
              else
                call extrapolateToCorner(ibc,i,jto,kto,i,jmn,kmn
     .                  ,i,jmx,kmx,array(i,jmn:jmx,kmn:kmx),bcond)
              endif
            enddo
          endif

        endif

      case (4) !On face #4, fill edges z=kmin-1,z=kmax+1, vary x

        if (jmax == nybc) then

          jto   = jmax+1
          jfrom = max(jmax-order,1)

          jmn   = min(jto,jfrom)
          jmx   = max(jto,jfrom)

          if (kmin == 1) then

            kto   = kmin-1
            kfrom = min(kmin+order,nzbc)

            kmn   = min(kto,kfrom)
            kmx   = max(kto,kfrom)

            do i=imin,imax
              if ((bcond(4)==EQU.or.bcond(5)==EQU).and. order == 0)then
                array(i,jto,kto) = array0(i,jto,kto)
              else
                call extrapolateToCorner(ibc,i,jto,kto,i,jmn,kmn
     .                  ,i,jmx,kmx,array(i,jmn:jmx,kmn:kmx),bcond)
              endif
            enddo
          endif

          if (kmax == nzbc) then

            kto   = kmax+1
            kfrom = max(kmax-order,1)

            kmn   = min(kto,kfrom)
            kmx   = max(kto,kfrom)

            do i=imin,imax
              if ((bcond(4)==EQU.or.bcond(6)==EQU).and. order == 0)then
                array(i,jto,kto) = array0(i,jto,kto)
              else
                call extrapolateToCorner(ibc,i,jto,kto,i,jmn,kmn
     .                  ,i,jmx,kmx,array(i,jmn:jmx,kmn:kmx),bcond)
              endif
            enddo
          endif

        endif

      case (5) !On face #5, fill edges x=imin-1,x=imax+1, vary y

        if (kmin == 1) then

          kto   = kmin-1
          kfrom = min(kmin+order,nzbc)

          kmn   = min(kto,kfrom)
          kmx   = max(kto,kfrom)

          if (imin == 1) then

            ito   = imin-1
            ifrom = min(imin+order,nxbc)

            imn   = min(ito,ifrom)
            imx   = max(ito,ifrom)

            do j=jmin,jmax
              if ((bcond(5)==EQU.or.bcond(1)==EQU).and. order == 0)then
                array(ito,j,kto) = array0(ito,j,kto)
              else
                call extrapolateToCorner(ibc,ito,j,kto,imn,j,kmn
     .                  ,imx,j,kmx,array(imn:imx,j,kmn:kmx),bcond)
              endif
            enddo
          endif

          if (imax == nxbc) then

            ito   = imax+1
            ifrom = max(imax-order,1)

            imn   = min(ito,ifrom)
            imx   = max(ito,ifrom)

            do j=jmin,jmax
              if ((bcond(5)==EQU.or.bcond(2)==EQU).and. order == 0)then
                array(ito,j,kto) = array0(ito,j,kto)
              else
                call extrapolateToCorner(ibc,ito,j,kto,imn,j,kmn
     .                  ,imx,j,kmx,array(imn:imx,j,kmn:kmx),bcond)
              endif
            enddo
          endif

        endif

      case (6) !On face #6, fill edges x=imin-1,x=imax+1, vary y

        if (kmax == nzbc) then

          kto   = kmax+1
          kfrom = max(kmax-order,1)

          kmn   = min(kto,kfrom)
          kmx   = max(kto,kfrom)

          if (imin == 1) then

            ito   = imin-1
            ifrom = min(imin+order,nxbc)

            imn   = min(ito,ifrom)
            imx   = max(ito,ifrom)

            do j=jmin,jmax
              if ((bcond(6)==EQU.or.bcond(1)==EQU).and. order == 0)then
                array(ito,j,kto) = array0(ito,j,kto)
              else
                call extrapolateToCorner(ibc,ito,j,kto,imn,j,kmn
     .                  ,imx,j,kmx,array(imn:imx,j,kmn:kmx),bcond)
              endif
            enddo
          endif

          if (imax == nxbc) then

            ito   = imax+1
            ifrom = max(imax-order,1)

            imn   = min(ito,ifrom)
            imx   = max(ito,ifrom)

            do j=jmin,jmax
              if ((bcond(6)==EQU.or.bcond(2)==EQU).and. order == 0)then
                array(ito,j,kto) = array0(ito,j,kto)
              else
                call extrapolateToCorner(ibc,ito,j,kto,imn,j,kmn
     .                  ,imx,j,kmx,array(imn:imx,j,kmn:kmx),bcond)
              endif
            enddo
          endif

        endif

      end select

c     End program

      end subroutine fillEdges

c     fillCorners
c     #################################################################
      subroutine fillCorners(array,array0,dim,loc,bcond)
c     -----------------------------------------------------------------
c     Imposes ghost cells in corners (8 in 3D). On input:
c        * array -> 3D array to impose BC's on.
c        * array0 -> 3D array containing equilibrium BCs
c        * dim    -> dimension (1 -> X, 2 -> Y, 3 -> Z)
c        * loc    -> location in dimension (0 -> right, 1 -> left)
c     -----------------------------------------------------------------

      implicit none

c     Call variables

      integer(4) :: dim,loc,bcond(6)
      real(8)    :: array (0:nxbc+1,0:nybc+1,0:nzbc+1)
     .             ,array0(0:nxbc+1,0:nybc+1,0:nzbc+1)

c     Local variables

      integer(4) :: ito,jto,kto,ifrom,jfrom,kfrom
     .             ,imn,imx,jmn,jmx,kmn,kmx,ibc

c     Begin program

      ibc = (1+loc)+2*(dim-1)

c     Determine boundary limits

      call findLoopLimits(dim,loc,iimin,iimax,jjmin,jjmax,kkmin,kkmax
     .                           ,imin ,imax ,jmin ,jmax ,kmin ,kmax)

      if (do_nothing) return

c     Select edge

      if     (imin == 1  .and. jmin == 1  .and. kmin == 1) then
        ifrom = min(imin+order,nxbc)
        jfrom = min(jmin+order,nybc)
        kfrom = min(kmin+order,nzbc)
        ito   = imin-1
        jto   = jmin-1
        kto   = kmin-1

        imn   = min(ito,ifrom)
        jmn   = min(jto,jfrom)
        kmn   = min(kto,kfrom)
        imx   = max(ito,ifrom)
        jmx   = max(jto,jfrom)
        kmx   = max(kto,kfrom)

        call extrapolateToCorner(ibc,ito,jto,kto,imn,jmn,kmn
     .          ,imx,jmx,kmx,array(imn:imx,jmn:jmx,kmn:kmx),bcond)
      endif

      if (imin == 1  .and. jmax == nybc .and. kmin == 1) then

        ifrom = min(imin+order,nxbc)
        jfrom = max(jmax-order,1   )
        kfrom = min(kmin+order,nzbc)
        ito   = imin-1
        jto   = jmax+1
        kto   = kmin-1

        imn   = min(ito,ifrom)
        jmn   = min(jto,jfrom)
        kmn   = min(kto,kfrom)
        imx   = max(ito,ifrom)
        jmx   = max(jto,jfrom)
        kmx   = max(kto,kfrom)

        call extrapolateToCorner(ibc,ito,jto,kto,imn,jmn,kmn
     .           ,imx,jmx,kmx,array(imn:imx,jmn:jmx,kmn:kmx),bcond)
      endif

      if (imax == nxbc .and. jmin == 1  .and. kmin == 1) then

        ifrom = max(imax-order,1   )
        jfrom = min(jmin+order,nybc)
        kfrom = min(kmin+order,nzbc)
        ito   = imax+1
        jto   = jmin-1
        kto   = kmin-1

        imn   = min(ito,ifrom)
        jmn   = min(jto,jfrom)
        kmn   = min(kto,kfrom)
        imx   = max(ito,ifrom)
        jmx   = max(jto,jfrom)
        kmx   = max(kto,kfrom)

        call extrapolateToCorner(ibc,ito,jto,kto,imn,jmn,kmn
     .            ,imx,jmx,kmx,array(imn:imx,jmn:jmx,kmn:kmx),bcond)
      endif

      if (imax == nxbc .and. jmax == nybc .and. kmin == 1) then

        ifrom = max(imax-order,1   )
        jfrom = max(jmax-order,1   )
        kfrom = min(kmin+order,nzbc)
        ito   = imax+1
        jto   = jmax+1
        kto   = kmin-1

        imn   = min(ito,ifrom)
        jmn   = min(jto,jfrom)
        kmn   = min(kto,kfrom)
        imx   = max(ito,ifrom)
        jmx   = max(jto,jfrom)
        kmx   = max(kto,kfrom)

        call extrapolateToCorner(ibc,ito,jto,kto,imn,jmn,kmn
     .            ,imx,jmx,kmx,array(imn:imx,jmn:jmx,kmn:kmx),bcond)
      endif

      if (imin == 1  .and. jmin == 1  .and. kmax == nzbc) then

        ifrom = min(imin+order,nxbc)
        jfrom = min(jmin+order,nybc)
        kfrom = max(kmax-order,1   )
        ito   = imin-1
        jto   = jmin-1
        kto   = kmax+1

        imn   = min(ito,ifrom)
        jmn   = min(jto,jfrom)
        kmn   = min(kto,kfrom)
        imx   = max(ito,ifrom)
        jmx   = max(jto,jfrom)
        kmx   = max(kto,kfrom)

        call extrapolateToCorner(ibc,ito,jto,kto,imn,jmn,kmn
     .            ,imx,jmx,kmx,array(imn:imx,jmn:jmx,kmn:kmx),bcond)
      endif

      if (imin == 1  .and. jmax == nybc .and. kmax == nzbc) then

        ifrom = min(imin+order,nxbc)
        jfrom = max(jmax-order,1   )
        kfrom = max(kmax-order,1   )
        ito   = imin-1
        jto   = jmax+1
        kto   = kmax+1

        imn   = min(ito,ifrom)
        jmn   = min(jto,jfrom)
        kmn   = min(kto,kfrom)
        imx   = max(ito,ifrom)
        jmx   = max(jto,jfrom)
        kmx   = max(kto,kfrom)

        call extrapolateToCorner(ibc,ito,jto,kto,imn,jmn,kmn
     .            ,imx,jmx,kmx,array(imn:imx,jmn:jmx,kmn:kmx),bcond)
      endif

      if (imax == nxbc .and. jmin == 1  .and. kmax == nzbc) then

        ifrom = max(imax-order,1   )
        jfrom = min(jmin+order,nybc)
        kfrom = max(kmax-order,1   )
        ito   = imax+1
        jto   = jmin-1
        kto   = kmax+1

        imn   = min(ito,ifrom)
        jmn   = min(jto,jfrom)
        kmn   = min(kto,kfrom)
        imx   = max(ito,ifrom)
        jmx   = max(jto,jfrom)
        kmx   = max(kto,kfrom)

        call extrapolateToCorner(ibc,ito,jto,kto,imn,jmn,kmn
     .            ,imx,jmx,kmx,array(imn:imx,jmn:jmx,kmn:kmx),bcond)
      endif

      if (imax == nxbc .and. jmax == nybc .and. kmax == nzbc) then

        ifrom = max(imax-order,1 )
        jfrom = max(jmax-order,1 )
        kfrom = max(kmax-order,1 )
        ito   = imax+1
        jto   = jmax+1
        kto   = kmax+1

        imn   = min(ito,ifrom)
        jmn   = min(jto,jfrom)
        kmn   = min(kto,kfrom)
        imx   = max(ito,ifrom)
        jmx   = max(jto,jfrom)
        kmx   = max(kto,kfrom)

        call extrapolateToCorner(ibc,ito,jto,kto,imn,jmn,kmn
     .            ,imx,jmx,kmx,array(imn:imx,jmn:jmx,kmn:kmx),bcond)

      endif

c     End program

      end subroutine fillCorners

c     extrapolateToCorner
c     #################################################################
      subroutine extrapolateToCorner(ibc,ito,jto,kto,imn,jmn,kmn
     .                              ,imx,jmx,kmx,array,bcond)
c     -----------------------------------------------------------------
c     Routine that extrapolates information to edge or vertex corners.
c     It is the workhorse of the fillEdge and fillCorner routines.
c     In the call sequence:
c        * ibc: boundary identifier (integer)
c        * ito,jto,kto: coordinates of ghost cell to be filled (integer)
c        * imn,jmn,kmn: node coordinates identifying start of
c                       interpolation region (integer)
c        * imx,jmx,kmx: node coordinates identifying end of interpolation
c                       region (integer)
c        * array: local array defined on interpolation region, including
c                 ghost cell (ito,jto,kto) (real)
c     The order of extrapolation 'order' is defined in the module header.
c     -----------------------------------------------------------------

      implicit none

c     Call variables

      integer(4) :: ibc,ito,jto,kto,imn,jmn,kmn,imx,jmx,kmx,bcond(6)
      real(8)    :: array(imn:imx,jmn:jmx,kmn:kmx)

c     Local variables

      integer(4) :: i,j,k,ig,jg,kg,inx,iny,inz
     .             ,orderx,ordery,orderz,offset
     .             ,imni,imxi,jmni,jmxi,kmni,kmxi,num

      real(8)    :: valx(1),valy(1),valz(1),xto(1),yto(1),zto(1)

      logical    :: flgx,flgy,flgz

c     Extrapolation

      real(8),allocatable,dimension(:) :: xx,yy,zz

c     Begin program

c     Identify directions for extrapolation.

      inx = imx-imn
      iny = jmx-jmn
      inz = kmx-kmn

      flgx = .true.
      flgy = .true.
      flgz = .true.

      !Do not interpolate **along** direction if
      !   1) not enough points,
      if (inx == 0) flgx = .false.
      if (iny == 0) flgy = .false.
      if (inz == 0) flgz = .false.

      !Interpolate **across** boundary if
      !   1) periodic boundary
      if (bcond(1) == PER .and. flgx) then
        flgy = .false.
        flgz = .false.
      endif
      if (bcond(3) == PER .and. flgy) then
        flgx = .false.
        flgz = .false.
      endif
      if (bcond(5) == PER .and. flgz) then
        flgx = .false.
        flgy = .false.
      endif

cc      select case(ibc)
cc      case(1,2)
cc        if (bcond(ibc) == PER .and. flgx) then
cc          flgy = .false.
cc          flgz = .false.
cc        endif
cc      case(3,4)
cc        if (bcond(ibc) == PER .and. flgy) then
cc          flgx = .false.
cc          flgz = .false.
cc        endif
cc      case(5,6)
cc        if (bcond(ibc) == PER .and. flgz) then
cc          flgx = .false.
cc          flgy = .false.
cc        endif
cc      end select

      !Interpolate **along** boundary if
      !   1) equilibrium boundary
      select case(ibc)
      case(1,2)
        if (bcond(ibc) == EQU .and. (flgy.or.flgz)) then
          flgx = .false.
        endif
      case(3,4)
        if (bcond(ibc) == EQU .and. (flgx.or.flgz)) then
          flgy = .false.
        endif
      case(5,6)
        if (bcond(ibc) == EQU .and. (flgx.or.flgy)) then
          flgz = .false.
        endif
      end select

      !Check some direction for interpolation exists
      if ((.not.flgx).and.(.not.flgy).and.(.not.flgz)) then
        write (*,*) 'Error in extrapolateToCorner:'
        write (*,*) 'No extrapolation directions left'
        write (*,*) 'Aborting...'
        stop
      endif

c     Allocate spatial coordinate vectors.

      allocate(xx(inx),yy(iny),zz(inz))

c     Initialize accumulators

      valx = 0d0
      valy = 0d0
      valz = 0d0

      num = 0

c     X extrapolation

      if (flgx) then
        if (ito == imn) then
          call getMGmap(imn+1,jmn,kmn,igxbc,igybc,igzbc,ig,jg,kg)
          xx = grid_params%xx(ig:ig+inx-1)
          imni = imn+1
          imxi = imx
        else
          call getMGmap(imx-1,jmn,kmn,igxbc,igybc,igzbc,ig,jg,kg)
          xx = grid_params%xx(ig-(inx-1):ig)
          imni = imn
          imxi = imx-1
        endif

        call getMGmap(ito,jto,kto,igxbc,igybc,igzbc,ig,jg,kg)
        xto = grid_params%xx(ig)

        orderx = min(order,inx-1)
        call IntDriver1d(inx,xx,array(imni:imxi,jto,kto)
     .                  ,1,xto,valx,orderx)

        num = num + 1
      endif

c     Y extrapolation

      if (flgy) then

        if (jto == jmn) then
          call getMGmap(imn,jmn+1,kmn,igxbc,igybc,igzbc,ig,jg,kg)
          yy = grid_params%yy(jg:jg+iny-1)
          jmni = jmn+1
          jmxi = jmx
        else
          call getMGmap(imn,jmx-1,kmn,igxbc,igybc,igzbc,ig,jg,kg)
          yy = grid_params%yy(jg-(iny-1):jg)
          jmni = jmn
          jmxi = jmx-1
        endif

        call getMGmap(ito,jto,kto,igxbc,igybc,igzbc,ig,jg,kg)
        yto = grid_params%yy(jg)

        ordery = min(order,iny-1)
        call IntDriver1d(iny,yy,array(ito,jmni:jmxi,kto)
     .                  ,1,yto,valy,ordery)
        
        num = num + 1

      endif

c     Z extrapolation

      if (flgz) then

        if (kto == kmn) then
          call getMGmap(imn,jmn,kmn+1,igxbc,igybc,igzbc,ig,jg,kg)
          zz = grid_params%zz(kg:kg+inz-1)
          kmni = kmn+1
          kmxi = kmx
        else
          call getMGmap(imn,jmn,kmx-1,igxbc,igybc,igzbc,ig,jg,kg)
          zz = grid_params%zz(kg-(inz-1):kg)
          kmni = kmn
          kmxi = kmx-1
        endif

        call getMGmap(ito,jto,kto,igxbc,igybc,igzbc,ig,jg,kg)
        zto = grid_params%zz(kg)

        orderz = min(order,inz-1)
        call IntDriver1d(inz,zz,array(ito,jto,kmni:kmxi)
     .                  ,1,zto,valz,orderz)

        num = num + 1
      endif

c     Average extrapolation results and fill ghost cell

      array(ito,jto,kto) = (valx(1) + valy(1) + valz(1))/num

c     End program

      deallocate(xx,yy,zz)

      end subroutine extrapolateToCorner

      end module imposeBCinterface
