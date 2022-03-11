c module bc_ops
c #####################################################################
      module bc_ops

c     -----------------------------------------------------------------
c     This module defines various BC-related operators
c     -----------------------------------------------------------------

        use grid_def

      contains

c     isBdry
c     #################################################################
      function isBdry(gl_def,i,igr,ibc)

c     -----------------------------------------------------------------
c     Detects whether we are at a PHYSICAL boundary or not. Boundary is
c     identified by ibc, with the following convention:
c        * ibc = 1, 2 (xmin,xmax)
c        * ibc = 3, 4 (ymin,ymax)
c        * ibc = 5, 6 (zmin,zmax)
c     -----------------------------------------------------------------

        implicit none

c     Call variables

        integer :: i,igr,ibc
        logical :: isBdry

        type(grid_mg_def),pointer :: gl_def

c     Local variables

        integer :: ilog,itst,diml,locl

c     Begin program

cc        locl = abs(mod(ibc,2)-2)-1
cc        diml = 1 + (ibc-1-locl)/2
        diml = (ibc+1)/2
        locl = ibc + 1 - 2*diml

        isBdry = isBdry2(gl_def,i,igr,diml,locl)

      end function isBdry

c     isBdry2
c     #################################################################
      function isBdry2(gl_def,i,igr,diml,locl)

c     -----------------------------------------------------------------
c     Detects whether we are at a PHYSICAL boundary or not. Boundary is
c     identified by dim, loc.
c     -----------------------------------------------------------------

        implicit none

c     Call variables

        integer :: i,igr
        logical :: isBdry2

        integer :: diml,locl

        type(grid_mg_def),pointer :: gl_def

c     Local variables

        integer :: ilog,itst

c     Begin program

        itst = 1

        select case(diml)
        case(1)
          ilog = gl_def%ilo(igr)
          if (locl == 1) itst = gl_def%nxgl(igr)
        case(2)
          ilog = gl_def%jlo(igr)
          if (locl == 1) itst = gl_def%nygl(igr)
        case(3)
          ilog = gl_def%klo(igr)
          if (locl == 1) itst = gl_def%nzgl(igr)
        end select

        isBdry2 = (i+ilog-1 == itst)

      end function isBdry2

c     isBC
c     #################################################################
      function isBC(g_def,i,igr,ibc,BC,bcnd)

c     -----------------------------------------------------------------
c     Checks whether point i at grid level igr, boundary location
c     ibc=(1,6) (see isBdry for definition) is of boundary type BC.
c     -----------------------------------------------------------------

        implicit none

c     Call variables

        integer :: i,igr,ibc,BC
        logical :: isBC
        integer,optional :: bcnd(6)

        type(grid_mg_def),pointer :: g_def

c     Local variables

        integer :: bcc(6)

c     Begin program

        if (ibc < 1 .or. ibc > 6) then
          isBC = .false.
          return
        endif

        if (PRESENT(bcnd)) then
          bcc = bcnd
        else
          bcc = bcond
        endif

        isBC = (isBdry(g_def,i,igr,ibc).and.abs(bcc(ibc)) == BC)

      end function isBC

c     isSP
c     #################################################################
      function isSP(gl_def,i,j,k,igx,igy,igz) result(sing_point)

        implicit none

        integer    :: i,j,k,igx,igy,igz
        logical    :: sing_point

        type(grid_mg_def),pointer :: gl_def

c     Local variables

c     Begin program

        if (.not.associated(gl_def)) then
          write (*,*) 'Grid def struct not associated'
          write (*,*) 'Aborting in routine isSP'
          stop
        endif

        sing_point = (i+gl_def%ilo(igx)-1==1).and.(abs(bcond(1)) == SP)

      end function isSP

c$$$c     isSP2
c$$$c     #################################################################
c$$$      function isSP2(g_def,i,igr,dim,loc,ibc)
c$$$
c$$$c     -----------------------------------------------------------------
c$$$c     Checks whether point at grid level igr, dimension dim(=1,2,3)
c$$$c     and location loc(=0,1) is at SP boundary.
c$$$c     -----------------------------------------------------------------
c$$$
c$$$        implicit none
c$$$
c$$$        integer    :: i,igr
c$$$        logical    :: isSP2
c$$$        integer,optional :: dim,loc,ibc
c$$$
c$$$        type(grid_mg_def),pointer :: g_def
c$$$
c$$$        integer    :: ibcl,ilog,itst,diml,locl
c$$$
c$$$        if (PRESENT(dim) .and.PRESENT(loc)) then
c$$$          diml = dim
c$$$          locl = loc
c$$$          ibcl = (1+locl)+2*(diml-1)
c$$$        elseif (PRESENT(ibc)) then
c$$$          ibcl = ibc
c$$$        else
c$$$          ibcl = 1
c$$$        endif
c$$$
c$$$        if (ibcl < 1 .or. ibcl > 6) then
c$$$          isSP2 = .false.
c$$$          return
c$$$        else
c$$$          isSP2 = isBC(g_def,i,igr,ibcl,SP)
c$$$        endif
c$$$
c$$$      end function isSP2

c     isSP2
c     #################################################################
      function isSP2(g_def,i,igr)

c     -----------------------------------------------------------------
c     Checks whether point at grid level igr, dimension dim(=1,2,3)
c     and location loc(=0,1) is at SP boundary.
c     -----------------------------------------------------------------

        implicit none

        integer    :: i,igr
        logical    :: isSP2

        type(grid_mg_def),pointer :: g_def

        isSP2 = (i+g_def%ilo(igr)-1==1).and.(abs(bcond(1)) == SP)

      end function isSP2

c     isSYM
c     #################################################################
      function isSYM(g_def,i,igr,dim,loc,ibc)

c     -----------------------------------------------------------------
c     Checks whether point at grid level igr, dimension dim(=1,2,3)
c     and location loc(=0,1) is at SYM boundary
c     -----------------------------------------------------------------

        implicit none

        type(grid_mg_def),pointer :: g_def

        integer    :: i,igr
        logical    :: isSYM
        integer,optional :: dim,loc,ibc

        integer    :: ibcl,ilog,itst,diml,locl

        if (PRESENT(dim) .and.PRESENT(loc)) then
          diml = dim
          locl = loc
          ibcl = (1+locl)+2*(diml-1)
        elseif (PRESENT(ibc)) then
          ibcl = ibc
cc          locl = abs(mod(ibcl,2)-2)-1
cc          diml = 1 + (ibcl-1-locl)/2
          diml = (ibcl+1)/2
          locl =  ibcl+1 - 2*diml
        else
          write (*,*) 'Error in isSYM'
          write (*,*) 'Aborting...'
          stop
        endif

        isSYM = isBC(g_def,i,igr,ibcl,SYM)

      end function isSYM

c$$$c     isSYMBdry
c$$$c     #################################################################
c$$$      function isSYMBdry(gl_def,dim,loc,ibc)
c$$$
c$$$c     -----------------------------------------------------------------
c$$$c     Checks whether boundary at dimension dim(=1,2,3)
c$$$c     and location loc(=0,1) is a SYM boundary
c$$$c     -----------------------------------------------------------------
c$$$
c$$$        implicit none
c$$$
c$$$        logical    :: isSYMBdry
c$$$        integer,optional :: dim,loc,ibc
c$$$
c$$$        integer    :: ibcl,ilog,itst,diml,locl,igr
c$$$
c$$$        type(grid_mg_def),pointer :: gl_def
c$$$
c$$$        if (PRESENT(dim) .and.PRESENT(loc)) then
c$$$          diml = dim
c$$$          locl = loc
c$$$          ibcl = (1+locl)+2*(diml-1)
c$$$        elseif (PRESENT(ibc)) then
c$$$          ibcl = ibc
c$$$          locl = abs(mod(ibcl,2)-2)-1
c$$$          diml = 1 + (ibcl-1-locl)/2
c$$$        else
c$$$          write (*,*) 'Error in isSYMBdry'
c$$$          write (*,*) 'Aborting...'
c$$$          stop
c$$$        endif
c$$$
c$$$        igr = 1 !Use fine grid for test
c$$$
c$$$        select case(diml)
c$$$        case(1)
c$$$          select case(locl)
c$$$          case (0)
c$$$            itst = 1
c$$$          case (1)
c$$$            itst = gl_def%nxv(igr)
c$$$          end select
c$$$        case(2)
c$$$          select case(locl)
c$$$          case (0)
c$$$            itst = 1
c$$$          case (1)
c$$$            itst = gl_def%nyv(igr)
c$$$          end select
c$$$        case(3)
c$$$          select case(locl)
c$$$          case (0)
c$$$            itst = 1
c$$$          case (1)
c$$$            itst = gl_def%nzv(igr)
c$$$          end select
c$$$        end select
c$$$
c$$$        isSYMBdry = isBC(gl_def,itst,igr,ibcl,SYM)
c$$$
c$$$      end function isSYMBdry

      end module bc_ops

c module grid_operations_cell
c #####################################################################
      module grid_operations_cell

c ---------------------------------------------------------------------
c     This module packs routines that perform operations on grid
c     quantities, such as coordinate transformation of vector components,
c     vector norms and scalar products. It contains the following
c     routines:
c          * transformVectorToCartesian 
c          * transformVectorToCurvilinear 
c          * transformFromCurvToCurv 
c          * vectorNorm
c          * scalarProduct
c     It is assumed that the grid metric structure gmetric is allocated 
c     and filled.
c ---------------------------------------------------------------------

        use grid_packing

        use bc_ops

        implicit none

        INTERFACE crossProduct
          module procedure crossProduct_ijk,crossProduct_mesh
        end INTERFACE

        INTERFACE vectorNorm
          module procedure vectorNorm_ijk,vectorNorm_mesh
        end INTERFACE

        INTERFACE scalarProduct
          module procedure scalarProduct_ijk,scalarProduct_mesh
        end INTERFACE

        INTERFACE XformToCar
          module procedure XformToCar_ijk,XformToCar_mesh
        end INTERFACE

        INTERFACE XformToCurv
          module procedure XformToCurv_ijk,XformToCurv_mesh
        end INTERFACE

        INTERFACE XformToCnv
          module procedure XformToCnv_ijk,XformToCnv_mesh
        end INTERFACE

        INTERFACE XformToCov
          module procedure XformToCov_ijk,XformToCov_mesh
        end INTERFACE

        INTERFACE XformToCov_vrtx
          module procedure XformToCov_vrtx_ijk,XformToCov_vrtx_mesh
        end INTERFACE

        INTERFACE XformVector
          module procedure XformVector_mesh!,XformVector_mesh2
        end INTERFACE

cc        private :: transformFromCurvToCurv
cc     .             transformVectorToCartesian
cc     .            ,transformVectorToCurvilinear

        logical,private :: no_map=.false.

      contains

c     have_map
c     #################################################################
      function have_map()

      logical :: have_map

      have_map = .not.no_map
      
      end function have_map

c     set_map
c     #################################################################
      subroutine set_map(is_map)

      logical :: is_map

      no_map = .not.is_map
      
      end subroutine set_map

c     getMGmap
c     #################################################################
      subroutine getMGmap(g_def,i,j,k,igx,igy,igz,ig,jg,kg)

c     -----------------------------------------------------------------
c     Gets MG vector components (ig,jg,kg) for grid quantities
c     corresponding to node position (i,j,k) in grid levels igx,igy,igz
c     -----------------------------------------------------------------

        implicit none

c     Input variables

        integer :: i,j,k,igx,igy,igz,ig,jg,kg

        type(grid_mg_def),pointer :: g_def

c     Local variables

c     Begin program

        ig = i + g_def%istartx(igx)
        jg = j + g_def%istarty(igy)
        kg = k + g_def%istartz(igz)

      end subroutine getMGmap

c     getCurvilinearCoordinates
c     #################################################################
      subroutine getCurvilinearCoordinates(g_def,i,j,k,igx,igy,igz
     .                                    ,ig,jg,kg,x1,y1,z1,vertex)

c     -----------------------------------------------------------------
c     Finds curvilinear coordinates for position (i,j,k)
c     -----------------------------------------------------------------

        implicit none

c     Input variables

        type(grid_mg_def),pointer :: g_def

        integer :: i,j,k,igx,igy,igz,ig,jg,kg
        real(8) :: x1,y1,z1
        logical,optional :: vertex

c     Local variables

        logical :: vrtx

c     Begin program

        call getMGmap(g_def,i,j,k,igx,igy,igz,ig,jg,kg)

        if (PRESENT(vertex)) then
          vrtx = vertex
        else
          vrtx = .false.
        endif

        if (vrtx) then
          x1 = 0.5*(g_def%xx(ig)+g_def%xx(ig+1))
          y1 = 0.5*(g_def%yy(jg)+g_def%yy(jg+1))
          z1 = 0.5*(g_def%zz(kg)+g_def%zz(kg+1))
        else
          x1 = g_def%xx(ig)
          y1 = g_def%yy(jg)
          z1 = g_def%zz(kg)
        endif

      end subroutine getCurvilinearCoordinates

c     getCartesianCoordinates
c     #################################################################
      subroutine getCartesianCoordinates(g_def,i,j,k,igx,igy,igz
     .                                  ,ig,jg,kg,x1,y1,z1)

c     -----------------------------------------------------------------
c     Returns Cartesian coordinates.
c     -----------------------------------------------------------------

        implicit none

c     Input variables

        integer :: i,j,k,igx,igy,igz,ig,jg,kg
        real(8) :: x1,y1,z1

        type(grid_mg_def),pointer :: g_def

c     Local variables

c     Begin program

        !Find coordinates
        x1 = g_def%gmetric%grid(igx)%car(i,j,k,1)
        y1 = g_def%gmetric%grid(igy)%car(i,j,k,2)
        z1 = g_def%gmetric%grid(igz)%car(i,j,k,3)

      end subroutine getCartesianCoordinates

c$$$c     vecPos
c$$$c     #################################################################
c$$$      function vecPos(neq,i,j,k,igx,igy,igz) result (ijk)
c$$$
c$$$        implicit none
c$$$
c$$$        integer    :: ijk,i,j,k,neq,igx,igy,igz
c$$$
c$$$        integer    :: nxl,nyl,nzl,ilo,jlo,klo
c$$$
c$$$        ilo = 1
c$$$        jlo = 1
c$$$        klo = 1
c$$$
c$$$        nxl = grid_params%nxv(igx)
c$$$        nyl = grid_params%nyv(igy)
c$$$        nzl = grid_params%nzv(igz)
c$$$
c$$$        ijk = neq*(i-ilo + nxl*(j-jlo) + nxl*nyl*(k-klo))
c$$$      
c$$$      end function vecPos

c     vPos
c     #################################################################
      function vPos(g_def,igr,ilv) result (ijk)

        implicit none

        integer :: igr,ilv(3),ijk

        integer :: nxl,nyl,nzl,ilo,jlo,klo

        type(grid_mg_def),pointer :: g_def

        ilo = 1
        jlo = 1
        klo = 1

        nxl = g_def%nxv(igr)
        nyl = g_def%nyv(igr)
        nzl = g_def%nzv(igr)

        ijk = ilv(1)-ilo + nxl*(ilv(2)-jlo) + nxl*nyl*(ilv(3)-klo)
      
      end function vPos

c     fromGlobalToLocalLimits
c     #################################################################
      subroutine fromGlobalToLocalLimits(g_def,igr,ig,jg,kg,il,jl,kl)

        implicit none

        type(grid_mg_def),pointer :: g_def

        integer :: ig,jg,kg,il,jl,kl,igr

        il = ig - g_def%ilo(igr) + 1
        jl = jg - g_def%jlo(igr) + 1
        kl = kg - g_def%klo(igr) + 1

      end subroutine fromGlobalToLocalLimits

c     fromLocalToGlobalLimits
c     #################################################################
      subroutine fromLocalToGlobalLimits(g_def,igr,il,jl,kl,ig,jg,kg)

        implicit none

        type(grid_mg_def),pointer :: g_def

        integer :: ig,jg,kg,il,jl,kl,igr

        ig = il + g_def%ilo(igr) - 1
        jg = jl + g_def%jlo(igr) - 1
        kg = kl + g_def%klo(igr) - 1

      end subroutine fromLocalToGlobalLimits

c     localIndex
c     #################################################################
      function localIndex(g_def,igrid,gindex) result(lindex)

        implicit none

        integer  :: igrid,lindex(3),gindex(3)

        type(grid_mg_def),pointer :: g_def

        lindex(1) = gindex(1) - g_def%ilo(igrid) + 1
        lindex(2) = gindex(2) - g_def%jlo(igrid) + 1
        lindex(3) = gindex(3) - g_def%klo(igrid) + 1

      end function localIndex

c     globalIndex
c     #################################################################
      function globalIndex(g_def,igrid,lindex) result(gindex)

        implicit none

        integer  :: igrid,lindex(3),gindex(3)

        type(grid_mg_def),pointer :: g_def

        gindex(1) = lindex(1) + g_def%ilo(igrid) - 1
        gindex(2) = lindex(2) + g_def%jlo(igrid) - 1
        gindex(3) = lindex(3) + g_def%klo(igrid) - 1

      end function globalIndex

c     XformToCar_ijk
c     #################################################################
      function XformToCar_ijk(g_def,i,j,k,igrid,curv,covariant)
     .         result(car)

c     -----------------------------------------------------------------
c     Transforms a curvilinear vector (c1,c2,c3) to Cartesian (cx,cy,cz)
c     at grid coordinates (i,j,k). Curvilinear vector is covariant if
c     covariant=.true., and contravariant otherwise.
c     -----------------------------------------------------------------

        implicit none

c     Input variables

        type(grid_mg_def),pointer :: g_def

        integer    :: i,j,k,igrid
        real(8)    :: curv(3),car(3)
        logical    :: covariant

c     Local variables

        type(MG_grid_metrics),pointer :: gmetric 

c     Begin program

        gmetric => g_def%gmetric

        if (covariant) then
          car = matmul(curv,gmetric%grid(igrid)%cov(i,j,k,:,:))
        else
          car = matmul(curv,gmetric%grid(igrid)%cnv(i,j,k,:,:))
        endif

      end function XformToCar_ijk

c     XformToCar_vrtx_ijk
c     #################################################################
      function XformToCar_vrtx_ijk(g_def,i,j,k,igrid,curv,covariant)
     .         result(car)

c     -----------------------------------------------------------------
c     Transforms a curvilinear vector (c1,c2,c3) to Cartesian (cx,cy,cz)
c     at grid coordinates (i,j,k). Curvilinear vector is covariant if
c     covariant=.true., and contravariant otherwise.
c     -----------------------------------------------------------------

        implicit none

c     Input variables

        type(grid_mg_def),pointer :: g_def

        integer    :: i,j,k,igrid
        real(8)    :: curv(3),car(3)
        logical    :: covariant

c     Local variables

        type(MG_grid_metrics),pointer :: gmetric 

c     Begin program

        gmetric => g_def%gmetric

        if (covariant) then
          car = matmul(curv,gmetric%grid(igrid)%cov_v(i,j,k,:,:))
        else
          car = matmul(curv,gmetric%grid(igrid)%cnv_v(i,j,k,:,:))
        endif

      end function XformToCar_vrtx_ijk

c     XformToCar_mesh
c     #################################################################
      function XformToCar_mesh(g_def,igrid,curv,covariant) result(car)

c     -----------------------------------------------------------------
c     Transforms a curvilinear vector to Cartesian 
c     over the whole mesh. Curvilinear vector is covariant if
c     covariant=.true., and contravariant otherwise.
c     -----------------------------------------------------------------

        implicit none

c     Input variables

        type(grid_mg_def),pointer :: g_def

        integer :: i,j,k,igrid
        real(8) :: curv(0:,0:,0:,:)
     .            ,car(0:size(curv,1)-1
     .                ,0:size(curv,2)-1
     .                ,0:size(curv,3)-1,size(curv,4))
        logical :: covariant

c     Local variables

        integer :: nxl,nyl,nzl

c     Begin program

        nxl = g_def%nxv(igrid)
        nyl = g_def%nyv(igrid)
        nzl = g_def%nzv(igrid)

        do k=0,nzl+1
          do j=0,nyl+1
            do i=0,nxl+1
              car(i,j,k,:) = XformToCar_ijk(g_def,i,j,k,igrid
     .                                     ,curv(i,j,k,:)
     .                                     ,covariant)
            enddo
          enddo
        enddo

      end function XformToCar_mesh

c     XformToCurv_ijk
c     #################################################################
      function XformToCurv_ijk(g_def,i,j,k,igrid,car,tocnv,nojac)
     .         result(curv)

c     -----------------------------------------------------------------
c     Transforms a Cartesian vector to curvilinear representation at
c     grid coordinates (i,j,k). Curvilinear vector is covariant if
c     tocnv=.false., and contravariant otherwise.
c     -----------------------------------------------------------------

        implicit none

c     Input variables

        type(grid_mg_def),pointer :: g_def

        integer :: i,j,k,igrid
        real(8) :: car(3),curv(3)
        logical :: tocnv

        logical,optional :: nojac

c     Local variables

        real(8) :: jac
        logical :: nj
        
        type(MG_grid_metrics),pointer :: gmetric 

c     Begin program

        gmetric => g_def%gmetric

        if (PRESENT(nojac)) then
          nj = nojac
        else
          nj = .false.
        endif

        if (nj) then
          jac = 1d0
        else
          jac = gmetric%grid(igrid)%jac(i,j,k)
        endif

        if (.not.tocnv) then
          curv = jac*matmul(gmetric%grid(igrid)%cnv(i,j,k,:,:),car)
        else
          curv = jac*matmul(gmetric%grid(igrid)%cov(i,j,k,:,:),car)
        endif

      end function XformToCurv_ijk

c     XformToCurv_vrtx_ijk
c     #################################################################
      function XformToCurv_vrtx_ijk(g_def,i,j,k,igrid,car,tocnv,nojac)
     .         result(curv)

c     -----------------------------------------------------------------
c     Transforms a Cartesian vector to curvilinear representation at
c     vertex grid coordinates (i+1/2,j+1/2,k+1/2). Returned curvilinear
c     vector is covariant if tocnv=.false., and contravariant otherwise.
c     -----------------------------------------------------------------

        implicit none

c     Input variables

        type(grid_mg_def),pointer :: g_def

        integer :: i,j,k,igrid
        real(8) :: car(3),curv(3)
        logical :: tocnv

        logical,optional :: nojac

c     Local variables

        real(8) :: jac
        logical :: nj
        
        type(MG_grid_metrics),pointer :: gmetric 

c     Begin program

        gmetric => g_def%gmetric

        if (PRESENT(nojac)) then
          nj = nojac
        else
          nj = .false.
        endif

        if (nj) then
          jac = 1d0
        else
          jac = gmetric%grid(igrid)%jac_v(i,j,k)
          if (jac == 0d0) jac = SP_flsv
        endif

        if (.not.tocnv) then
          curv = jac*matmul(gmetric%grid(igrid)%cnv_v(i,j,k,:,:),car)
        else
          curv = jac*matmul(gmetric%grid(igrid)%cov_v(i,j,k,:,:),car)
        endif

c$$$        if (nj) then
c$$$          jac = 1d0
c$$$        else
c$$$          if (gmetric%grid(igrid)%jac_v(i,j,k) == 0d0.and.tocnv) then
c$$$            jac(1) = 0d0
c$$$            jac(2) = SP_flsv
c$$$            jac(3) = 0d0
c$$$          elseif (gmetric%grid(igrid)%jac_v(i,j,k) == 0d0) then
c$$$            jac(1) = SP_flsv
c$$$            jac(2) = 0d0
c$$$            jac(3) = SP_flsv
c$$$          else
c$$$            jac = gmetric%grid(igrid)%jac_v(i,j,k)            
c$$$          endif
c$$$        endif
c$$$
c$$$        if (.not.tocnv) then
c$$$          vec(1,:) = jac(1)*gmetric%grid(igrid)%cnv_v(i,j,k,1,:)
c$$$          vec(2,:) = jac(2)*gmetric%grid(igrid)%cnv_v(i,j,k,2,:)
c$$$          vec(3,:) = jac(3)*gmetric%grid(igrid)%cnv_v(i,j,k,3,:)
c$$$c$$$          if (gmetric%grid(igrid)%jac_v(i,j,k) == 0d0)
c$$$c$$$     .         write (*,*) vec(1,:)
c$$$        else
c$$$          vec(1,:) = jac(1)*gmetric%grid(igrid)%cov_v(i,j,k,1,:)
c$$$          vec(2,:) = jac(2)*gmetric%grid(igrid)%cov_v(i,j,k,2,:)
c$$$          vec(3,:) = jac(3)*gmetric%grid(igrid)%cov_v(i,j,k,3,:)
c$$$        endif
c$$$
c$$$        curv = matmul(vec,car)
        
      end function XformToCurv_vrtx_ijk

c     XformToCurv_mesh
c     #################################################################
      function XformToCurv_mesh(g_def,igrid,car,tocnv) result(curv)

c     -----------------------------------------------------------------
c     Transforms a Cartesian vector to curvilinear
c     over the whole mesh. Curvilinear vector is contravariant if
c     tocnv=.true., and covariant otherwise.
c     -----------------------------------------------------------------

        implicit none

c     Input variables

        type(grid_mg_def),pointer :: g_def

        integer :: i,j,k,igrid
        real(8) :: car(0:,0:,0:,:)
     .            ,curv(0:size(car,1)-1
     .                 ,0:size(car,2)-1
     .                 ,0:size(car,3)-1,size(car,4))
        logical :: tocnv

c     Local variables

        integer :: nxl,nyl,nzl

c     Begin program

        nxl = g_def%nxv(igrid)
        nyl = g_def%nyv(igrid)
        nzl = g_def%nzv(igrid)

        do k=0,nzl+1
          do j=0,nyl+1
            do i=0,nxl+1
              curv(i,j,k,:) = XformToCurv_ijk(g_def,i,j,k,igrid
     .                                       ,car(i,j,k,:)
     .                                       ,tocnv)
            enddo
          enddo
        enddo

      end function XformToCurv_mesh

c     XformToCurv_mesh2
c     #################################################################
      subroutine XformToCurv_mesh2(g_def,igrid,vec,tocnv)

c     -----------------------------------------------------------------
c     Transforms a Cartesian vector vec to curvilinear
c     over the whole mesh. Curvilinear vector is contravariant if
c     tocnv=.true., and covariant otherwise. Result is returned in vec.
c     -----------------------------------------------------------------

        implicit none

c     Input variables

        type(grid_mg_def),pointer :: g_def

        integer :: i,j,k,igrid
        real(8) :: vec(0:,0:,0:,:)
        logical :: tocnv

c     Local variables

        integer :: nxl,nyl,nzl
        real(8) :: dum(size(vec,4))

c     Begin program

        nxl = g_def%nxv(igrid)
        nyl = g_def%nyv(igrid)
        nzl = g_def%nzv(igrid)

        do k=0,nzl+1
          do j=0,nyl+1
            do i=0,nxl+1
              dum = vec(i,j,k,:)
              vec(i,j,k,:) =XformToCurv_ijk(g_def,i,j,k,igrid,dum,tocnv)
            enddo
          enddo
        enddo

      end subroutine XformToCurv_mesh2

c     transformVectorToCurvilinear
c     #################################################################
      subroutine transformVectorToCurvilinear(g_def,i,j,k,igx,igy,igz
     .                                       ,cx,cy,cz,tocnv
     .                                       ,c1,c2,c3)

c     -----------------------------------------------------------------
c     Transforms a Cartesian vector (cx,cy,cz) to curvilinear (c1,c2,c3)
c     (contravariant if tocnv=.true., contravariant otherwise)
c     at grid coordinates (i,j,k). 
c     -----------------------------------------------------------------

        implicit none

c     Input variables

        type(grid_mg_def),pointer :: g_def

        integer    :: i,j,k,igx,igy,igz
        real(8)    :: c1,c2,c3,cx,cy,cz
        logical    :: tocnv

c     Local variables

        integer    :: ic
        real(8)    :: T_to_curv(3,3),vec(3),jac

        type(MG_grid_metrics),pointer :: gmetric 

c     Begin program

        gmetric => g_def%gmetric

        jac = gmetric%grid(igx)%jac(i,j,k)

        vec = (/ cx,cy,cz /)

        if (.not.tocnv) then

          c1 = jac*dot_product(vec,gmetric%grid(igx)%cnv(i,j,k,1,:))
          c2 = jac*dot_product(vec,gmetric%grid(igx)%cnv(i,j,k,2,:))
          c3 = jac*dot_product(vec,gmetric%grid(igx)%cnv(i,j,k,3,:))
            
        else

          c1 = jac*dot_product(vec,gmetric%grid(igx)%cov(i,j,k,1,:))
          c2 = jac*dot_product(vec,gmetric%grid(igx)%cov(i,j,k,2,:))
          c3 = jac*dot_product(vec,gmetric%grid(igx)%cov(i,j,k,3,:))

        endif

      end subroutine transformVectorToCurvilinear

c     transformFromCurvToCurv
c     #################################################################
      subroutine transformFromCurvToCurv(g_def,i,j,k,igx,igy,igz
     .             ,cov1,cov2,cov3,cnv1,cnv2,cnv3,tocnv,half_elem
     .             ,hex,hey,hez)
c     -----------------------------------------------------------------
c     Transforms a curvilinear vector from covariant to contravariant 
c     (tocnv=.true.) or viceversa at cell centers, faces, and vertices
c     of grid coordinate (i,j,k).
c     -----------------------------------------------------------------

        implicit none

c     Input variables

        type(grid_mg_def),pointer :: g_def

        integer :: i,j,k,igx,igy,igz
        integer,optional,intent(IN) :: half_elem,hex,hey,hez

        real(8) :: cov1,cov2,cov3,cnv1,cnv2,cnv3
        logical :: tocnv

c     Local variables

        integer :: he,he_x,he_y,he_z,ip,jp,kp
        real(8) :: cov(3),cnv(3),tensor(3,3),coeff
        real(8),pointer,dimension(:,:,:,:,:) :: gtnsr

        type(MG_grid_metrics),pointer :: gmetric 

c     Begin program

        gmetric => g_def%gmetric

c     Special case of cartesian coords

        if (no_map) then

          if (tocnv) then
            cnv1 = cov1
            cnv2 = cov2
            cnv3 = cov3
          else
            cov1 = cnv1
            cov2 = cnv2
            cov3 = cnv3
          endif

          return
        endif

c     Optional parameters

        if (PRESENT(half_elem)) then
          he = half_elem
        else
          he = 0
        endif

        if (PRESENT(hex)) then
          he_x = hex
        else
          he_x = 0
        endif

        if (PRESENT(hey)) then
          he_y = hey
        else
          he_y = 0
        endif

        if (PRESENT(hez)) then
          he_z = hez
        else
          he_z = 0
        endif

c     Consistency check

        if (he == 1) then
          he_x = 1; he_y = 0 ; he_z = 0
        endif

        if (he == 2) then
          he_x = 0; he_y = 1 ; he_z = 0
        endif

        if (he == 3) then
          he_x = 0; he_y = 0 ; he_z = 1
        endif

c     Find metric tensor

        ip = i+1
        jp = j+1
        kp = k+1

        if (tocnv) then
          gtnsr => gmetric%grid(igx)%gsup
        else
          gtnsr => gmetric%grid(igx)%gsub
        endif

        tensor = gtnsr(i,j,k,:,:)
        coeff = 1d0

        if (he_x == 1) then

          tensor = tensor + gtnsr(ip,j,k,:,:)
          coeff = 0.5

          if (he_y == 1) then
            tensor = tensor + gtnsr(i ,jp,k,:,:)
     .                      + gtnsr(ip,jp,k,:,:)
            coeff = 0.25
          endif

          if (he_z == 1) then
            tensor = tensor + gtnsr(i ,j ,kp,:,:)
     .                      + gtnsr(ip,j ,kp,:,:)
     .                      + gtnsr(i ,jp,kp,:,:)
     .                      + gtnsr(ip,jp,kp,:,:)
            coeff = 0.125
          endif

        elseif (he_y == 1) then

          tensor = tensor + gtnsr(i,jp,k,:,:)
          coeff = 0.5

          if (he_z == 1) then
            tensor = tensor + gtnsr(i,j ,kp,:,:)
     .                      + gtnsr(i,jp,kp,:,:)
            coeff = 0.25
          endif

        elseif (he_z == 1) then

          tensor = tensor + gtnsr(i,j,kp,:,:)
          coeff = 0.5

        endif

        tensor = coeff*tensor

c     Transform vector

        if (tocnv) then
          cnv1 = tensor(1,1)*cov1+tensor(1,2)*cov2+tensor(1,3)*cov3
          cnv2 = tensor(2,1)*cov1+tensor(2,2)*cov2+tensor(2,3)*cov3
          cnv3 = tensor(3,1)*cov1+tensor(3,2)*cov2+tensor(3,3)*cov3
        else
          cov1 = tensor(1,1)*cnv1+tensor(1,2)*cnv2+tensor(1,3)*cnv3
          cov2 = tensor(2,1)*cnv1+tensor(2,2)*cnv2+tensor(2,3)*cnv3
          cov3 = tensor(3,1)*cnv1+tensor(3,2)*cnv2+tensor(3,3)*cnv3
        endif

      end subroutine transformFromCurvToCurv

c     transformVector
c     ######################################################################
      subroutine transformVector(g_def,igx,igy,igz
     .                          ,imin,imax,jmin,jmax,kmin,kmax
     .                          ,arr1,arr2,arr3,in_repr,out_repr)

c     ----------------------------------------------------------------------
c     Transforms vectors components in arrays arr1,arr2,arr3 from
c     Cartesian to curvilinear and viceversa, or from curvilinear to
c     curvilinear, depending on the value of in_repr (which indicates the
c     representation of the input vector) and out_repr (which indicates
c     representation of the input vector). Possible representations are:
c     'car','cnv',and 'cov'.
c     ----------------------------------------------------------------------

        implicit none

c     Input variables

        type(grid_mg_def),pointer :: g_def

        integer    :: imin,imax,jmin,jmax,kmin,kmax
        integer    :: igx,igy,igz
        real(8)    :: arr1(imin:imax,jmin:jmax,kmin:kmax)
     .               ,arr2(imin:imax,jmin:jmax,kmin:kmax)
     .               ,arr3(imin:imax,jmin:jmax,kmin:kmax)

        character(3) :: in_repr,out_repr

c     Local variables

        integer    :: i,j,k
        real(8)    :: vec(3)

        logical    :: covariant,cnvariant,cartesian

c     Begin program

        cartesian = (in_repr == 'car')
        covariant = (in_repr == 'cov')
        cnvariant = (.not. cartesian) .and. (.not. covariant)

        select case(out_repr)
        case('car')

          if (cartesian) return

          do k=kmin,kmax
            do j=jmin,jmax
              do i=imin,imax

                vec = XformToCar(g_def,i,j,k,igx
     .               ,(/arr1(i,j,k),arr2(i,j,k),arr3(i,j,k)/)
     .               ,covariant)

                arr1(i,j,k) = vec(1)
                arr2(i,j,k) = vec(2)
                arr3(i,j,k) = vec(3)
                
              enddo
            enddo
          enddo

        case('cnv')

          if (cnvariant) return

          if (cartesian) then
            do k=kmin,kmax
              do j=jmin,jmax
                do i=imin,imax

                  call transformVectorToCurvilinear
     .                 (g_def,i,j,k,igx,igy,igz
     .                 ,arr1(i,j,k),arr2(i,j,k),arr3(i,j,k)
     .                 ,.true.
     .                 ,vec(1),vec(2),vec(3))

                  arr1(i,j,k) = vec(1)
                  arr2(i,j,k) = vec(2)
                  arr3(i,j,k) = vec(3)
                  
                enddo
              enddo
            enddo
          else
            do k=kmin,kmax
              do j=jmin,jmax
                do i=imin,imax

                  call transformFromCurvToCurv
     .                 (g_def,i,j,k,igx,igy,igz
     .                 ,arr1(i,j,k),arr2(i,j,k),arr3(i,j,k)
     .                 ,vec(1),vec(2),vec(3)
     .                 ,.true.,0)

                  arr1(i,j,k) = vec(1)
                  arr2(i,j,k) = vec(2)
                  arr3(i,j,k) = vec(3)
                  
                enddo
              enddo
            enddo
          endif

        case('cov')

          if (covariant) return

          if (cartesian) then
            do k=kmin,kmax
              do j=jmin,jmax
                do i=imin,imax

                  call transformVectorToCurvilinear
     .                 (g_def,i,j,k,igx,igy,igz
     .                 ,arr1(i,j,k),arr2(i,j,k),arr3(i,j,k)
     .                 ,.false.
     .                 ,vec(1),vec(2),vec(3))

                  arr1(i,j,k) = vec(1)
                  arr2(i,j,k) = vec(2)
                  arr3(i,j,k) = vec(3)
                  
                enddo
              enddo
            enddo
          else
            do k=kmin,kmax
              do j=jmin,jmax
                do i=imin,imax

                  call transformFromCurvToCurv
     .                 (g_def,i,j,k,igx,igy,igz
     .                 ,vec(1),vec(2),vec(3)
     .                 ,arr1(i,j,k),arr2(i,j,k),arr3(i,j,k)
     .                 ,.false.,0)

                  arr1(i,j,k) = vec(1)
                  arr2(i,j,k) = vec(2)
                  arr3(i,j,k) = vec(3)
                  
                enddo
              enddo
            enddo
          endif

        end select

      end subroutine transformVector

c     XformVector_mesh2
c     #################################################################
      subroutine XformVector_mesh2(g_def,igrid,vec,in_repr,out_repr)

c     -----------------------------------------------------------------
c     Transforms vector vec from Cartesian to curvilinear and viceversa,
c     or from curvilinear to curvilinear, depending on the value of
c     in_repr (which indicates the representation of the input vector)
c     and out_repr (which indicates representation of the input
c     vector). Possible representations are: 'car','cnv',and 'cov'.
c     -----------------------------------------------------------------

        implicit none

c     Input variables

        type(grid_mg_def),pointer :: g_def

        integer :: igrid
        real(8) :: vec(:,:,:,:)

        character(3) :: in_repr,out_repr

c     Local variables

        integer :: imin,imax,jmin,jmax,kmin,kmax

c     Begin program

        imin = 0
        imax = g_def%nxv(igrid)+1

        jmin = 0
        jmax = g_def%nyv(igrid)+1

        kmin = 0
        kmax = g_def%nzv(igrid)+1

        call  transformVector(g_def,igrid,igrid,igrid
     .                       ,imin,imax,jmin,jmax,kmin,kmax
     .                       ,vec(:,:,:,1)
     .                       ,vec(:,:,:,2)
     .                       ,vec(:,:,:,3),in_repr,out_repr)

c     End program

      end subroutine XformVector_mesh2

c     XformVector_mesh
c     #################################################################
      function XformVector_mesh(g_def,igrid,vec1,tocnv) result(vec2)
c     -----------------------------------------------------------------
c     Transforms a curvilinear vector from covariant to contravariant 
c     (tocnv=.true.) or viceversa at cell centers on all mesh.
c     -----------------------------------------------------------------

      implicit none

c     Input variables

      type(grid_mg_def),pointer :: g_def

      integer    :: igrid
      real(8)    :: vec1(0:,0:,0:,:)
      real(8)    :: vec2(0:size(vec1,1)-1
     .                  ,0:size(vec1,2)-1
     .                  ,0:size(vec1,3)-1
     .                  ,size(vec1,4))

      logical    :: tocnv

c     Local variables

      integer    :: i,j,k,nxl,nyl,nzl
      real(8),pointer,dimension(:,:,:,:,:) :: tensor

      type(MG_grid_metrics),pointer :: gmetric 

c     Begin program

      gmetric => g_def%gmetric

      nxl = g_def%nxv(igrid)
      nyl = g_def%nyv(igrid)
      nzl = g_def%nzv(igrid)

      if (tocnv) then
        tensor => gmetric%grid(igrid)%gsup
      else
        tensor => gmetric%grid(igrid)%gsub
      endif

      do k=0,nzl+1
        do j=0,nyl+1
          do i=0,nxl+1
            vec2(i,j,k,:) = matmul(vec1(i,j,k,:),tensor(i,j,k,:,:))
          enddo
        enddo
      enddo

      nullify(tensor)

c     End program

      end function XformVector_mesh

c     XformVector_BC
c     ###################################################################
      subroutine XformVector_BC(g_def,igrid,vec1,vec2,tocnv)

      implicit none

c     -------------------------------------------------------------------
c     Routine to transfer ghost cell values form vec1 to vec2
c     -------------------------------------------------------------------

c     Input variables

      type(grid_mg_def),pointer :: g_def

      integer :: igrid
      real(8) :: vec1(0:,0:,0:,:)
      real(8) :: vec2(0:,0:,0:,:)
      logical :: tocnv

c     Local variables

      integer :: nxp,nyp,nzp,i,j,k
      real(8),pointer,dimension(:,:,:,:,:) :: tensor

      type(MG_grid_metrics),pointer :: gmetric 

c     Begin program

      gmetric => g_def%gmetric

      nxp = size(vec1,1)-1
      nyp = size(vec1,2)-1
      nzp = size(vec1,3)-1

      if (tocnv) then
        tensor => gmetric%grid(igrid)%gsup
      else
        tensor => gmetric%grid(igrid)%gsub
      endif

      !X boundaries
      do k = 0,nzp
        do j= 0,nyp
          vec2(0  ,j,k,:) = matmul(tensor(0  ,j,k,:,:),vec1(0  ,j,k,:))
          vec2(nxp,j,k,:) = matmul(tensor(nxp,j,k,:,:),vec1(nxp,j,k,:))
        enddo
      enddo

      !Y boundaries
      do k= 0,nzp
        do i=0,nxp
          vec2(i,0  ,k,:) = matmul(tensor(i,0  ,k,:,:),vec1(i,0  ,k,:))
          vec2(i,nyp,k,:) = matmul(tensor(i,nyp,k,:,:),vec1(i,nyp,k,:))
        enddo
      enddo

      !Z boundaries
      do j= 0,nyp
        do i= 0,nxp
          vec2(i,j,0  ,:) = matmul(tensor(i,j,0  ,:,:),vec1(i,j,0  ,:))
          vec2(i,j,nzp,:) = matmul(tensor(i,j,nzp,:,:),vec1(i,j,nzp,:))
        enddo
      enddo

c     End program

      end subroutine XformVector_BC

c     XformToCnv_ijk
c     #################################################################
      function XformToCnv_ijk(g_def,i,j,k,igrid,cov) result(cnv)

c     -----------------------------------------------------------------
c     Transforms a curvilinear vector (c1,c2,c3) to Cartesian (cx,cy,cz)
c     at grid coordinates (i,j,k). Curvilinear vector is covariant if
c     covariant=.true., and contravariant otherwise.
c     -----------------------------------------------------------------

        implicit none

c     Input variables

        type(grid_mg_def),pointer :: g_def

        integer    :: i,j,k,igrid
        real(8)    :: cov(3),cnv(3)

c     Local variables

        type(MG_grid_metrics),pointer :: gmetric 

c     Begin program

        gmetric => g_def%gmetric

        cnv = matmul(cov,gmetric%grid(igrid)%gsup(i,j,k,:,:))

      end function XformToCnv_ijk

c     XformToCov_ijk
c     #################################################################
      function XformToCov_ijk(g_def,i,j,k,igrid,cnv) result(cov)

c     -----------------------------------------------------------------
c     Transforms a curvilinear vector (c1,c2,c3) to Cartesian (cx,cy,cz)
c     at grid coordinates (i,j,k). Curvilinear vector is covariant if
c     covariant=.true., and contravariant otherwise.
c     -----------------------------------------------------------------

        implicit none

c     Input variables

        type(grid_mg_def),pointer :: g_def

        integer    :: i,j,k,igrid
        real(8)    :: cov(3),cnv(3)

c     Local variables

        type(MG_grid_metrics),pointer :: gmetric 

c     Begin program

        gmetric => g_def%gmetric

        cov = matmul(cnv,gmetric%grid(igrid)%gsub(i,j,k,:,:))

      end function XformToCov_ijk

c     XformToCov_vrtx_ijk
c     #################################################################
      function XformToCov_vrtx_ijk(g_def,i,j,k,igrid,cnv)
     .         result(cov)

c     -----------------------------------------------------------------
c     Transforms a contravariant vector into covariant at a vertex.
c     -----------------------------------------------------------------

        implicit none

c     Input variables

        type(grid_mg_def),pointer :: g_def

        integer    :: i,j,k,igrid
        real(8)    :: cov(3),cnv(3)

c     Local variables

        real(8)    :: car(3)

c     Begin program

        cov = matmul(cnv,g_def%gmetric%grid(igrid)%gsub_v(i,j,k,:,:))
        
      end function XformToCov_vrtx_ijk

c     XformToCov_ijk_face
c     #################################################################
      function XformToCov_ijk_face(g_def,i,j,k,igrid,cnv,he)
     .         result(cov)

c     -----------------------------------------------------------------
c     Transforms a curvilinear vector (c1,c2,c3) to Cartesian (cx,cy,cz)
c     at grid coordinates (i,j,k). Curvilinear vector is covariant if
c     covariant=.true., and contravariant otherwise.
c     -----------------------------------------------------------------

        implicit none

c     Input variables

        type(grid_mg_def),pointer :: g_def

        integer    :: i,j,k,igrid,he,hex,hey,hez
        real(8)    :: cov(3),cnv(3)

c     Local variables

c     Begin program

        call transformFromCurvToCurv(g_def,i,j,k,igrid,igrid,igrid
     .                              ,cov(1),cov(2),cov(3)
     .                              ,cnv(1),cnv(2),cnv(3)
     .                              ,.false.,half_elem=he)

      end function XformToCov_ijk_face

c     XformToCnv_ijk_face
c     #################################################################
      function XformToCnv_ijk_face(g_def,i,j,k,igrid,cov,he)
     .         result(cnv)

c     -----------------------------------------------------------------
c     Transforms covariant vector into contravariant at a face.
c     -----------------------------------------------------------------

        implicit none

c     Input variables

        type(grid_mg_def),pointer :: g_def

        integer    :: i,j,k,igrid,he
        real(8)    :: cov(3),cnv(3)

c     Local variables

c     Begin program

        call transformFromCurvToCurv(g_def,i,j,k,igrid,igrid,igrid
     .                              ,cov(1),cov(2),cov(3)
     .                              ,cnv(1),cnv(2),cnv(3)
     .                              ,.true.,half_elem=he)

      end function XformToCnv_ijk_face

c     XformToCnv_mesh
c     #################################################################
      function XformToCnv_mesh(g_def,igrid,cov) result(cnv)

c     -----------------------------------------------------------------
c     Transforms a covariant vector to a contravariant representation
c     everywhere on the mesh.
c     -----------------------------------------------------------------

        implicit none

c     Input variables

        type(grid_mg_def),pointer :: g_def

        integer :: igrid
        real(8),INTENT(IN) :: cov(0:,0:,0:,:)
        real(8) :: cnv(0:size(cov,1)-1
     .                ,0:size(cov,2)-1
     .                ,0:size(cov,3)-1,3)

c     Local variables

        integer :: i,j,k,nx,ny,nz

c     Begin program

        nx = size(cov,1)-2
        ny = size(cov,2)-2
        nz = size(cov,3)-2

        do k=0,nz+1
          do j=0,ny+1
            do i=0,nx+1
              cnv(i,j,k,:) = XformToCnv_ijk(g_def,i,j,k,igrid
     .                                     ,cov(i,j,k,:))
            enddo
          enddo
        enddo

      end function XformToCnv_mesh

c     XformToCov_mesh
c     #################################################################
      function XformToCov_mesh(g_def,igrid,cnv) result(cov)

c     -----------------------------------------------------------------
c     Transforms a contravariant vector to a covariant representation
c     everywhere on the mesh.
c     -----------------------------------------------------------------

        implicit none

c     Input variables

        type(grid_mg_def),pointer :: g_def

        integer :: igrid
        real(8),INTENT(IN) :: cnv(0:,0:,0:,:)
        real(8) :: cov(0:size(cnv,1)-1
     .                ,0:size(cnv,2)-1
     .                ,0:size(cnv,3)-1,3)

c     Local variables

        integer :: i,j,k,nx,ny,nz

c     Begin program

        nx = size(cnv,1)-2
        ny = size(cnv,2)-2
        nz = size(cnv,3)-2

        do k=0,nz+1
          do j=0,ny+1
            do i=0,nx+1
              cov(i,j,k,:) = XformToCov_ijk(g_def,i,j,k,igrid
     .                                     ,cnv(i,j,k,:))
            enddo
          enddo
        enddo

      end function XformToCov_mesh

c     XformToCov_vrtx_mesh
c     #################################################################
      function XformToCov_vrtx_mesh(g_def,igrid,cnv) result(cov)

c     -----------------------------------------------------------------
c     Transforms a contravariant vector to a covariant representation
c     everywhere on the mesh.
c     -----------------------------------------------------------------

        use oned_int

        use xdraw_io
     
        implicit none

c     Input variables

        type(grid_mg_def),pointer :: g_def

        integer :: igrid
        real(8),INTENT(IN) :: cnv(0:,0:,0:,:)
        real(8) :: cov(0:size(cnv,1)-1
     .                ,0:size(cnv,2)-1
     .                ,0:size(cnv,3)-1,3)

c     Local variables

        integer :: i,j,k,nx,ny,nz,ig,jg,kg,oint
        real(8) :: xx(10),car(3,3),car0(3),cov0(0:size(cnv,2)-1,3)

        integer :: jlog,jglobal,jj,mm

c     Begin program

        nx = size(cnv,1)-2
        ny = size(cnv,2)-2
        nz = size(cnv,3)-2

        do k=0,nz
          do j=0,ny
            do i=0,nx
              cov(i,j,k,:) = XformToCov_vrtx_ijk(g_def,i,j,k,igrid
     .                                          ,cnv(i,j,k,:))
            enddo
          enddo
        enddo
        
c$$$        if (isSP2(g_def,1,igrid)) then
c$$$          mm = 1
c$$$          if (coords=='hel') mm = nint(g_def%params(1))
c$$$          jlog = g_def%jlo(igrid)
c$$$
c$$$          !Order of interpolation
c$$$          oint = 3
c$$$          !Vertex positions
c$$$          call getMGmap(g_def,1,1,1,igrid,igrid,igrid,ig,jg,kg)
c$$$          xx(1:oint+1) = 0.5*(g_def%xx(ig  :ig+oint+1)
c$$$     .                       +g_def%xx(ig+1:ig+oint+2))
c$$$
c$$$          do k=0,nz
c$$$            !Compute extrapolation at r=0
c$$$            do j=0,ny
c$$$              call IntDriver1d(oint+1,xx(1:oint+1),cov(1:oint+1,j,k,1),1
c$$$     .                        ,(/0d0/),cov0(j:j,1),oint)
c$$$              call IntDriver1d(oint+1,xx(1:oint+1),cov(1:oint+1,j,k,2),1
c$$$     .                        ,(/0d0/),cov0(j:j,2),oint)
c$$$              call IntDriver1d(oint+1,xx(1:oint+1),cov(1:oint+1,j,k,3),1
c$$$     .                        ,(/0d0/),cov0(j:j,3),oint)
c$$$            enddo
c$$$c$$$            !Assign value
c$$$c$$$            cov(0,:,k,:) = cov0(:,:)
c$$$            !Average extrapolation at r=0 across rays
c$$$            do j=0,ny
c$$$              jglobal = j + jlog - 1
c$$$              jj = mod(jglobal+mm*ny/2,ny)
c$$$              if (jj == 0) jj = ny
c$$$              cov(0,j,k,1) = 0.5*(cov0(j,1)-cov0(jj,1)) !Average r+ and r- to origin
c$$$              cov(0,j,k,2) = 0.5*(cov0(j,2)+cov0(jj,2)) !Average r+ and r- to origin
c$$$              cov(0,j,k,3) = 0.5*(cov0(j,3)+cov0(jj,3)) !Average r+ and r- to origin
c$$$            enddo
c$$$          enddo
c$$$        endif

c diag ****
c$$$        write (*,*)
c$$$c$$$        write (*,*) "Cartesian comps at r=0",car0
c$$$        do i=0,nx
c$$$          write (*,*) "Bcnv_vrtx"
c$$$     .         ,0.5*(g_def%xg(i)+g_def%xg(i+1)),cnv(i,jplot,1,:)
c$$$        enddo
c$$$        write (*,*)
c$$$        do i=0,nx
c$$$          write (*,*) "Bcov_vrtx"
c$$$     .         ,0.5*(g_def%xg(i)+g_def%xg(i+1)),cov(i,jplot,1,:)
c$$$        enddo
c$$$        
c$$$        call createDrawInCfile(3,'Bcov_vrtx.bin','Solution','t','x','y'
c$$$     .        ,(/'B_1','B_2','B_3'/),'-c -X0 -L57','drawBcov_vrtx.in')
c$$$
c$$$        open(unit=110,file='Bcov_vrtx.bin',form='unformatted'
c$$$     .      ,status='replace')
c$$$
c$$$        call contour(cov(0:nx,0:ny,1,1),nx+1,ny+1
c$$$     .            ,g_def%gxmin,g_def%gxmax
c$$$     .            ,g_def%gymin,g_def%gymax
c$$$     .            ,0,110)
c$$$        call contour(cov(0:nx,0:ny,1,2),nx+1,ny+1
c$$$     .            ,g_def%gxmin,g_def%gxmax
c$$$     .            ,g_def%gymin,g_def%gymax
c$$$     .            ,1,110)
c$$$        call contour(cov(0:nx,0:ny,1,3),nx+1,ny+1
c$$$     .            ,g_def%gxmin,g_def%gxmax
c$$$     .            ,g_def%gymin,g_def%gymax
c$$$     .            ,1,110)
c$$$        close(110)
c diag ****

      end function XformToCov_vrtx_mesh

c     vectorNorm_ijk
c     ################################################################
      function vectorNorm_ijk(g_def,i,j,k,igrid,vec,covar) result(norm)

c     ---------------------------------------------------------------
c     Finds L2-norm of vector A given its curvilinear components.
c     ---------------------------------------------------------------

      implicit none

c     Call variables

      type(grid_mg_def),pointer :: g_def

      integer    :: i,j,k,igrid
      real(8)    :: vec(3),norm
      logical    :: covar

c     Local variables

      real(8)    :: tensor(3,3)

      type(MG_grid_metrics),pointer :: gmetric 

c     Begin program

      gmetric => g_def%gmetric

      if (covar) then
        tensor = gmetric%grid(igrid)%gsup(i,j,k,:,:)
      else
        tensor = gmetric%grid(igrid)%gsub(i,j,k,:,:)
      endif

      norm = dot_product(vec,matmul(tensor,vec))
     .      /gmetric%grid(igrid)%jac(i,j,k)

c     End 

      end function vectorNorm_ijk

c     vectorNorm_mesh
c     ################################################################
      function vectorNorm_mesh(g_def,igrid,vec,covar) result(norm)

c     ---------------------------------------------------------------
c     Finds L2-norm of vector A given its curvilinear components.
c     ---------------------------------------------------------------

      implicit none

c     Call variables

      type(grid_mg_def),pointer :: g_def

      integer :: igrid
      real(8) :: vec(0:,0:,0:,:)
     .          ,norm(0:size(vec,1)-1,0:size(vec,2)-1,0:size(vec,3)-1)
      logical :: covar

c     Local variables

      integer :: i,j,k,nxl,nyl,nzl

c     Begin program

      nxl = g_def%nxv(igrid)
      nyl = g_def%nyv(igrid)
      nzl = g_def%nzv(igrid)

      do k=0,nzl+1
        do j=0,nyl+1
          do i=0,nxl+1
            norm(i,j,k) = vectorNorm_ijk(g_def,i,j,k,igrid,vec(i,j,k,:)
     .                                  ,covar)
          enddo
        enddo
      enddo

c     End 

      end function vectorNorm_mesh

c     scalarProduct_ijk
c     ################################################################
      function scalarProduct_ijk(g_def,i,j,k,igrid,cov,cnv) result (dot)

c     ---------------------------------------------------------------
c     Finds scalar product of two vectors, one covariant and the
c     other contravariant.
c     ---------------------------------------------------------------

      implicit none

c     Call variables

      type(grid_mg_def),pointer :: g_def

      integer :: i,j,k,igrid
      real(8) :: dot,cnv(3),cov(3)

c     Local variables

      type(MG_grid_metrics),pointer :: gmetric 

c     Begin program

      gmetric => g_def%gmetric

      dot = dot_product(cov,cnv)/gmetric%grid(igrid)%jac(i,j,k)

c     End 

      end function scalarProduct_ijk

c     scalarProduct_mesh
c     ################################################################
      function scalarProduct_mesh(g_def,igrid,vec1,vec2) result (sp)

c     ---------------------------------------------------------------
c     Finds scalar product of two vectors, one covariant and the
c     other contravariant.
c     ---------------------------------------------------------------

      implicit none

c     Call variables

      type(grid_mg_def),pointer :: g_def

      integer :: igrid
      real(8) :: vec1(:,:,:,:),vec2(:,:,:,:)
      real(8) :: sp(size(vec1,1),size(vec1,2),size(vec1,3))

c     Local variables

      type(MG_grid_metrics),pointer :: gmetric 

c     Begin program

      gmetric => g_def%gmetric

      sp = (vec1(:,:,:,1)*vec2(:,:,:,1)
     .     +vec1(:,:,:,2)*vec2(:,:,:,2)
     .     +vec1(:,:,:,3)*vec2(:,:,:,3))/gmetric%grid(igrid)%jac

c     End 

      end function scalarProduct_mesh


c     tensorScalarProduct
c     ################################################################
      function tensorScalarProduct(g_def,i,j,k,igrid,tnsr1,tnsr2)
     .         result (dot)

c     ---------------------------------------------------------------
c     Finds tensor scalar product of two tensors, one covariant
c     and the other contravariant.
c     ---------------------------------------------------------------

      implicit none

c     Call variables

      type(grid_mg_def),pointer :: g_def

      integer    :: i,j,k,igrid
      real(8)    :: dot,tnsr1(3,3),tnsr2(3,3)

c     Local variables

      type(MG_grid_metrics),pointer :: gmetric 

c     Begin program

      gmetric => g_def%gmetric

      dot = sum(tnsr1*tnsr2)/gmetric%grid(igrid)%jac(i,j,k)**2

c     End 

      end function tensorScalarProduct

c     crossProduct_ijk
c     ################################################################
      function crossProduct_ijk(g_def,i,j,k,igrid,vec1,vec2,cnv_in)
     .         result (vec3)

c     ---------------------------------------------------------------
c     Finds cross product of two vectors, either both covariant
c     (cnv_in=.false.) or both contravariant (cnv_in=.true.) Result
c     is covariant if cnv_in=.true., and contravariant otherwise.
c     ---------------------------------------------------------------

      implicit none

c     Call variables

      type(grid_mg_def),pointer :: g_def

      integer    :: i,j,k,igrid
      real(8)    :: vec1(3),vec2(3),vec3(3)
      logical    :: cnv_in

c     Local variables

      real(8)    :: ijac

c     Begin program

      if (cnv_in) then  !Result is covariant
        ijac = 1d0/g_def%gmetric%grid(igrid)%jac(i,j,k)
        vec3(1) = (vec1(2)*vec2(3)
     .            -vec1(3)*vec2(2))*ijac
        vec3(2) = (vec1(3)*vec2(1)
     .            -vec1(1)*vec2(3))*ijac
        vec3(3) = (vec1(1)*vec2(2)
     .            -vec1(2)*vec2(1))*ijac
      else           !Result is contravariant
        vec3(1) = (vec1(2)*vec2(3)
     .            -vec1(3)*vec2(2))
        vec3(2) = (vec1(3)*vec2(1)
     .            -vec1(1)*vec2(3))
        vec3(3) = (vec1(1)*vec2(2)
     .            -vec1(2)*vec2(1))
      endif

c     End 

      end function crossProduct_ijk

c     crossProduct_mesh
c     ################################################################
      function crossProduct_mesh(g_def,igrid,vec1,vec2,cnv_in)
     .         result (vec3)

c     ---------------------------------------------------------------
c     Finds cross product of two vectors, either both covariant
c     (cnv_in=.false.) or both contravariant (cnv_in=.true.) Result
c     is covariant if cnv_in=.true., and contravariant otherwise.
c     ---------------------------------------------------------------

      implicit none

c     Call variables

      type(grid_mg_def),pointer :: g_def

      integer    :: igrid
      real(8)    :: vec1(:,:,:,:),vec2(:,:,:,:)
      real(8)    :: vec3(size(vec1,1),size(vec1,2)
     .                  ,size(vec1,3),size(vec1,4))
      logical    :: cnv_in

c     Local variables

c     Begin program

      if (cnv_in) then  !Result is covariant
        vec3(:,:,:,1) = (vec1(:,:,:,2)*vec2(:,:,:,3)
     .                  -vec1(:,:,:,3)*vec2(:,:,:,2))
     .                  /g_def%gmetric%grid(igrid)%jac
        vec3(:,:,:,2) = (vec1(:,:,:,3)*vec2(:,:,:,1)
     .                  -vec1(:,:,:,1)*vec2(:,:,:,3))
     .                  /g_def%gmetric%grid(igrid)%jac
        vec3(:,:,:,3) = (vec1(:,:,:,1)*vec2(:,:,:,2)
     .                  -vec1(:,:,:,2)*vec2(:,:,:,1))
     .                  /g_def%gmetric%grid(igrid)%jac
      else           !Result is contravariant
        vec3(:,:,:,1) = (vec1(:,:,:,2)*vec2(:,:,:,3)
     .                  -vec1(:,:,:,3)*vec2(:,:,:,2))
        vec3(:,:,:,2) = (vec1(:,:,:,3)*vec2(:,:,:,1)
     .                  -vec1(:,:,:,1)*vec2(:,:,:,3))
        vec3(:,:,:,3) = (vec1(:,:,:,1)*vec2(:,:,:,2)
     .                  -vec1(:,:,:,2)*vec2(:,:,:,1))
      endif

c     End 

      end function crossProduct_mesh

      end module grid_operations_cell

c module grid_operations_stg
c #####################################################################
      module grid_operations_stg

c ---------------------------------------------------------------------
c     This module packs routines that perform operations on staggered
c     vectors, such as coordinate transformation of vector components,
c     vector norms and scalar products.
c ---------------------------------------------------------------------

      use grid_packing

      implicit none

      INTERFACE XformVector_stg
        module procedure XformVector_stg_mesh
      end INTERFACE

      contains

c     XformToCar_stg_ijk
c     #################################################################
      function XformToCar_stg_ijk(g_def,i,j,k,igrid,curv,covariant
     .                           ,mul_by_jac) result(car)

c     -----------------------------------------------------------------
c     Transforms a STAGGERED curvilinear vector (c1,c2,c3) to Cartesian
c     (cx,cy,cz) at grid coordinates (i,j,k). Curvilinear vector is
c     covariant if covariant=.true., and contravariant otherwise.
c     -----------------------------------------------------------------

        implicit none

c     Input variables

        type(grid_mg_def),pointer :: g_def

        integer    :: i,j,k,igrid
        real(8)    :: curv(0:,0:,0:,:),car(3)
        logical    :: covariant

        logical,optional :: mul_by_jac

c     Local variables

        integer    :: im,jm,km
        real(8)    :: cov(3),cnv(3),jac,curv_l(3)
        logical    :: muljac

        type(MG_grid_metrics),pointer :: gmetric 

c     Begin program

        gmetric => g_def%gmetric

        if (PRESENT(mul_by_jac)) then
          muljac = mul_by_jac
        else
          muljac = .false.
        endif

        im = max(i-1,0)
        jm = max(j-1,0)
        km = max(k-1,0)

        if (covariant) then  !Covariant
          !i+1/2,j,k
          curv_l(1) = curv(i,j,k,1)
          curv_l(2) = 0.25*(curv(i  ,jm,k,2)
     $                     +curv(i  ,j ,k,2)
     $                     +curv(i+1,jm,k,2)
     $                     +curv(i+1,j ,k,2))
          curv_l(3) = 0.25*(curv(i  ,j,km,3)
     $                     +curv(i  ,j,k ,3)
     $                     +curv(i+1,j,km,3)
     $                     +curv(i+1,j,k ,3))

          cov = 0.25*(gmetric%grid(igrid)%cov_v(i,jm,k ,:,1)
     $               +gmetric%grid(igrid)%cov_v(i,j ,k ,:,1)
     $               +gmetric%grid(igrid)%cov_v(i,jm,km,:,1)
     $               +gmetric%grid(igrid)%cov_v(i,j ,km,:,1))
          if (muljac) then
            jac = 0.25*(gmetric%grid(igrid)%jac_v(i,jm,k )
     $                 +gmetric%grid(igrid)%jac_v(i,j ,k )
     $                 +gmetric%grid(igrid)%jac_v(i,jm,km)
     $                 +gmetric%grid(igrid)%jac_v(i,j ,km))
            cov = cov*jac
          endif
          car(1) = dot_product(curv_l,cov)

          !i,j+1/2,k
          curv_l(1) = 0.25*(curv(i ,j+1,k,1)
     $                     +curv(i ,j  ,k,1)
     $                     +curv(im,j+1,k,1)
     $                     +curv(im,j  ,k,1))
          curv_l(2) = curv(i,j,k,2)
          curv_l(3) = 0.25*(curv(i,j+1,k ,3)
     $                     +curv(i,j  ,k ,3)
     $                     +curv(i,j+1,km,3)
     $                     +curv(i,j  ,km,3))
          cov = 0.25*(gmetric%grid(igrid)%cov_v(im,j,k ,:,2)
     $               +gmetric%grid(igrid)%cov_v(i ,j,k ,:,2)
     $               +gmetric%grid(igrid)%cov_v(im,j,km,:,2)
     $               +gmetric%grid(igrid)%cov_v(i ,j,km,:,2))
          if (muljac) then
            jac = 0.25*(gmetric%grid(igrid)%jac_v(im,j,k )
     $                 +gmetric%grid(igrid)%jac_v(i ,j,k )
     $                 +gmetric%grid(igrid)%jac_v(im,j,km)
     $                 +gmetric%grid(igrid)%jac_v(i ,j,km))
            cov = cov*jac
          endif
          car(2) = dot_product(curv_l,cov)

          !i,j,k+1/2
          curv_l(1) = 0.25*(curv(i ,j,k+1,1)
     $                     +curv(i ,j,k  ,1)
     $                     +curv(im,j,k+1,1)
     $                     +curv(im,j,k  ,1))
          curv_l(2) = 0.25*(curv(i,j ,k+1,2)
     $                     +curv(i,j ,k  ,2)
     $                     +curv(i,jm,k+1,2)
     $                     +curv(i,jm,k  ,2))
          curv_l(3) = curv(i,j,k,3)
          cov = 0.25*(gmetric%grid(igrid)%cov_v(im,j ,k,:,3)
     $               +gmetric%grid(igrid)%cov_v(i ,j ,k,:,3)
     $               +gmetric%grid(igrid)%cov_v(im,jm,k,:,3)
     $               +gmetric%grid(igrid)%cov_v(i ,jm,k,:,3))
          if (muljac) then
            jac = 0.25*(gmetric%grid(igrid)%jac_v(im,j ,k )
     $                 +gmetric%grid(igrid)%jac_v(i ,j ,k )
     $                 +gmetric%grid(igrid)%jac_v(im,jm,k )
     $                 +gmetric%grid(igrid)%jac_v(i ,jm,k ))
            cov = cov*jac
          endif
          car(3) = dot_product(curv_l,cov)
        else  !Contravariant
          !i+1/2,j,k
          curv_l(1) = curv(i,j,k,1)
          curv_l(2) = 0.25*(curv(i  ,jm,k,2)
     $                     +curv(i  ,j ,k,2)
     $                     +curv(i+1,jm,k,2)
     $                     +curv(i+1,j ,k,2))
          curv_l(3) = 0.25*(curv(i  ,j,km,3)
     $                     +curv(i  ,j,k ,3)
     $                     +curv(i+1,j,km,3)
     $                     +curv(i+1,j,k ,3))
          cnv = 0.25*(gmetric%grid(igrid)%cnv_v(i,jm,k ,:,1)
     $               *gmetric%grid(igrid)%jac_v(i,jm,k  )
     $               +gmetric%grid(igrid)%cnv_v(i,j ,k ,:,1)
     $               *gmetric%grid(igrid)%jac_v(i,j ,k  )
     $               +gmetric%grid(igrid)%cnv_v(i,jm,km,:,1)
     $               *gmetric%grid(igrid)%jac_v(i,jm,km)
     $               +gmetric%grid(igrid)%cnv_v(i,j ,km,:,1)
     $               *gmetric%grid(igrid)%jac_v(i,j ,km))
          if (.not.muljac) then
            jac = 0.25*(gmetric%grid(igrid)%jac_v(i,jm,k )
     $                 +gmetric%grid(igrid)%jac_v(i,j ,k )
     $                 +gmetric%grid(igrid)%jac_v(i,jm,km)
     $                 +gmetric%grid(igrid)%jac_v(i,j ,km))
            cnv = cnv/jac
          endif
          car(1) = dot_product(curv_l,cnv)

          !i,j+1/2,k
          curv_l(1) = 0.25*(curv(i ,j+1,k,1)
     $                     +curv(i ,j  ,k,1)
     $                     +curv(im,j+1,k,1)
     $                     +curv(im,j  ,k,1))
          curv_l(2) = curv(i,j,k,2)
          curv_l(3) = 0.25*(curv(i,j+1,k ,3)
     $                     +curv(i,j  ,k ,3)
     $                     +curv(i,j+1,km,3)
     $                     +curv(i,j  ,km,3))
          cnv = 0.25*(gmetric%grid(igrid)%cnv_v(im,j,k  ,:,2)
     $               *gmetric%grid(igrid)%jac_v(im,j,k )  
     $               +gmetric%grid(igrid)%cnv_v(i ,j,k  ,:,2)
     $               *gmetric%grid(igrid)%jac_v(i ,j,k )  
     $               +gmetric%grid(igrid)%cnv_v(im,j,km ,:,2)
     $               *gmetric%grid(igrid)%jac_v(im,j,km)  
     $               +gmetric%grid(igrid)%cnv_v(i ,j,km ,:,2)
     $               *gmetric%grid(igrid)%jac_v(i ,j,km))
          if (.not.muljac) then
            jac = 0.25*(gmetric%grid(igrid)%jac_v(im,j,k )
     $                 +gmetric%grid(igrid)%jac_v(i ,j,k )
     $                 +gmetric%grid(igrid)%jac_v(im,j,km)
     $                 +gmetric%grid(igrid)%jac_v(i ,j,km))
            cnv = cnv/jac
          endif
          car(2) = dot_product(curv_l,cnv)

          !i,j,k+1/2
          curv_l(1) = 0.25*(curv(i ,j,k+1,1)
     $                     +curv(i ,j,k  ,1)
     $                     +curv(im,j,k+1,1)
     $                     +curv(im,j,k  ,1))
          curv_l(2) = 0.25*(curv(i,j ,k+1,2)
     $                     +curv(i,j ,k  ,2)
     $                     +curv(i,jm,k+1,2)
     $                     +curv(i,jm,k  ,2))
          curv_l(3) = curv(i,j,k,3)
          cnv = 0.25*(gmetric%grid(igrid)%cnv_v(im,j ,k,:,3)
     $               *gmetric%grid(igrid)%jac_v(im,j ,k    )
     $               +gmetric%grid(igrid)%cnv_v(i ,j ,k,:,3)
     $               *gmetric%grid(igrid)%jac_v(i ,j ,k    )
     $               +gmetric%grid(igrid)%cnv_v(im,jm,k,:,3)
     $               *gmetric%grid(igrid)%jac_v(im,jm,k    )
     $               +gmetric%grid(igrid)%cnv_v(i ,jm,k,:,3)
     $               *gmetric%grid(igrid)%jac_v(i ,jm,k    ) )
          if (.not.muljac) then
            jac = 0.25*(gmetric%grid(igrid)%jac_v(im,j ,k )
     $                 +gmetric%grid(igrid)%jac_v(i ,j ,k )
     $                 +gmetric%grid(igrid)%jac_v(im,jm,k )
     $                 +gmetric%grid(igrid)%jac_v(i ,jm,k ))
            cnv = cnv/jac
          endif
          car(3) = dot_product(curv_l,cnv)
        endif

      end function XformToCar_stg_ijk

c     XformToCar_stg_mesh
c     #################################################################
      function XformToCar_stg_mesh(g_def,igrid,curv,covariant
     .                            ,mul_by_jac) result(car)

c     -----------------------------------------------------------------
c     Transforms a staggered curvilinear vector to Cartesian 
c     over the whole mesh. Curvilinear vector is covariant if
c     covariant=.true., and contravariant otherwise.
c     -----------------------------------------------------------------

        implicit none

c     Input variables

        type(grid_mg_def),pointer :: g_def

        integer :: i,j,k,igrid
        real(8) :: curv(0:,0:,0:,:)
     .            ,car(0:size(curv,1)-1
     .                ,0:size(curv,2)-1
     .                ,0:size(curv,3)-1,size(curv,4))
        logical :: covariant

        logical,optional :: mul_by_jac

c     Local variables

        integer :: nxl,nyl,nzl
        real(8) :: curv_l(3),cov_v(3,3),car_l(3)

c     Begin program

        nxl = size(curv,1)-2
        nyl = size(curv,2)-2
        nzl = size(curv,3)-2

        do k = 0,nzl
          do j = 0,nyl
            do i = 0,nxl
              car(i,j,k,:) = XformToCar_stg_ijk(g_def,i,j,k,igrid,curv
     .                                         ,covariant
     .                                         ,mul_by_jac=mul_by_jac)
            enddo
          enddo
        enddo

      end function XformToCar_stg_mesh

c     XformToCurv_stg_ijk
c     #################################################################
      function XformToCurv_stg_ijk(g_def,i,j,k,igrid,car,tocnv)
     .        result(curv)

c     -----------------------------------------------------------------
c     Transforms a cartesian STAGGERED vector (cx,cy,cz) to curvilinear
c     (c1,c2,c3) at grid coordinates (i,j,k). Cartesian vector is
c     contravariant if tocnv=.true., and covariant otherwise.
c     -----------------------------------------------------------------

        implicit none

c     Input variables

        type(grid_mg_def),pointer :: g_def

        integer :: i,j,k,igrid
        real(8) :: curv(3),car(0:,0:,0:,:)
        logical :: tocnv

c     Local variables

        integer    :: im,jm,km
        real(8)    :: cov(3),cnv(3),car_l(3)
        
        type(MG_grid_metrics),pointer :: gmetric 

c     Begin program

        gmetric => g_def%gmetric

        im = max(i-1,0)
        jm = max(j-1,0)
        km = max(k-1,0)

        if (tocnv) then
          !i+1/2,j,k
          car_l(1) = car(i,j,k,1)
          car_l(2) = 0.25*(car(i  ,jm,k,2)
     $                    +car(i  ,j ,k,2)
     $                    +car(i+1,jm,k,2)
     $                    +car(i+1,j ,k,2))
          car_l(3) = 0.25*(car(i  ,j,km,3)
     $                    +car(i  ,j,k ,3)
     $                    +car(i+1,j,km,3)
     $                    +car(i+1,j,k ,3))
c$$$          cov = 0.25*(gmetric%grid(igrid)%cov_v(i,jm,k ,1,:)
c$$$     $               *gmetric%grid(igrid)%jac_v(i,jm,k  )
c$$$     $               +gmetric%grid(igrid)%cov_v(i,j ,k ,1,:)
c$$$     $               *gmetric%grid(igrid)%jac_v(i,j ,k  )
c$$$     $               +gmetric%grid(igrid)%cov_v(i,jm,km,1,:)
c$$$     $               *gmetric%grid(igrid)%jac_v(i,jm,km)
c$$$     $               +gmetric%grid(igrid)%cov_v(i,j ,km,1,:)
c$$$     $               *gmetric%grid(igrid)%jac_v(i,j ,km))
          cov = 0.25*(gmetric%grid(igrid)%cov_v(i,jm,k ,1,:)
     $               +gmetric%grid(igrid)%cov_v(i,j ,k ,1,:)
     $               +gmetric%grid(igrid)%cov_v(i,jm,km,1,:)
     $               +gmetric%grid(igrid)%cov_v(i,j ,km,1,:))
     $         *0.25*(gmetric%grid(igrid)%jac_v(i,jm,k )
     $               +gmetric%grid(igrid)%jac_v(i,j ,k )
     $               +gmetric%grid(igrid)%jac_v(i,jm,km)
     $               +gmetric%grid(igrid)%jac_v(i,j ,km))

          curv(1) = dot_product(cov,car_l)

          !i,j+1/2,k
          car_l(1) = 0.25*(car(i ,j+1,k,1)
     $                    +car(i ,j  ,k,1)
     $                    +car(im,j+1,k,1)
     $                    +car(im,j  ,k,1))
          car_l(2) = car(i,j,k,2)
          car_l(3) = 0.25*(car(i,j+1,k ,3)
     $                    +car(i,j  ,k ,3)
     $                    +car(i,j+1,km,3)
     $                    +car(i,j  ,km,3))
c$$$          cov = 0.25*(gmetric%grid(igrid)%cov_v(im,j,k ,2,:)
c$$$     $               *gmetric%grid(igrid)%jac_v(im,j,k )
c$$$     $               +gmetric%grid(igrid)%cov_v(i ,j,k ,2,:)
c$$$     $               *gmetric%grid(igrid)%jac_v(i ,j,k )
c$$$     $               +gmetric%grid(igrid)%cov_v(im,j,km,2,:)
c$$$     $               *gmetric%grid(igrid)%jac_v(im,j,km)
c$$$     $               +gmetric%grid(igrid)%cov_v(i ,j,km,2,:)
c$$$     $               *gmetric%grid(igrid)%jac_v(i ,j,km))
          cov = 0.25*(gmetric%grid(igrid)%cov_v(im,j,k ,2,:)
     $               +gmetric%grid(igrid)%cov_v(i ,j,k ,2,:)
     $               +gmetric%grid(igrid)%cov_v(im,j,km,2,:)
     $               +gmetric%grid(igrid)%cov_v(i ,j,km,2,:))
     $         *0.25*(gmetric%grid(igrid)%jac_v(im,j,k )
     $               +gmetric%grid(igrid)%jac_v(i ,j,k )
     $               +gmetric%grid(igrid)%jac_v(im,j,km)
     $               +gmetric%grid(igrid)%jac_v(i ,j,km))
          curv(2) = dot_product(cov,car_l)

          !i,j,k+1/2
          car_l(1) = 0.25*(car(i ,j,k+1,1)
     $                    +car(i ,j,k  ,1)
     $                    +car(im,j,k+1,1)
     $                    +car(im,j,k  ,1))
          car_l(2) = 0.25*(car(i,j ,k+1,2)
     $                    +car(i,j ,k  ,2)
     $                    +car(i,jm,k+1,2)
     $                    +car(i,jm,k  ,2))
          car_l(3) = car(i,j,k,3)
c$$$          cov = 0.25*(gmetric%grid(igrid)%cov_v(im,j ,k,3,:)
c$$$     $               *gmetric%grid(igrid)%jac_v(im,j ,k )
c$$$     $               +gmetric%grid(igrid)%cov_v(i ,j ,k,3,:)
c$$$     $               *gmetric%grid(igrid)%jac_v(i ,j ,k )
c$$$     $               +gmetric%grid(igrid)%cov_v(im,jm,k,3,:)
c$$$     $               *gmetric%grid(igrid)%jac_v(im,jm,k )
c$$$     $               +gmetric%grid(igrid)%cov_v(i ,jm,k,3,:))
c$$$     $               *gmetric%grid(igrid)%jac_v(i ,jm,k )
          cov = 0.25*(gmetric%grid(igrid)%cov_v(im,j ,k,3,:)
     $               +gmetric%grid(igrid)%cov_v(i ,j ,k,3,:)
     $               +gmetric%grid(igrid)%cov_v(im,jm,k,3,:)
     $               +gmetric%grid(igrid)%cov_v(i ,jm,k,3,:))
     $         *0.25*(gmetric%grid(igrid)%jac_v(im,j ,k )
     $               +gmetric%grid(igrid)%jac_v(i ,j ,k )
     $               +gmetric%grid(igrid)%jac_v(im,jm,k )
     $               +gmetric%grid(igrid)%jac_v(i ,jm,k ))
          curv(3) = dot_product(cov,car_l)
        else
          !i+1/2,j,k
          car_l(1) = car(i,j,k,1)
          car_l(2) = 0.25*(car(i  ,jm,k,2)
     $                    +car(i  ,j ,k,2)
     $                    +car(i+1,jm,k,2)
     $                    +car(i+1,j ,k,2))
          car_l(3) = 0.25*(car(i  ,j,km,3)
     $                    +car(i  ,j,k ,3)
     $                    +car(i+1,j,km,3)
     $                    +car(i+1,j,k ,3))
c$$$          cnv = 0.25*(gmetric%grid(igrid)%cnv_v(i,jm,k ,1,:)
c$$$     $               *gmetric%grid(igrid)%jac_v(i,jm,k )
c$$$     $               +gmetric%grid(igrid)%cnv_v(i,j ,k ,1,:)
c$$$     $               *gmetric%grid(igrid)%jac_v(i,j ,k )
c$$$     $               +gmetric%grid(igrid)%cnv_v(i,jm,km,1,:)
c$$$     $               *gmetric%grid(igrid)%jac_v(i,jm,km)
c$$$     $               +gmetric%grid(igrid)%cnv_v(i,j ,km,1,:)
c$$$     $               *gmetric%grid(igrid)%jac_v(i,j ,km))
          cnv = 0.25*(gmetric%grid(igrid)%cnv_v(i,jm,k ,1,:)
     $               +gmetric%grid(igrid)%cnv_v(i,j ,k ,1,:)
     $               +gmetric%grid(igrid)%cnv_v(i,jm,km,1,:)
     $               +gmetric%grid(igrid)%cnv_v(i,j ,km,1,:))
     $         *0.25*(gmetric%grid(igrid)%jac_v(i,jm,k )
     $               +gmetric%grid(igrid)%jac_v(i,j ,k )
     $               +gmetric%grid(igrid)%jac_v(i,jm,km)
     $               +gmetric%grid(igrid)%jac_v(i,j ,km))

          curv(1) = dot_product(cnv,car_l)

          !i,j+1/2,k
          car_l(1) = 0.25*(car(i ,j+1,k,1)
     $                    +car(i ,j  ,k,1)
     $                    +car(im,j+1,k,1)
     $                    +car(im,j  ,k,1))
          car_l(2) = car(i,j,k,2)
          car_l(3) = 0.25*(car(i,j+1,k ,3)
     $                    +car(i,j  ,k ,3)
     $                    +car(i,j+1,km,3)
     $                    +car(i,j  ,km,3))
c$$$          cnv = 0.25*(gmetric%grid(igrid)%cnv_v(im,j,k ,2,:)
c$$$     $               *gmetric%grid(igrid)%jac_v(im,j,k )
c$$$     $               +gmetric%grid(igrid)%cnv_v(i ,j,k ,2,:)
c$$$     $               *gmetric%grid(igrid)%jac_v(i ,j,k )
c$$$     $               +gmetric%grid(igrid)%cnv_v(im,j,km,2,:)
c$$$     $               *gmetric%grid(igrid)%jac_v(im,j,km)
c$$$     $               +gmetric%grid(igrid)%cnv_v(i ,j,km,2,:)
c$$$     $               *gmetric%grid(igrid)%jac_v(i ,j,km))
          cnv = 0.25*(gmetric%grid(igrid)%cnv_v(im,j,k ,2,:)
     $               +gmetric%grid(igrid)%cnv_v(i ,j,k ,2,:)
     $               +gmetric%grid(igrid)%cnv_v(im,j,km,2,:)
     $               +gmetric%grid(igrid)%cnv_v(i ,j,km,2,:))
     $         *0.25*(gmetric%grid(igrid)%jac_v(im,j,k )
     $               +gmetric%grid(igrid)%jac_v(i ,j,k )
     $               +gmetric%grid(igrid)%jac_v(im,j,km)
     $               +gmetric%grid(igrid)%jac_v(i ,j,km))

          curv(2) = dot_product(cnv,car_l)

          !i,j,k+1/2
          car_l(1) = 0.25*(car(i ,j,k+1,1)
     $                    +car(i ,j,k  ,1)
     $                    +car(im,j,k+1,1)
     $                    +car(im,j,k  ,1))
          car_l(2) = 0.25*(car(i,j ,k+1,2)
     $                    +car(i,j ,k  ,2)
     $                    +car(i,jm,k+1,2)
     $                    +car(i,jm,k  ,2))
          car_l(3) = car(i,j,k,3)
c$$$          cnv = 0.25*(gmetric%grid(igrid)%cnv_v(im,j ,k,3,:)
c$$$     $               *gmetric%grid(igrid)%jac_v(im,j ,k  )
c$$$     $               +gmetric%grid(igrid)%cnv_v(i ,j ,k,3,:)
c$$$     $               *gmetric%grid(igrid)%jac_v(i ,j ,k  )
c$$$     $               +gmetric%grid(igrid)%cnv_v(im,jm,k,3,:)
c$$$     $               *gmetric%grid(igrid)%jac_v(im,jm,k  )
c$$$     $               +gmetric%grid(igrid)%cnv_v(i ,jm,k,3,:))
c$$$     $               *gmetric%grid(igrid)%jac_v(i ,jm,k  )
          cnv = 0.25*(gmetric%grid(igrid)%cnv_v(im,j ,k,3,:)
     $               +gmetric%grid(igrid)%cnv_v(i ,j ,k,3,:)
     $               +gmetric%grid(igrid)%cnv_v(im,jm,k,3,:)
     $               +gmetric%grid(igrid)%cnv_v(i ,jm,k,3,:))
     $         *0.25*(gmetric%grid(igrid)%jac_v(im,j ,k  )
     $               +gmetric%grid(igrid)%jac_v(i ,j ,k  )
     $               +gmetric%grid(igrid)%jac_v(im,jm,k  )
     $               +gmetric%grid(igrid)%jac_v(i ,jm,k  ))
          curv(3) = dot_product(cnv,car_l)
        endif

      end function XformToCurv_stg_ijk

c     XformToCurv_stg_mesh
c     #################################################################
      function XformToCurv_stg_mesh(g_def,igrid,car,tocnv) result(curv)

c     -----------------------------------------------------------------
c     Transforms a staggered Cartesian vector to curvilinear
c     over the whole mesh. Curvilinear vector is contravariant if
c     tocnv=.true., and covariant otherwise.
c     -----------------------------------------------------------------

        implicit none

c     Input variables

        type(grid_mg_def),pointer :: g_def

        integer :: igrid
        real(8) :: car(0:,0:,0:,:)
     .            ,curv(0:size(car,1)-1
     .                 ,0:size(car,2)-1
     .                 ,0:size(car,3)-1,size(car,4))
        logical :: tocnv

c     Local variables

        integer :: i,j,k,nxl,nyl,nzl

c     Begin program

        nxl = size(car,1)-2
        nyl = size(car,2)-2
        nzl = size(car,3)-2

        curv = 0d0
        
        do k = 0,nzl
          do j = 0,nyl
            do i = 0,nxl
              curv(i,j,k,:) = XformToCurv_stg_ijk(g_def,i,j,k,igrid,car
     .                                           ,tocnv)
            enddo
          enddo
        enddo

      end function XformToCurv_stg_mesh

c$$$c     XformToCurv_stg_mesh2
c$$$c     #################################################################
c$$$      subroutine XformToCurv_stg_mesh2(g_def,igrid,vec,tocnv)
c$$$
c$$$c     -----------------------------------------------------------------
c$$$c     Transforms a staggered Cartesian vector to curvilinear
c$$$c     over the whole mesh. Curvilinear vector is contravariant if
c$$$c     tocnv=.true., and covariant otherwise.
c$$$c     -----------------------------------------------------------------
c$$$
c$$$        implicit none
c$$$
c$$$c     Input variables
c$$$
c$$$        type(grid_mg_def),pointer :: g_def
c$$$
c$$$        integer :: igrid
c$$$        real(8) :: vec(0:,0:,0:,:)
c$$$        logical :: tocnv
c$$$
c$$$c     Local variables
c$$$
c$$$c     Begin program
c$$$
c$$$        vec = XformToCurv_stg_mesh(g_def,igrid,vec,tocnv)
c$$$
c$$$      end subroutine XformToCurv_stg_mesh2

c     XformToCnv_stg_ijk
c     #################################################################
      function XformToCnv_stg_ijk(g_def,i,j,k,igrid,cov) result(cnv)

c     -----------------------------------------------------------------
c     Transforms a staggered covariant vector to contravariant
c     representation at faces for cell (i,j,k).
c     -----------------------------------------------------------------

        implicit none

c     Input variables

        type(grid_mg_def),pointer :: g_def

        integer :: i,j,k,igrid
        real(8) :: cov(0:,0:,0:,:),cnv(3)

c     Local variables

        integer :: nx,ny,nz,ip,im,jp,jm,kp,km

        real(8) :: d_xx_ip,d_xy_ipjp,d_xy_ipjm,d_xz_ipkp,d_xz_ipkm
        real(8) :: d_yy_jp,d_xy_imjp,d_yz_jpkp,d_yz_jpkm
        real(8) :: d_zz_kp,d_xz_imkp,d_yz_jmkp

c     Begin program

        cnv = 0d0

        ip = min(i+1,size(cov,1)-1)
        im = i-1

        jp = min(j+1,size(cov,2)-1)
        jm = j-1

        kp = min(k+1,size(cov,3)-1)
        km = k-1

        if (jm >= 0.and.km >= 0) then
          !i+1/2,j,k
          d_xx_ip = 0.25*(g_def%gmetric%grid(igrid)%gsup_v(i,j ,k ,1,1)
     .                   +g_def%gmetric%grid(igrid)%gsup_v(i,jm,k ,1,1)
     .                   +g_def%gmetric%grid(igrid)%gsup_v(i,j ,km,1,1)
     .                   +g_def%gmetric%grid(igrid)%gsup_v(i,jm,km,1,1))

          d_xy_ipjp =0.5*(g_def%gmetric%grid(igrid)%gsup_v(i,j ,k ,1,2)
     .                   +g_def%gmetric%grid(igrid)%gsup_v(i,j ,km,1,2))
          d_xy_ipjm =0.5*(g_def%gmetric%grid(igrid)%gsup_v(i,jm,k ,1,2)
     .                   +g_def%gmetric%grid(igrid)%gsup_v(i,jm,km,1,2))

          d_xz_ipkp =0.5*(g_def%gmetric%grid(igrid)%gsup_v(i,j ,k ,1,3)
     .                   +g_def%gmetric%grid(igrid)%gsup_v(i,jm,k ,1,3))
          d_xz_ipkm =0.5*(g_def%gmetric%grid(igrid)%gsup_v(i,j ,km,1,3)
     .                   +g_def%gmetric%grid(igrid)%gsup_v(i,jm,km,1,3))

          cnv(1) = d_xx_ip*cov(i,j,k,1)
     .           +0.25*d_xy_ipjp*(cov(i ,j ,k ,2)
     .                           +cov(ip,j ,k ,2))
     .           +0.25*d_xy_ipjm*(cov(i ,jm,k ,2)
     .                           +cov(ip,jm,k ,2))
     .           +0.25*d_xz_ipkp*(cov(i ,j ,k ,3)
     .                           +cov(ip,j ,k ,3))
     .           +0.25*d_xz_ipkm*(cov(i ,j ,km,3)
     .                           +cov(ip,j ,km,3))
        endif

        if (im >= 0.and.km >= 0) then
          !i,j+1/2,k
          d_xy_ipjp =0.5*(g_def%gmetric%grid(igrid)%gsup_v(i ,j,k ,1,2)
     .                   +g_def%gmetric%grid(igrid)%gsup_v(i ,j,km,1,2))
          d_xy_imjp =0.5*(g_def%gmetric%grid(igrid)%gsup_v(im,j,k ,1,2)
     .                   +g_def%gmetric%grid(igrid)%gsup_v(im,j,km,1,2))

          d_yy_jp = 0.25*(g_def%gmetric%grid(igrid)%gsup_v(i ,j,k ,2,2)
     .                   +g_def%gmetric%grid(igrid)%gsup_v(im,j,k ,2,2)
     .                   +g_def%gmetric%grid(igrid)%gsup_v(i ,j,km,2,2)
     .                   +g_def%gmetric%grid(igrid)%gsup_v(im,j,km,2,2))

          d_yz_jpkp =0.5*(g_def%gmetric%grid(igrid)%gsup_v(i ,j,k ,2,3)
     .                   +g_def%gmetric%grid(igrid)%gsup_v(im,j,k ,2,3))
          d_yz_jpkm =0.5*(g_def%gmetric%grid(igrid)%gsup_v(i ,j,km,2,3)
     .                   +g_def%gmetric%grid(igrid)%gsup_v(im,j,km,2,3))

          cnv(2) = d_yy_jp*cov(i,j,k,2)
     .        +0.25*d_xy_ipjp*(cov(i ,j ,k ,1)
     .                        +cov(i ,jp,k ,1))
     .        +0.25*d_xy_imjp*(cov(im,j ,k ,1)
     .                        +cov(im,jp,k ,1))
     .        +0.25*d_yz_jpkp*(cov(i ,j ,k ,3)
     .                        +cov(i ,jp,k ,3))
     .        +0.25*d_yz_jpkm*(cov(i ,j ,km,3)
     .                        +cov(i ,jp,km,3))
        endif

        if (im >= 0.and.jm >= 0) then
          !i,j,k+1/2
          d_xz_ipkp =0.5*(g_def%gmetric%grid(igrid)%gsup_v(i ,j ,k,1,3)
     .                   +g_def%gmetric%grid(igrid)%gsup_v(i ,jm,k,1,3))
          d_xz_imkp =0.5*(g_def%gmetric%grid(igrid)%gsup_v(im,j ,k,1,3)
     .                   +g_def%gmetric%grid(igrid)%gsup_v(im,jm,k,1,3))

          d_yz_jpkp =0.5*(g_def%gmetric%grid(igrid)%gsup_v(i ,j ,k,2,3)
     .                   +g_def%gmetric%grid(igrid)%gsup_v(im,j ,k,2,3))
          d_yz_jmkp =0.5*(g_def%gmetric%grid(igrid)%gsup_v(i ,jm,k,2,3)
     .                   +g_def%gmetric%grid(igrid)%gsup_v(im,jm,k,2,3))

          d_zz_kp = 0.25*(g_def%gmetric%grid(igrid)%gsup_v(i ,j ,k,3,3)
     .                   +g_def%gmetric%grid(igrid)%gsup_v(im,j ,k,3,3)
     .                   +g_def%gmetric%grid(igrid)%gsup_v(i ,jm,k,3,3)
     .                   +g_def%gmetric%grid(igrid)%gsup_v(im,jm,k,3,3))

          cnv(3) = d_zz_kp*cov(i,j,k,3)
     .        +0.25*d_xz_ipkp*(cov(i ,j ,k ,1)
     .                        +cov(i ,j ,kp,1))
     .        +0.25*d_xz_imkp*(cov(im,j ,k ,1)
     .                        +cov(im,j ,kp,1))
     .        +0.25*d_yz_jpkp*(cov(i ,j ,k ,2)
     .                        +cov(i ,j ,kp,2))
     .        +0.25*d_yz_jmkp*(cov(i ,jm,k ,2)
     .                        +cov(i ,jm,kp,2))
        endif

      end function XformToCnv_stg_ijk

c     XformToCov_stg_ijk
c     #################################################################
      function XformToCov_stg_ijk(g_def,i,j,k,igrid,cnv) result(cov)

c     -----------------------------------------------------------------
c     Transforms a staggered contrariant vector to covariant
c     representation at faces for cell (i,j,k).
c     -----------------------------------------------------------------

        implicit none

c     Input variables

        type(grid_mg_def),pointer :: g_def

        integer :: i,j,k,igrid
        real(8) :: cnv(0:,0:,0:,:),cov(3)

c     Local variables

        integer :: nx,ny,nz,ip,im,jp,jm,kp,km

        real(8) :: d_xx_ip,d_xy_ipjp,d_xy_ipjm,d_xz_ipkp,d_xz_ipkm
        real(8) :: d_yy_jp,d_xy_imjp,d_yz_jpkp,d_yz_jpkm
        real(8) :: d_zz_kp,d_xz_imkp,d_yz_jmkp

c     Begin program

        cov = 0d0

        ip = min(i+1,size(cnv,1)-1)
        im = i-1

        jp = min(j+1,size(cnv,2)-1)
        jm = j-1

        kp = min(k+1,size(cnv,3)-1)
        km = k-1

        if (jm >= 0.and.km >= 0) then
          !i+1/2,j,k
          d_xx_ip = 0.25*(g_def%gmetric%grid(igrid)%gsub_v(i,j ,k ,1,1)
     .                   +g_def%gmetric%grid(igrid)%gsub_v(i,jm,k ,1,1)
     .                   +g_def%gmetric%grid(igrid)%gsub_v(i,j ,km,1,1)
     .                   +g_def%gmetric%grid(igrid)%gsub_v(i,jm,km,1,1))

          d_xy_ipjp =0.5*(g_def%gmetric%grid(igrid)%gsub_v(i,j ,k ,1,2)
     .                   +g_def%gmetric%grid(igrid)%gsub_v(i,j ,km,1,2))
          d_xy_ipjm =0.5*(g_def%gmetric%grid(igrid)%gsub_v(i,jm,k ,1,2)
     .                   +g_def%gmetric%grid(igrid)%gsub_v(i,jm,km,1,2))

          d_xz_ipkp =0.5*(g_def%gmetric%grid(igrid)%gsub_v(i,j ,k ,1,3)
     .                   +g_def%gmetric%grid(igrid)%gsub_v(i,jm,k ,1,3))
          d_xz_ipkm =0.5*(g_def%gmetric%grid(igrid)%gsub_v(i,j ,km,1,3)
     .                   +g_def%gmetric%grid(igrid)%gsub_v(i,jm,km,1,3))

          cov(1) = d_xx_ip*cnv(i,j,k,1)
     .           +0.25*d_xy_ipjp*(cnv(i ,j ,k ,2)
     .                           +cnv(ip,j ,k ,2))
     .           +0.25*d_xy_ipjm*(cnv(i ,jm,k ,2)
     .                           +cnv(ip,jm,k ,2))
     .           +0.25*d_xz_ipkp*(cnv(i ,j ,k ,3)
     .                           +cnv(ip,j ,k ,3))
     .           +0.25*d_xz_ipkm*(cnv(i ,j ,km,3)
     .                           +cnv(ip,j ,km,3))
        endif

        if (im >= 0.and.km >= 0) then
          !i,j+1/2,k
          d_xy_ipjp =0.5*(g_def%gmetric%grid(igrid)%gsub_v(i ,j,k ,1,2)
     .                   +g_def%gmetric%grid(igrid)%gsub_v(i ,j,km,1,2))
          d_xy_imjp =0.5*(g_def%gmetric%grid(igrid)%gsub_v(im,j,k ,1,2)
     .                   +g_def%gmetric%grid(igrid)%gsub_v(im,j,km,1,2))

          d_yy_jp = 0.25*(g_def%gmetric%grid(igrid)%gsub_v(i ,j,k ,2,2)
     .                   +g_def%gmetric%grid(igrid)%gsub_v(im,j,k ,2,2)
     .                   +g_def%gmetric%grid(igrid)%gsub_v(i ,j,km,2,2)
     .                   +g_def%gmetric%grid(igrid)%gsub_v(im,j,km,2,2))

          d_yz_jpkp =0.5*(g_def%gmetric%grid(igrid)%gsub_v(i ,j,k ,2,3)
     .                   +g_def%gmetric%grid(igrid)%gsub_v(im,j,k ,2,3))
          d_yz_jpkm =0.5*(g_def%gmetric%grid(igrid)%gsub_v(i ,j,km,2,3)
     .                   +g_def%gmetric%grid(igrid)%gsub_v(im,j,km,2,3))

          cov(2) = d_yy_jp*cnv(i,j,k,2)
     .        +0.25*d_xy_ipjp*(cnv(i ,j ,k ,1)
     .                        +cnv(i ,jp,k ,1))
     .        +0.25*d_xy_imjp*(cnv(im,j ,k ,1)
     .                        +cnv(im,jp,k ,1))
     .        +0.25*d_yz_jpkp*(cnv(i ,j ,k ,3)
     .                        +cnv(i ,jp,k ,3))
     .        +0.25*d_yz_jpkm*(cnv(i ,j ,km,3)
     .                        +cnv(i ,jp,km,3))
        endif

        if (im >= 0.and.jm >= 0) then
          !i,j,k+1/2
          d_xz_ipkp =0.5*(g_def%gmetric%grid(igrid)%gsub_v(i ,j ,k,1,3)
     .                   +g_def%gmetric%grid(igrid)%gsub_v(i ,jm,k,1,3))
          d_xz_imkp =0.5*(g_def%gmetric%grid(igrid)%gsub_v(im,j ,k,1,3)
     .                   +g_def%gmetric%grid(igrid)%gsub_v(im,jm,k,1,3))

          d_yz_jpkp =0.5*(g_def%gmetric%grid(igrid)%gsub_v(i ,j ,k,2,3)
     .                   +g_def%gmetric%grid(igrid)%gsub_v(im,j ,k,2,3))
          d_yz_jmkp =0.5*(g_def%gmetric%grid(igrid)%gsub_v(i ,jm,k,2,3)
     .                   +g_def%gmetric%grid(igrid)%gsub_v(im,jm,k,2,3))

          d_zz_kp = 0.25*(g_def%gmetric%grid(igrid)%gsub_v(i ,j ,k,3,3)
     .                   +g_def%gmetric%grid(igrid)%gsub_v(im,j ,k,3,3)
     .                   +g_def%gmetric%grid(igrid)%gsub_v(i ,jm,k,3,3)
     .                   +g_def%gmetric%grid(igrid)%gsub_v(im,jm,k,3,3))

          cov(3) = d_zz_kp*cnv(i,j,k,3)
     .        +0.25*d_xz_ipkp*(cnv(i ,j ,k ,1)
     .                        +cnv(i ,j ,kp,1))
     .        +0.25*d_xz_imkp*(cnv(im,j ,k ,1)
     .                        +cnv(im,j ,kp,1))
     .        +0.25*d_yz_jpkp*(cnv(i ,j ,k ,2)
     .                        +cnv(i ,j ,kp,2))
     .        +0.25*d_yz_jmkp*(cnv(i ,jm,k ,2)
     .                        +cnv(i ,jm,kp,2))
        endif

      end function XformToCov_stg_ijk

c$$$c     XformToCov_stg_ijk
c$$$c     #################################################################
c$$$      function XformToCov_stg_ijk(g_def,i,j,k,igrid,cnv) result(cov)
c$$$
c$$$c     -----------------------------------------------------------------
c$$$c     Transforms a staggered contravariant vector to a covariant
c$$$c     representation at faces for cell (i,j,k).
c$$$c     -----------------------------------------------------------------
c$$$
c$$$        implicit none
c$$$
c$$$c     Input variables
c$$$
c$$$        type(grid_mg_def),pointer :: g_def
c$$$
c$$$        integer :: i,j,k,igrid
c$$$        real(8) :: cnv(0:,0:,0:,:),cov(3)
c$$$
c$$$c     Local variables
c$$$
c$$$        integer :: nx,ny,nz,ip,im,jp,jm,kp,km
c$$$        real(8) :: gsub11,gsub22,gsub33
c$$$        real(8),pointer,dimension(:,:,:,:,:) :: gsub_v
c$$$
c$$$        type(MG_grid_metrics),pointer :: gmetric 
c$$$
c$$$c     Begin program
c$$$
c$$$        gmetric => g_def%gmetric
c$$$
c$$$        cov = 0d0
c$$$
c$$$        ip = min(i+1,size(cnv,1)-1)
c$$$        im = i-1
c$$$
c$$$        jp = min(j+1,size(cnv,2)-1)
c$$$        jm = j-1
c$$$
c$$$        kp = min(k+1,size(cnv,3)-1)
c$$$        km = k-1
c$$$
c$$$        gsub_v => gmetric%grid(igrid)%gsub_v
c$$$
c$$$c$$$        if (jm >= 0) then
c$$$c$$$          gsub11 = 0.5*(gsub_v(i,jm,k,1,1)
c$$$c$$$     .                 +gsub_v(i,j ,k,1,1))
c$$$c$$$
c$$$c$$$          cov(1) = gsub11*cnv(i,j,k,1)
c$$$c$$$     .                +0.25*(gsub_v(i,j ,k,1,2)*(cnv(i ,j ,k,2)
c$$$c$$$     .                                          +cnv(ip,j ,k,2))
c$$$c$$$     .                      +gsub_v(i,jm,k,1,2)*(cnv(i ,jm,k,2)
c$$$c$$$     .                                          +cnv(ip,jm,k,2)))
c$$$c$$$        endif
c$$$        if (jm >= 0.and.km >= 0) then
c$$$          !i+1/2,j,k
c$$$          gsub11 = 0.25*(gsub_v(i,jm,k ,1,1)
c$$$     .                  +gsub_v(i,j ,k ,1,1)
c$$$     .                  +gsub_v(i,jm,km,1,1)
c$$$     .                  +gsub_v(i,j ,km,1,1))
c$$$
c$$$          cov(1) = gsub11*cnv(i,j,k,1)
c$$$     .            +0.25*(0.25*gsub_v(i,j ,k ,1,2)*(cnv(i ,j ,k ,2)
c$$$     .                                            +cnv(ip,j ,k ,2)
c$$$     .                                            +cnv(i ,j ,kp,2)
c$$$     .                                            +cnv(ip,j ,kp,2))
c$$$     .                  +0.25*gsub_v(i,jm,k ,1,2)*(cnv(i ,jm,k ,2)
c$$$     .                                            +cnv(ip,jm,k ,2)
c$$$     .                                            +cnv(i ,jm,kp,2)
c$$$     .                                            +cnv(ip,jm,kp,2))
c$$$     .                  +0.25*gsub_v(i,j ,km,1,2)*(cnv(i ,j ,k ,2)
c$$$     .                                            +cnv(ip,j ,k ,2)
c$$$     .                                            +cnv(i ,j ,km,2)
c$$$     .                                            +cnv(ip,j ,km,2))
c$$$     .                  +0.25*gsub_v(i,jm,km,1,2)*(cnv(i ,jm,k ,2)
c$$$     .                                            +cnv(ip,jm,k ,2)
c$$$     .                                            +cnv(i ,jm,km,2)
c$$$     .                                            +cnv(ip,jm,km,2)))
c$$$     .            +0.25*(0.25*gsub_v(i,j ,k ,1,3)*(cnv(i ,j ,k ,3)
c$$$     .                                            +cnv(ip,j ,k ,3)
c$$$     .                                            +cnv(i ,jp,k ,3)
c$$$     .                                            +cnv(ip,jp,k ,3))
c$$$     .                  +0.25*gsub_v(i,j ,km,1,3)*(cnv(i ,j ,km,3)
c$$$     .                                            +cnv(ip,j ,km,3)
c$$$     .                                            +cnv(i ,jp,km,3)
c$$$     .                                            +cnv(ip,jp,km,3))
c$$$     .                  +0.25*gsub_v(i,jm,k ,1,3)*(cnv(i ,j ,k ,3)
c$$$     .                                            +cnv(ip,j ,k ,3)
c$$$     .                                            +cnv(i ,jm,k ,3)
c$$$     .                                            +cnv(ip,jm,k ,3))
c$$$     .                  +0.25*gsub_v(i,jm,km,1,3)*(cnv(i ,j ,km,3)
c$$$     .                                            +cnv(ip,j ,km,3)
c$$$     .                                            +cnv(i ,jm,km,3)
c$$$     .                                            +cnv(ip,jm,km,3)))
c$$$        endif
c$$$
c$$$c$$$        if (im >= 0) then
c$$$c$$$          gsub22 = 0.5*(gsub_v(im,j,k,2,2)
c$$$c$$$     .                 +gsub_v(i ,j,k,2,2))
c$$$c$$$
c$$$c$$$          cov(2) = gsub22*cnv(i,j,k,2)
c$$$c$$$     .                +0.25*(gsub_v(i ,j,k,2,1)*(cnv(i ,j ,k,1)
c$$$c$$$     .                                          +cnv(i ,jp,k,1))
c$$$c$$$     .                      +gsub_v(im,j,k,2,1)*(cnv(im,j ,k,1)
c$$$c$$$     .                                          +cnv(im,jp,k,1)))
c$$$c$$$        endif
c$$$        if (im >= 0.and.km >= 0) then
c$$$          !i,j+1/2,k
c$$$          gsub22 = 0.25*(gsub_v(im,j,k ,2,2)
c$$$     .                  +gsub_v(i ,j,k ,2,2)
c$$$     .                  +gsub_v(im,j,km,2,2)
c$$$     .                  +gsub_v(i ,j,km,2,2))
c$$$
c$$$          cov(2) = gsub22*cnv(i,j,k,2)
c$$$     .            +0.25*(0.25*gsub_v(i ,j,k ,2,1)*(cnv(i ,j ,k ,1)
c$$$     .                                            +cnv(i ,jp,k ,1)
c$$$     .                                            +cnv(i ,j ,kp,1)
c$$$     .                                            +cnv(i ,jp,kp,1))
c$$$     .                  +0.25*gsub_v(im,j,k ,2,1)*(cnv(im,j ,k ,1)
c$$$     .                                            +cnv(im,jp,k ,1)
c$$$     .                                            +cnv(im,j ,kp,1)
c$$$     .                                            +cnv(im,jp,kp,1))
c$$$     .                  +0.25*gsub_v(i ,j,km,2,1)*(cnv(i ,j ,k ,1)
c$$$     .                                            +cnv(i ,jp,k ,1)
c$$$     .                                            +cnv(i ,j ,km,1)
c$$$     .                                            +cnv(i ,jp,km,1))
c$$$     .                  +0.25*gsub_v(im,j,km,2,1)*(cnv(im,j ,k ,1)
c$$$     .                                            +cnv(im,jp,k ,1)
c$$$     .                                            +cnv(im,j ,km,1)
c$$$     .                                            +cnv(im,jp,km,1)))
c$$$     .            +0.25*(0.25*gsub_v(i ,j,k ,2,3)*(cnv(i ,j ,k ,3)
c$$$     .                                            +cnv(ip,j ,k ,3)
c$$$     .                                            +cnv(i ,jp,k ,3)
c$$$     .                                            +cnv(ip,jp,k ,3))
c$$$     .                  +0.25*gsub_v(i ,j,km,2,3)*(cnv(i ,j ,km,3)
c$$$     .                                            +cnv(ip,j ,km,3)
c$$$     .                                            +cnv(i ,jp,km,3)
c$$$     .                                            +cnv(ip,jp,km,3))
c$$$     .                  +0.25*gsub_v(im,j,k ,2,3)*(cnv(i ,j ,k ,3)
c$$$     .                                            +cnv(im,j ,k ,3)
c$$$     .                                            +cnv(i ,jp,k ,3)
c$$$     .                                            +cnv(im,jp,k ,3))
c$$$     .                  +0.25*gsub_v(im,j,km,2,3)*(cnv(i ,j ,km,3)
c$$$     .                                            +cnv(im,j ,km,3)
c$$$     .                                            +cnv(i ,jp,km,3)
c$$$     .                                            +cnv(im,jp,km,3)))
c$$$        endif
c$$$
c$$$cc        cov(3) = gmetric%grid(igrid)%gsub(i,j,k,3,3)*cnv(i,j,k,3)
c$$$        if (im >= 0.and.jm >= 0) then
c$$$          !i,j,k+1/2
c$$$          gsub33 = 0.25*(gsub_v(im,j ,k,3,3)
c$$$     .                  +gsub_v(i ,j ,k,3,3)
c$$$     .                  +gsub_v(im,jm,k,3,3)
c$$$     .                  +gsub_v(i ,jm,k,3,3))
c$$$
c$$$          cov(3) = gsub33*cnv(i,j,k,3)
c$$$     .            +0.25*(0.25*gsub_v(i ,j ,k,3,1)*(cnv(i ,j ,k ,1)
c$$$     .                                            +cnv(i ,j ,kp,1)
c$$$     .                                            +cnv(i ,jp,k ,1)
c$$$     .                                            +cnv(i ,jp,kp,1))
c$$$     .                  +0.25*gsub_v(im,j ,k,3,1)*(cnv(im,j ,k ,1)
c$$$     .                                            +cnv(im,j ,kp,1)
c$$$     .                                            +cnv(im,jp,k ,1)
c$$$     .                                            +cnv(im,jp,kp,1))
c$$$     .                  +0.25*gsub_v(i ,jm,k,3,1)*(cnv(i ,j ,k ,1)
c$$$     .                                            +cnv(i ,j ,kp,1)
c$$$     .                                            +cnv(i ,jm,k ,1)
c$$$     .                                            +cnv(i ,jm,kp,1))
c$$$     .                  +0.25*gsub_v(im,jm,k,3,1)*(cnv(im,j ,k ,1)
c$$$     .                                            +cnv(im,j ,kp,1)
c$$$     .                                            +cnv(im,jm,k ,1)
c$$$     .                                            +cnv(im,jm,kp,1)))
c$$$     .            +0.25*(0.25*gsub_v(i ,j ,k,3,2)*(cnv(i ,j ,k ,2)
c$$$     .                                            +cnv(ip,j ,k ,2)
c$$$     .                                            +cnv(i ,j ,kp,2)
c$$$     .                                            +cnv(ip,j ,kp,2))
c$$$     .                  +0.25*gsub_v(i ,jm,k,3,2)*(cnv(i ,jm,k ,2)
c$$$     .                                            +cnv(ip,jm,k ,2)
c$$$     .                                            +cnv(i ,jm,kp,2)
c$$$     .                                            +cnv(ip,jm,kp,2))
c$$$     .                  +0.25*gsub_v(im,j ,k,3,2)*(cnv(i ,j ,k ,2)
c$$$     .                                            +cnv(im,j ,k ,2)
c$$$     .                                            +cnv(i ,j ,kp,2)
c$$$     .                                            +cnv(im,j ,kp,2))
c$$$     .                  +0.25*gsub_v(im,jm,k,3,2)*(cnv(i ,jm,k ,2)
c$$$     .                                            +cnv(im,jm,k ,2)
c$$$     .                                            +cnv(i ,jm,kp,2)
c$$$     .                                            +cnv(im,jm,kp,2)))
c$$$        endif
c$$$
c$$$      end function XformToCov_stg_ijk

c     XformToCov_vrtx_stg_ijk
c     #################################################################
      function XformToCov_vrtx_stg_ijk(g_def,i,j,k,igrid,cnv)result(cov)

c     -----------------------------------------------------------------
c     Transforms a vertex staggered contravariant vector to a covariant
c     vertex staggered representation at vertex (i+1/2,j+1/2,k+1/2).
c     -----------------------------------------------------------------

        implicit none

c     Input variables

        type(grid_mg_def),pointer :: g_def

        integer :: i,j,k,igrid
        real(8) :: cnv(0:,0:,0:,:),cov(3)

c     Local variables

        integer :: nx,ny,nz,ip,im,jp,jm,kp,km
        real(8) :: gsub11,gsub12,gsub13
     .            ,gsub21,gsub22,gsub23
     .            ,gsub31,gsub32,gsub33
        real(8),pointer,dimension(:,:,:,:,:) :: gsub_v

c     Begin program

        cov = 0d0

        ip = min(i+1,size(cnv,1)-1)
        im = i-1

        jp = min(j+1,size(cnv,2)-1)
        jm = j-1

        kp = min(k+1,size(cnv,3)-1)
        km = k-1

        gsub_v => g_def%gmetric%grid(igrid)%gsub_v

        if (im >= 0) then
          !i,j+1/2,k+1/2
          gsub11 = 0.5*(gsub_v(i ,j,k,1,1)
     .                 +gsub_v(im,j,k,1,1))
c$$$          gsub12 = 0.5*(gsub_v(i ,j,k,1,2)
c$$$     .                 +gsub_v(im,j,k,1,2))
c$$$          gsub13 = 0.5*(gsub_v(i ,j,k,1,3)
c$$$     .                 +gsub_v(im,j,k,1,3))

          cov(1) = gsub11      *cnv(i ,j ,k,1)
c$$$     .            +gsub12*0.25*(cnv(i ,j ,k,2)
c$$$     .                         +cnv(i ,jp,k,2)
c$$$     .                         +cnv(im,j ,k,2)
c$$$     .                         +cnv(im,jp,k,2))
c$$$     .            +gsub13*0.25*(cnv(i ,j,k ,3)
c$$$     .                         +cnv(i ,j,kp,3)
c$$$     .                         +cnv(im,j,k ,3)
c$$$     .                         +cnv(im,j,kp,3))
     .                +0.25*(gsub_v(i ,j,k,1,2)*(cnv(i ,j ,k,2)
     .                                          +cnv(i ,jp,k,2))
     .                      +gsub_v(im,j,k,1,2)*(cnv(im,j ,k,2)
     .                                          +cnv(im,jp,k,2)))
     .                +0.25*(gsub_v(i ,j,k,1,3)*(cnv(i ,j,k ,3)
     .                                          +cnv(i ,j,kp,3))
     .                      +gsub_v(im,j,k,1,3)*(cnv(im,j,k ,3)
     .                                          +cnv(im,j,kp,3)))
        endif

        if (jm >= 0) then
          !i+1/2,j,k+1/2
c$$$          gsub21 = 0.5*(gsub_v(i,jm,k,2,1)
c$$$     .                 +gsub_v(i,j ,k,2,1))
          gsub22 = 0.5*(gsub_v(i,jm,k,2,2)
     .                 +gsub_v(i,j ,k,2,2))
c$$$          gsub23 = 0.5*(gsub_v(i,jm,k,2,3)
c$$$     .                 +gsub_v(i,j ,k,2,3))

          cov(2) = gsub22      *cnv(i ,j ,k,2)
c$$$     .            +gsub21*0.25*(cnv(i ,j ,k,1)
c$$$     .                         +cnv(ip,j ,k,1)
c$$$     .                         +cnv(i ,jm,k,1)
c$$$     .                         +cnv(ip,jm,k,1))
c$$$     .            +gsub23*0.25*(cnv(i,j ,k ,3)
c$$$     .                         +cnv(i,j ,kp,3)
c$$$     .                         +cnv(i,jm,k ,3)
c$$$     .                         +cnv(i,jm,kp,3))
     .                +0.25*(gsub_v(i,j ,k,2,1)*(cnv(i ,j ,k,1)
     .                                          +cnv(ip,j ,k,1))
     .                      +gsub_v(i,jm,k,2,1)*(cnv(i ,jm,k,1)
     .                                          +cnv(ip,jm,k,1)))
     .                +0.25*(gsub_v(i,j ,k,2,3)*(cnv(i,j ,k ,3)
     .                                          +cnv(i,j ,kp,3))
     .                      +gsub_v(i,jm,k,2,3)*(cnv(i,jm,k ,3)
     .                                          +cnv(i,jm,kp,3)))
        endif

        if (km >= 0) then
          !i+1/2,j+1/2,k
c$$$          gsub31 = 0.5*(gsub_v(i,j,km,3,1)
c$$$     .                 +gsub_v(i,j,k ,3,1))
c$$$          gsub32 = 0.5*(gsub_v(i,j,km,3,2)
c$$$     .                 +gsub_v(i,j,k ,3,2))
          gsub33 = 0.5*(gsub_v(i,j,km,3,3)
     .                 +gsub_v(i,j,k ,3,3))

          cov(3) = gsub33      *cnv(i ,j,k ,3)
c$$$     .            +gsub31*0.25*(cnv(i ,j,k ,1)
c$$$     .                         +cnv(ip,j,k ,1)
c$$$     .                         +cnv(i ,j,km,1)
c$$$     .                         +cnv(ip,j,km,1))
c$$$     .            +gsub32*0.25*(cnv(i,j ,k ,2)
c$$$     .                         +cnv(i,jp,k ,2)
c$$$     .                         +cnv(i,j ,km,2)
c$$$     .                         +cnv(i,jp,km,2))
     .                +0.25*(gsub_v(i,j,k ,3,1)*(cnv(i ,j,k ,1)
     .                                          +cnv(ip,j,k ,1))
     .                      +gsub_v(i,j,km,3,1)*(cnv(i ,j,km,1)
     .                                          +cnv(ip,j,km,1)))
     .                +0.25*(gsub_v(i,j,k ,3,2)*(cnv(i,j ,k ,2)
     .                                          +cnv(i,jp,k ,2))
     .                      +gsub_v(i,j,km,3,2)*(cnv(i,j ,km,2)
     .                                          +cnv(i,jp,km,2)))
        endif

      end function XformToCov_vrtx_stg_ijk

c     Xform_From_VrtxCnv_To_CellCov_stg_ijk
c     #################################################################
      function Xform_From_VrtxCnv_To_CellCov_stg_ijk(g_def,i,j,k,igr
     .                  ,cnv_vrtx) result(cov_cell)

c     -----------------------------------------------------------------
c     Transforms a staggered VERTEX contravariant vector to a staggered
c     CELL covariant representation for cell (i,j,k).
c
c     CURRENTLY WORKS ONLY IN 2D
c     -----------------------------------------------------------------

        implicit none

c     Input variables

        type(grid_mg_def),pointer :: g_def

        integer :: i,j,k,igr
        real(8) :: cnv_vrtx(0:,0:,0:,:),cov_cell(3)

c     Local variables

        integer :: ip,im,jp,jm,kp,km
        real(8),pointer,dimension(:,:,:,:,:) :: gsub_v

c     Begin program

        ip = min(i+1,size(cnv_vrtx,1)-1)
        im = max(i-1,0)

        jp = min(j+1,size(cnv_vrtx,2)-1)
        jm = max(j-1,0)

        kp = min(k+1,size(cnv_vrtx,3)-1)
        km = max(k-1,0)

        gsub_v => g_def%gmetric%grid(igr)%gsub_v

        !i+1/2,j,k
        cov_cell(1) = 0.125*(gsub_v(i,j ,k ,1,1)*(cnv_vrtx(i ,j ,k ,1)
     .                                           +cnv_vrtx(ip,j ,k ,1))
     .                      +gsub_v(i,jm,k ,1,1)*(cnv_vrtx(i ,jm,k ,1)
     .                                           +cnv_vrtx(ip,jm,k ,1))
     .                      +gsub_v(i,j ,km,1,1)*(cnv_vrtx(i ,j ,km,1)
     .                                           +cnv_vrtx(ip,j ,km,1))
     .                      +gsub_v(i,jm,km,1,1)*(cnv_vrtx(i ,jm,km,1)
     .                                           +cnv_vrtx(ip,jm,km,1)))
     .               +0.25*((gsub_v(i,jm,k ,1,2)
     .                      +gsub_v(i,j ,k ,1,2))*cnv_vrtx(i,j,k ,2)
     .                     +(gsub_v(i,jm,km,1,2)
     .                      +gsub_v(i,j ,km,1,2))*cnv_vrtx(i,j,km,2))
     .               +0.25*((gsub_v(i,j ,km,1,3)
     .                      +gsub_v(i,j ,k ,1,3))*cnv_vrtx(i,j ,k,3)
     .                     +(gsub_v(i,jm,km,1,3)               
     .                      +gsub_v(i,jm,k ,1,3))*cnv_vrtx(i,jm,k,3))

        !i,j+1/2,k
        cov_cell(2) = 0.25*((gsub_v(im,j,k ,2,1)
     .                      +gsub_v(i ,j,k ,2,1))*cnv_vrtx(i,j,k ,1)
     .                     +(gsub_v(im,j,km,2,1)
     .                      +gsub_v(i ,j,km,2,1))*cnv_vrtx(i,j,km,1))
     .               +0.125*(gsub_v(i ,j,k ,2,2)*(cnv_vrtx(i ,j ,k ,2)
     .                                           +cnv_vrtx(i ,jp,k ,2))
     .                      +gsub_v(im,j,k ,2,2)*(cnv_vrtx(im,j ,k ,2)
     .                                           +cnv_vrtx(im,jp,k ,2))
     .                      +gsub_v(i ,j,km,2,2)*(cnv_vrtx(i ,j ,km,2)
     .                                           +cnv_vrtx(i ,jp,km,2))
     .                      +gsub_v(im,j,km,2,2)*(cnv_vrtx(im,j ,km,2)
     .                                           +cnv_vrtx(im,jp,km,2)))
     .               +0.25*((gsub_v(i ,j,km,2,3)
     .                      +gsub_v(i ,j,k ,2,3))*cnv_vrtx(i ,j,k,3)
     .                     +(gsub_v(im,j,km,2,3)               
     .                      +gsub_v(im,j,k ,2,3))*cnv_vrtx(im,j,k,3))

        !i,j,k+1/2
        cov_cell(3) = 0.25*((gsub_v(im,j ,k,3,1)
     .                      +gsub_v(i ,j ,k,3,1))*cnv_vrtx(i,j ,k,1)
     .                     +(gsub_v(im,jm,k,3,1)               
     .                      +gsub_v(i ,jm,k,3,1))*cnv_vrtx(i,jm,k,1))
     .               +0.25*((gsub_v(i ,jm,k,3,2)
     .                      +gsub_v(i ,j ,k,3,2))*cnv_vrtx(i ,j,k,2)
     .                     +(gsub_v(im,jm,k,3,2)             
     .                      +gsub_v(im,j ,k,3,2))*cnv_vrtx(im,j,k,2))
     .               +0.125*(gsub_v(i ,j ,k,3,3)*(cnv_vrtx(i ,j ,k ,3)
     .                                           +cnv_vrtx(i ,j ,kp,3))
     .                      +gsub_v(im,j ,k,3,3)*(cnv_vrtx(im,j ,k ,3)
     .                                           +cnv_vrtx(im,j ,kp,3))
     .                      +gsub_v(i ,jm,k,3,3)*(cnv_vrtx(i ,jm,k ,3)
     .                                           +cnv_vrtx(i ,jm,kp,3))
     .                      +gsub_v(im,jm,k,3,3)*(cnv_vrtx(im,jm,k ,3)
     .                                           +cnv_vrtx(im,jm,kp,3)))

      end function Xform_From_VrtxCnv_To_CellCov_stg_ijk

c     XformToCnv_stg_mesh
c     #################################################################
      function XformToCnv_stg_mesh(g_def,igrid,cov) result(cnv)

c     -----------------------------------------------------------------
c     Transforms a staggered covariant vector to a contravariant
c     representation at faces everywhere on the mesh.
c     -----------------------------------------------------------------

        implicit none

c     Input variables

        type(grid_mg_def),pointer :: g_def

        integer :: igrid
        real(8) :: cov(0:,0:,0:,:)
     .            ,cnv(0:size(cov,1)-1
     .                ,0:size(cov,2)-1
     .                ,0:size(cov,3)-1,3)

c     Local variables

        integer :: i,j,k,nx,ny,nz

c     Begin program

        nx = g_def%nxv(igrid)
        ny = g_def%nyv(igrid)
        nz = g_def%nzv(igrid)

        cnv = 0d0

        do k=0,nz+1
          do j=0,ny+1
            do i=0,nx+1
              cnv(i,j,k,:) = XformToCnv_stg_ijk(g_def,i,j,k,igrid,cov)
            enddo
          enddo
        enddo

      end function XformToCnv_stg_mesh

c     XformToCov_stg_mesh
c     #################################################################
      function XformToCov_stg_mesh(g_def,igrid,cnv,vertex_in,vertex_out)
     .         result(cov)

c     -----------------------------------------------------------------
c     Transforms a staggered contravariant vector (vertex-based if
c     vertex=.true.) to a cell-based covariant staggered representation
c     everywhere on the mesh.
c     -----------------------------------------------------------------

        implicit none

c     Input variables

        type(grid_mg_def),pointer :: g_def

        integer :: igrid
        real(8),INTENT(IN) :: cnv(0:,0:,0:,:)
        real(8) :: cov(0:size(cnv,1)-1
     .                ,0:size(cnv,2)-1
     .                ,0:size(cnv,3)-1,3)
        logical, optional :: vertex_in,vertex_out

c     Local variables

        integer :: i,j,k,nx,ny,nz
        logical :: vrtx_in,vrtx_out

c     Begin program

        if (PRESENT(vertex_in)) then
          vrtx_in = vertex_in
        else
          vrtx_in = .false.
        endif

        if (PRESENT(vertex_out)) then
          vrtx_out = vertex_out
        else
          vrtx_out = .false.
        endif

        nx = size(cnv,1)-2
        ny = size(cnv,2)-2
        nz = size(cnv,3)-2

        cov = 0d0

        if (vrtx_in.and.(.not.vrtx_out)) then
          do k=0,nz+1
            do j=0,ny+1
              do i=0,nx+1
                cov(i,j,k,:) =
     .            Xform_From_VrtxCnv_To_CellCov_stg_ijk(g_def,i,j,k
     .                                                 ,igrid,cnv)
              enddo
            enddo
          enddo
        elseif ((.not.vrtx_in).and.(.not.vrtx_out)) then
          do k=0,nz+1
            do j=0,ny+1
              do i=0,nx+1
                cov(i,j,k,:) = XformToCov_stg_ijk(g_def,i,j,k,igrid,cnv)
              enddo
            enddo
          enddo
        elseif (vrtx_in.and.vrtx_out) then
          do k=0,nz+1
            do j=0,ny+1
              do i=0,nx+1
                cov(i,j,k,:) = XformToCov_vrtx_stg_ijk(g_def,i,j,k,igrid
     .                                                ,cnv)
              enddo
            enddo
          enddo
        elseif ((.not.vrtx_in).and.vrtx_out) then
          write (*,*) 'Case not implemented yet'
          write (*,*) 'Aborting in XformToCov_stg_mesh'
          stop
        endif

      end function XformToCov_stg_mesh

ccc     XformVector_stg_BC
ccc     ###################################################################
cc      subroutine XformVector_stg_BC(igrid,vec1,vec2,tocnv)
cc
cc      implicit none
cc
ccc     -------------------------------------------------------------------
ccc     Routine to transfer ghost cell values form vec1 to vec2
ccc     -------------------------------------------------------------------
cc
ccc     Input variables
cc
cc      integer :: igrid
cc      real(8) :: vec1(0:,0:,0:,:)
cc      real(8) :: vec2(0:,0:,0:,:)
cc      logical :: tocnv
cc
ccc     Local variables
cc
cc      integer :: nxp,nyp,nzp,i,j,k
cc
ccc     Begin program
cc
cc      nxp = size(vec1,1)-1
cc      nyp = size(vec1,2)-1
cc      nzp = size(vec1,3)-1
cc
cc      if (tocnv) then
cc
cc        !X boundaries
cc        do k= 0,nzp
cc          do j= 0,nyp
cc            vec2(0  ,j,k,:) = XformToCnv_stg_ijk(0  ,j,k,igrid,vec1)
cc            vec2(nxp,j,k,:) = XformToCnv_stg_ijk(nxp,j,k,igrid,vec1)
cc          enddo
cc        enddo
cc
cc        !Y boundaries
cc        do k= 0,nzp
cc          do i=0,nxp
cc            vec2(i,0  ,k,:) = XformToCnv_stg_ijk(i,0  ,k,igrid,vec1)
cc            vec2(i,nyp,k,:) = XformToCnv_stg_ijk(i,nyp,k,igrid,vec1)
cc          enddo
cc        enddo
cc
cc        !Z boundaries
cc        do j= 0,nyp
cc          do i= 0,nxp
cc            vec2(i,j,0  ,:) = XformToCnv_stg_ijk(i,j,0  ,igrid,vec1)
cc            vec2(i,j,nzp,:) = XformToCnv_stg_ijk(i,j,nzp,igrid,vec1)
cc          enddo
cc        enddo
cc
cc      else
cc        call pstop("XformVector_stg_BC","Option not implemented")
cc      endif
cc
ccc     End program
cc
cc      end subroutine XformVector_stg_BC

c     XformVector_stg_mesh
c     #################################################################
      function XformVector_stg_mesh(g_def,igrid,vec1,tocnv) result(vec2)
c     -----------------------------------------------------------------
c     Transforms a STAGGERED curvilinear vector from covariant to
c     contravariant (tocnv=.true.) or viceversa
c     -----------------------------------------------------------------

      implicit none

c     Input variables

      type(grid_mg_def),pointer :: g_def

      integer    :: igrid
      real(8)    :: vec1(0:,0:,0:,:)
      real(8)    :: vec2(0:size(vec1,1)-1
     .                  ,0:size(vec1,2)-1
     .                  ,0:size(vec1,3)-1
     .                  ,size(vec1,4))

      logical    :: tocnv

c     Local variables

c     Begin program

      if (tocnv) then
        vec2 = XformToCnv_stg_mesh(g_def,igrid,vec1)
      else
        vec2 = XformToCov_stg_mesh(g_def,igrid,vec1)
      endif

c     End program

      end function XformVector_stg_mesh

c     XformVector_stg_mesh2
c     #################################################################
      subroutine XformVector_stg_mesh2(g_def,igrid,vec,in_repr,out_repr)
c     -----------------------------------------------------------------
c     Transforms a STAGGERED vector from one representation to another.
c     -----------------------------------------------------------------

      implicit none

c     Input variables

      type(grid_mg_def),pointer :: g_def

      integer    :: igrid
      real(8)    :: vec(0:,0:,0:,:)

      character(3) :: in_repr,out_repr

c     Local variables

      integer :: imin,imax,jmin,jmax,kmin,kmax

      logical :: covariant,cnvariant,cartesian

c     Begin program

      cartesian = (in_repr == 'car')
      covariant = (in_repr == 'cov')
      cnvariant = (.not. cartesian) .and. (.not. covariant)

      select case(out_repr)
      case('car')

        if (cartesian) return

        vec = XformToCar_stg_mesh(g_def,igrid,vec,covariant)

      case('cnv')

        if (cnvariant) return

        if (cartesian) then
          vec = XformToCurv_stg_mesh(g_def,igrid,vec,.true.)
        else
          vec = XformToCnv_stg_mesh (g_def,igrid,vec)
        endif

      case('cov')

        if (covariant) return

        if (cartesian) then
          vec = XformToCurv_stg_mesh(g_def,igrid,vec,.false.)
        else
          vec = XformToCov_stg_mesh (g_def,igrid,vec)
        endif

      end select

c     End program

      end subroutine XformVector_stg_mesh2

c     Xfer_stg_to_cell_mesh
c     #################################################################
      function Xfer_stg_to_cell_mesh(g_def,stg,vrtx) result(cell)

c     -----------------------------------------------------------------
c     Averages staggered vector to cell centers. The variable vrtx
c     determines whether vector is cell-based (components orthogonal to
c     cell-center volumes) or vertex-based (components orthogonal to
c     vertex-centered volumes)
c     -----------------------------------------------------------------

        implicit none

c     Input variables

        type(grid_mg_def),pointer :: g_def

        integer :: igrid
        real(8) :: stg(0:,0:,0:,:)
     .            ,cell(0:size(stg,1)-1
     .                 ,0:size(stg,2)-1
     .                 ,0:size(stg,3)-1
     .                 ,  size(stg,4)  )

        logical,optional :: vrtx

c     Local variables

        integer :: i,j,k,nx,ny,nz,neq
        logical :: vertex

c     Begin program

        if (PRESENT(vrtx)) then
          vertex = vrtx
        else
          vertex = .false.
        endif

        nx = size(stg,1)-2
        ny = size(stg,2)-2
        nz = size(stg,3)-2
        neq= size(stg,4)

        cell = 0d0

        if (.not.vertex) then
          cell(1:nx,:,:,1) = 0.5*(stg(0:nx-1,:,:,1)+stg(1:nx,:,:,1))
          cell(:,1:ny,:,2) = 0.5*(stg(:,0:ny-1,:,2)+stg(:,1:ny,:,2))
          cell(:,:,1:nz,3) = 0.5*(stg(:,:,0:nz-1,3)+stg(:,:,1:nz,3))
        else
          cell(1:nx,1:ny,1:nz,1) = 0.25*(stg(1:nx,0:ny-1,0:nz-1,1)
     .                                  +stg(1:nx,0:ny-1,1:nz  ,1)
     .                                  +stg(1:nx,1:ny  ,1:nz  ,1)
     .                                  +stg(1:nx,1:ny  ,0:nz-1,1))
          cell(1:nx,1:ny,1:nz,2) = 0.25*(stg(0:nx-1,1:ny,0:nz-1,2)
     .                                  +stg(0:nx-1,1:ny,1:nz  ,2)
     .                                  +stg(1:nx  ,1:ny,1:nz  ,2)
     .                                  +stg(1:nx  ,1:ny,0:nz-1,2))
          cell(1:nx,1:ny,1:nz,3) = 0.25*(stg(0:nx-1,0:ny-1,1:nz,3)
     .                                  +stg(0:nx-1,1:ny  ,1:nz,3)
     .                                  +stg(1:nx  ,1:ny  ,1:nz,3)
     .                                  +stg(1:nx  ,0:ny-1,1:nz,3))
        endif

      end function Xfer_stg_to_cell_mesh

c     crossProduct_stg_mesh
c     ###############################################################
      function crossProduct_stg_mesh(g_def,igrid,vec1,vertex1
     .                              ,vec2,vertex2,cnv_in)
     .         result(vec3)
c     ---------------------------------------------------------------
c     Calculates cross product of staggered vectors, and returns a
c     staggered vector at cell faces. Staggered input vectors can be
c     cell-based (vertex=.false.) or vertex-based (vertex=.true.). Also,
c     they can be contravariant (cnv_in=.true.) or covariant
c     (cnv_in=.false.). If the former, they return a covariant vector,
c     and contravariant otherwise.
c     ---------------------------------------------------------------

      implicit none           !For safe fortran

c     Call variables

      type(grid_mg_def),pointer :: g_def

      integer :: igrid

      real(8) :: vec1(0:,0:,0:,:)
     .          ,vec2(0:,0:,0:,:)
     .          ,vec3(0:size(vec1,1)-1
     .               ,0:size(vec1,2)-1
     .               ,0:size(vec1,3)-1
     .               ,  size(vec1,4))

      logical :: vertex1,vertex2,cnv_in

c     Local variables

      integer :: nx,ny,nz

c     Begin program

      nx = size(vec1,1)-2
      ny = size(vec1,2)-2
      nz = size(vec1,3)-2

      vec3 = 0d0

      if (.not.vertex1.and.vertex2) then
        !At (i+1/2,j,k)
        vec3(1:nx,1:ny,1:nz,1) =(0.25*(vec1(1:nx  ,1:ny  ,1:nz,2)
     .                                +vec1(2:nx+1,1:ny  ,1:nz,2)
     .                                +vec1(1:nx  ,0:ny-1,1:nz,2)
     .                                +vec1(2:nx+1,0:ny-1,1:nz,2))
     .                          * 0.5*(vec2(1:nx,1:ny  ,1:nz,3)
     .                                +vec2(1:nx,0:ny-1,1:nz,3))
     .                          -0.25*(vec1(1:nx  ,1:ny,1:nz  ,3)
     .                                +vec1(2:nx+1,1:ny,1:nz  ,3)
     .                                +vec1(1:nx  ,1:ny,0:nz-1,3)
     .                                +vec1(2:nx+1,1:ny,0:nz-1,3))
     .                          * 0.5*(vec2(1:nx  ,1:ny,1:nz,2)
     .                                +vec2(0:nx-1,1:ny,1:nz,2)))

        !At (i,j+1/2,k)
        vec3(1:nx,1:ny,1:nz,2) =(0.25*(vec1(1:nx,1:ny  ,1:nz  ,3)
     .                                +vec1(1:nx,2:ny+1,1:nz  ,3)
     .                                +vec1(1:nx,1:ny  ,0:nz-1,3)
     .                                +vec1(1:nx,2:ny+1,0:nz-1,3))
     .                           *0.5*(vec2(1:nx,1:ny  ,1:nz  ,1)
     .                                +vec2(1:nx,1:ny  ,0:nz-1,1))
     .                          -0.25*(vec1(1:nx  ,1:ny  ,1:nz,1)
     .                                +vec1(1:nx  ,2:ny+1,1:nz,1)
     .                                +vec1(0:nx-1,1:ny  ,1:nz,1)
     .                                +vec1(0:nx-1,2:ny+1,1:nz,1))
     .                           *0.5*(vec2(1:nx  ,1:ny,1:nz,3)
     .                                +vec2(0:nx-1,1:ny,1:nz,3)))

        !At (i,j,k+1/2)
        vec3(1:nx,1:ny,1:nz,3) =(0.25*(vec1(1:nx  ,1:ny,1:nz  ,1)
     .                                +vec1(1:nx  ,1:ny,1:nz  ,1)
     .                                +vec1(0:nx-1,1:ny,2:nz+1,1)
     .                                +vec1(0:nx-1,1:ny,2:nz+1,1))
     .                           *0.5*(vec2(1:nx  ,1:ny,1:nz,2)
     .                                +vec2(0:nx-1,1:ny,1:nz  ,2))
     .                          -0.25*(vec1(1:nx,1:ny  ,1:nz  ,2)
     .                                +vec1(1:nx,0:ny-1,1:nz  ,2)
     .                                +vec1(1:nx,1:ny  ,2:nz+1,2)
     .                                +vec1(1:nx,0:ny-1,2:nz+1,2))
     .                           *0.5*(vec2(1:nx,1:ny  ,1:nz,1)
     .                                +vec2(1:nx,0:ny-1,1:nz,1)))

        !Apply Jacobian factor
        if (cnv_in) then
          vec3(1:nx,1:ny,1:nz,1) = 4*vec3(1:nx,1:ny,1:nz,1)
     .            /(g_def%gmetric%grid(igrid)%jac_v(1:nx,1:ny  ,1:nz)
     .             +g_def%gmetric%grid(igrid)%jac_v(1:nx,0:ny-1,1:nz)
     .             +g_def%gmetric%grid(igrid)%jac_v(1:nx,1:ny  ,0:nz-1)
     .             +g_def%gmetric%grid(igrid)%jac_v(1:nx,0:ny-1,0:nz-1))
          vec3(1:nx,1:ny,1:nz,2) = 4*vec3(1:nx,1:ny,1:nz,2)
     .            /(g_def%gmetric%grid(igrid)%jac_v(1:nx  ,1:ny,1:nz)
     .             +g_def%gmetric%grid(igrid)%jac_v(0:nx-1,1:ny,1:nz)
     .             +g_def%gmetric%grid(igrid)%jac_v(1:nx  ,1:ny,0:nz-1)
     .             +g_def%gmetric%grid(igrid)%jac_v(0:nx-1,1:ny,0:nz-1))
          vec3(1:nx,1:ny,1:nz,3) = 4*vec3(1:nx,1:ny,1:nz,3)
     .            /(g_def%gmetric%grid(igrid)%jac_v(1:nx  ,1:ny  ,1:nz)
     .             +g_def%gmetric%grid(igrid)%jac_v(0:nx-1,1:ny  ,1:nz)
     .             +g_def%gmetric%grid(igrid)%jac_v(1:nx  ,0:ny-1,1:nz)
     .             +g_def%gmetric%grid(igrid)%jac_v(0:nx-1,0:ny-1,1:nz))
        endif

      else
        write (*,*) 'Case not implemented yet'
        write (*,*) 'Aborting in crossProduct_stg_mesh'
        stop
c$$$        call pstop("crossProduct_stg_mesh"
c$$$     .            ,"Case not implemented yet")
      endif

c     End program

      end function crossProduct_stg_mesh

c     scalarProduct_stg_mesh
c     ###############################################################
      function scalarProduct_stg_mesh(g_def,igrid,vec1,vertex1
     .                                           ,vec2,vertex2)
     .         result(scl)
c     ---------------------------------------------------------------
c     Calculates dot product of staggered vectors, and returns a
c     three scalar values at cell faces. Staggered input vectors can be
c     cell-based (vertex=.false.) or vertex-based (vertex=.true.).
c     ---------------------------------------------------------------

      implicit none           !For safe fortran

c     Call variables

      type(grid_mg_def),pointer :: g_def

      integer :: igrid

      real(8) :: vec1(0:,0:,0:,:)
     .          ,vec2(0:,0:,0:,:)
     .          ,scl (0:size(vec1,1)-1
     .               ,0:size(vec1,2)-1
     .               ,0:size(vec1,3)-1
     .               ,3)

      logical :: vertex1,vertex2

c     Local variables

      integer :: nx,ny,nz
      real(8) :: B1(1:size(vec1,1)-2
     .             ,1:size(vec1,2)-2
     .             ,1:size(vec1,3)-2)
     .          ,B2(1:size(vec1,1)-2
     .             ,1:size(vec1,2)-2
     .             ,1:size(vec1,3)-2)

      type(MG_grid_metrics),pointer :: gmetric 

c     Begin program

      gmetric => g_def%gmetric

      nx = size(vec1,1)-2
      ny = size(vec1,2)-2
      nz = size(vec1,3)-2

      scl = 0d0

      if (.not.vertex1.and.vertex2) then
!     At (i+1/2,j,k)
         scl(1:nx,1:ny,1:nz,1) =0.125d0*(
     .                    vec1(1:nx  ,1:ny  ,1:nz  ,1)
     .                  *(vec2(1:nx  ,1:ny  ,1:nz  ,1)+
     .                    vec2(1:nx  ,0:ny-1,1:nz  ,1)+         
     .                    vec2(2:nx+1,0:ny-1,1:nz  ,1)+         
     .                    vec2(2:nx+1,1:ny  ,1:nz  ,1)+
     .                    vec2(1:nx  ,1:ny  ,0:nz-1,1)+
     .                    vec2(1:nx  ,0:ny-1,0:nz-1,1)+         
     .                    vec2(2:nx+1,0:ny-1,0:nz-1,1)+         
     .                    vec2(2:nx+1,1:ny  ,0:nz-1,1))
     .                +  (vec1(1:nx  ,1:ny  ,1:nz,2)+
     .                    vec1(1:nx  ,0:ny-1,1:nz,2)+         
     .                    vec1(2:nx+1,0:ny-1,1:nz,2)+         
     .                    vec1(2:nx+1,1:ny  ,1:nz,2)) 
     .                  *(vec2(1:nx  ,1:ny  ,1:nz  ,2)
     .                   +vec2(1:nx  ,1:ny  ,0:nz-1,2))
     .                +  (vec1(1:nx  ,1:ny  ,1:nz  ,3)+
     .                    vec1(2:nx+1,1:ny  ,1:nz  ,3)+
     .                    vec1(1:nx  ,1:ny  ,0:nz-1,3)+
     .                    vec1(2:nx+1,1:ny  ,0:nz-1,3))
     .                  *(vec2(1:nx  ,1:ny  ,1:nz  ,3)+
     .                    vec2(1:nx  ,0:ny-1,1:nz  ,3)) )

!     At (i,j+1/2,k)
         scl(1:nx,1:ny,1:nz,2) = 0.125d0*(
     .                   (vec1(1:nx  ,1:ny  ,1:nz,1)+
     .                    vec1(1:nx  ,2:ny+1,1:nz,1)+
     .                    vec1(0:nx-1,2:ny+1,1:nz,1)+
     .                    vec1(0:nx-1,1:ny  ,1:nz,1))
     .                  *(vec2(1:nx  ,1:ny  ,1:nz  ,1)+
     .                    vec2(1:nx  ,1:ny  ,0:nz-1,1))
     .                 +  vec1(1:nx  ,1:ny  ,1:nz,2)
     .                  *(vec2(1:nx  ,1:ny  ,1:nz  ,2)+
     .                    vec2(1:nx  ,2:ny+1,1:nz  ,2)+ 
     .                    vec2(0:nx-1,2:ny+1,1:nz  ,2)+ 
     .                    vec2(0:nx-1,1:ny  ,1:nz  ,2)+
     .                    vec2(1:nx  ,1:ny  ,0:nz-1,2)+
     .                    vec2(1:nx  ,2:ny+1,0:nz-1,2)+ 
     .                    vec2(0:nx-1,2:ny+1,0:nz-1,2)+ 
     .                    vec2(0:nx-1,1:ny  ,0:nz-1,2))
     .                 + (vec1(1:nx  ,1:ny  ,1:nz  ,3)+
     .                    vec1(1:nx  ,2:ny+1,1:nz  ,3)+
     .                    vec1(1:nx  ,1:ny  ,0:nz-1,3)+
     .                    vec1(1:nx  ,2:ny+1,0:nz-1,3))
     .                  *(vec2(1:nx  ,1:ny  ,1:nz  ,3)+
     .                    vec2(0:nx-1,1:ny  ,1:nz  ,3)) )

!     At (i,j,k+1/2)
         scl(1:nx,1:ny,1:nz,3) = 0.125d0*(
     .                   (vec1(1:nx  ,1:ny  ,1:nz  ,1)+
     .                    vec1(0:nx-1,1:ny  ,1:nz  ,1)+
     .                    vec1(1:nx  ,1:ny  ,2:nz+1,1)+
     .                    vec1(0:nx-1,1:ny  ,2:nz+1,1))
     .                  *(vec2(1:nx  ,1:ny  ,1:nz,1)+
     .                    vec2(1:nx  ,0:ny-1,1:nz,1))
     .                 + (vec1(1:nx  ,1:ny  ,1:nz  ,2)+
     .                    vec1(1:nx  ,0:ny-1,1:nz  ,2)+
     .                    vec1(1:nx  ,1:ny  ,2:nz+1,2)+
     .                    vec1(1:nx  ,0:ny-1,2:nz+1,2))
     .                  *(vec2(1:nx  ,1:ny  ,1:nz,2)+
     .                    vec2(0:nx-1,1:ny  ,1:nz,2)) 
     .                 +  vec1(1:nx  ,1:ny  ,1:nz,3)
     .                  *(vec2(1:nx  ,1:ny  ,1:nz  ,3)+
     .                    vec2(0:nx-1,1:ny  ,1:nz  ,3)+ 
     .                    vec2(0:nx-1,0:ny-1,1:nz  ,3)+ 
     .                    vec2(1:nx  ,0:ny-1,1:nz  ,3)+
     .                    vec2(1:nx  ,1:ny  ,2:nz+1,3)+
     .                    vec2(0:nx-1,1:ny  ,2:nz+1,3)+ 
     .                    vec2(0:nx-1,0:ny-1,2:nz+1,3)+ 
     .                    vec2(1:nx  ,0:ny-1,2:nz+1,3)) )

      elseif (vertex1.and.vertex2) then
!     At (i+1/2,j,k)
!     Bx
         B1 = (vec1(1:nx  ,1:ny  ,1:nz  ,1)+
     .         vec1(2:nx+1,1:ny  ,1:nz  ,1)+
     .         vec1(2:nx+1,0:ny-1,1:nz  ,1)+
     .         vec1(1:nx  ,0:ny-1,1:nz  ,1)+
     .         vec1(1:nx  ,1:ny  ,0:nz-1,1)+
     .         vec1(2:nx+1,1:ny  ,0:nz-1,1)+
     .         vec1(2:nx+1,0:ny-1,0:nz-1,1)+
     .         vec1(1:nx  ,0:ny-1,0:nz-1,1))*0.125d0
         B2 = (vec2(1:nx  ,1:ny  ,1:nz,  1)+       
     .         vec2(2:nx+1,1:ny  ,1:nz  ,1)+       
     .         vec2(2:nx+1,0:ny-1,1:nz  ,1)+       
     .         vec2(1:nx  ,0:ny-1,1:nz  ,1)+
     .         vec2(1:nx  ,1:ny  ,0:nz-1,1)+       
     .         vec2(2:nx+1,1:ny  ,0:nz-1,1)+       
     .         vec2(2:nx+1,0:ny-1,0:nz-1,1)+       
     .         vec2(1:nx  ,0:ny-1,0:nz-1,1))*0.125d0 
         scl(1:nx,1:ny,1:nz,1) = scl(1:nx,1:ny,1:nz,1) + B1*B2
                                                 
!     By
         B1 = (vec1(1:nx  ,1:ny  ,1:nz  ,2)+
     .         vec1(1:nx  ,1:ny  ,0:nz-1,2))*0.5d0
         B2 = (vec2(1:nx  ,1:ny  ,1:nz  ,2)+
     .         vec2(1:nx  ,1:ny  ,0:nz-1,2))*0.5d0
         scl(1:nx,1:ny,1:nz,1) = scl(1:nx,1:ny,1:nz,1) + B1*B2

!     Bz
         B1 =(vec1(1:nx  ,1:ny  ,1:nz  ,3)+
     .        vec1(1:nx  ,0:ny-1,1:nz  ,3))*0.5d0
         B2 =(vec2(1:nx  ,1:ny  ,1:nz  ,3)+         
     .        vec2(1:nx  ,0:ny-1,1:nz  ,3))*0.5d0
         scl(1:nx,1:ny,1:nz,1) = scl(1:nx,1:ny,1:nz,1) + B1*B2

!     At (i,j+1/2,k)
!     Bx
         B1 =(vec1(1:nx  ,1:ny  ,1:nz  ,1)+
     .        vec1(1:nx  ,1:ny  ,0:nz-1,1))*0.5d0
         B2 =(vec2(1:nx  ,1:ny  ,1:nz  ,1)+
     .        vec2(1:nx  ,1:ny  ,0:nz-1,1))*0.5d0
         scl(1:nx,1:ny,1:nz,2) = scl(1:nx,1:ny,1:nz,2) + B1*B2
                                                 
!     By
         B1 =(vec1(1:nx  ,1:ny  ,1:nz  ,2)+
     .        vec1(1:nx  ,2:ny+1,1:nz  ,2)+ 
     .        vec1(0:nx-1,2:ny+1,1:nz  ,2)+ 
     .        vec1(0:nx-1,1:ny  ,1:nz  ,2)+
     .        vec1(1:nx  ,1:ny  ,0:nz-1,2)+
     .        vec1(1:nx  ,2:ny+1,0:nz-1,2)+ 
     .        vec1(0:nx-1,2:ny+1,0:nz-1,2)+ 
     .        vec1(0:nx-1,1:ny  ,0:nz-1,2))*0.125d0
         B2 =(vec2(1:nx  ,1:ny  ,1:nz  ,2)+
     .        vec2(1:nx  ,2:ny+1,1:nz  ,2)+ 
     .        vec2(0:nx-1,2:ny+1,1:nz  ,2)+ 
     .        vec2(0:nx-1,1:ny  ,1:nz  ,2)+
     .        vec2(1:nx  ,1:ny  ,0:nz-1,2)+
     .        vec2(1:nx  ,2:ny+1,0:nz-1,2)+ 
     .        vec2(0:nx-1,2:ny+1,0:nz-1,2)+ 
     .        vec2(0:nx-1,1:ny  ,0:nz-1,2))*0.125d0
         scl(1:nx,1:ny,1:nz,2) = scl(1:nx,1:ny,1:nz,2) + B1*B2

!     Bz
         B1 =(vec1(1:nx  ,1:ny  ,1:nz  ,3)+
     .        vec1(0:nx-1,1:ny  ,1:nz  ,3))*0.5d0
         B2 =(vec2(1:nx  ,1:ny  ,1:nz  ,3)+
     .        vec2(0:nx-1,1:ny  ,1:nz  ,3))*0.5d0
         scl(1:nx,1:ny,1:nz,2) = scl(1:nx,1:ny,1:nz,2) + B1*B2

!     At (i,j,k+1/2)
!     Bx
         B1 =(vec1(1:nx  ,1:ny  ,1:nz,1)+
     .        vec1(1:nx  ,0:ny-1,1:nz,1))*0.5d0
         B2 =(vec2(1:nx  ,1:ny  ,1:nz,1)+
     .        vec2(1:nx  ,0:ny-1,1:nz,1))*0.5d0
         scl(1:nx,1:ny,1:nz,3) = scl(1:nx,1:ny,1:nz,3) + B1*B2
                                                 
!     By
         B1 =(vec1(1:nx  ,1:ny  ,1:nz,2)+
     .        vec1(0:nx-1,1:ny  ,1:nz,2))*0.5d0
         B2 =(vec2(1:nx  ,1:ny  ,1:nz,2)+
     .        vec2(0:nx-1,1:ny  ,1:nz,2))*0.5d0
         scl(1:nx,1:ny,1:nz,3) = scl(1:nx,1:ny,1:nz,3) + B1*B2

!     Bz
         B1 =(vec1(1:nx  ,1:ny  ,1:nz  ,3)+
     .        vec1(0:nx-1,1:ny  ,1:nz  ,3)+ 
     .        vec1(0:nx-1,0:ny-1,1:nz  ,3)+ 
     .        vec1(1:nx  ,0:ny-1,1:nz  ,3)+
     .        vec1(1:nx  ,1:ny  ,2:nz+1,3)+
     .        vec1(0:nx-1,1:ny  ,2:nz+1,3)+ 
     .        vec1(0:nx-1,0:ny-1,2:nz+1,3)+ 
     .        vec1(1:nx  ,0:ny-1,2:nz+1,3))*0.125d0
         B2 =(vec2(1:nx  ,1:ny  ,1:nz  ,3)+
     .        vec2(0:nx-1,1:ny  ,1:nz  ,3)+ 
     .        vec2(0:nx-1,0:ny-1,1:nz  ,3)+ 
     .        vec2(1:nx  ,0:ny-1,1:nz  ,3)+
     .        vec2(1:nx  ,1:ny  ,2:nz+1,3)+
     .        vec2(0:nx-1,1:ny  ,2:nz+1,3)+ 
     .        vec2(0:nx-1,0:ny-1,2:nz+1,3)+ 
     .        vec2(1:nx  ,0:ny-1,2:nz+1,3))*0.125d0
         scl(1:nx,1:ny,1:nz,3) = scl(1:nx,1:ny,1:nz,3) + B1*B2
      else
        write (*,*) 'Case not implemented yet'
        write (*,*) 'Aborting in scalarProduct_stg_mesh'
        stop
c$$$        call pstop("scalarProduct_stg_mesh"
c$$$     .            ,"Case not implemented yet")
      endif                                       

      !Apply Jacobian factor
      scl(1:nx,1:ny,1:nz,1) = 4*scl(1:nx,1:ny,1:nz,1)
     .            /(g_def%gmetric%grid(igrid)%jac_v(1:nx,1:ny  ,1:nz)
     .             +g_def%gmetric%grid(igrid)%jac_v(1:nx,0:ny-1,1:nz)
     .             +g_def%gmetric%grid(igrid)%jac_v(1:nx,1:ny  ,0:nz-1)
     .             +g_def%gmetric%grid(igrid)%jac_v(1:nx,0:ny-1,0:nz-1))
      scl(1:nx,1:ny,1:nz,2) = 4*scl(1:nx,1:ny,1:nz,2)
     .            /(g_def%gmetric%grid(igrid)%jac_v(1:nx  ,1:ny,1:nz)
     .             +g_def%gmetric%grid(igrid)%jac_v(0:nx-1,1:ny,1:nz)
     .             +g_def%gmetric%grid(igrid)%jac_v(1:nx  ,1:ny,0:nz-1)
     .             +g_def%gmetric%grid(igrid)%jac_v(0:nx-1,1:ny,0:nz-1))
      scl(1:nx,1:ny,1:nz,3) = 4*scl(1:nx,1:ny,1:nz,3)
     .            /(g_def%gmetric%grid(igrid)%jac_v(1:nx  ,1:ny  ,1:nz)
     .             +g_def%gmetric%grid(igrid)%jac_v(0:nx-1,1:ny  ,1:nz)
     .             +g_def%gmetric%grid(igrid)%jac_v(1:nx  ,0:ny-1,1:nz)
     .             +g_def%gmetric%grid(igrid)%jac_v(0:nx-1,0:ny-1,1:nz))

c     End program

      end function scalarProduct_stg_mesh

      end module grid_operations_stg

c module grid_operations
c #####################################################################
      module grid_operations

      use grid_operations_cell

      use grid_operations_stg

      end module grid_operations
    
