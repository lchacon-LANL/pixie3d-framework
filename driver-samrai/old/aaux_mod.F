c module equilibrium
c ######################################################################
      module equilibrium

        integer    :: IRHO,IVX,IVY,IVZ,IBX,IBY,IBZ,ITMP,IJX,IJY,IJZ
     .               ,IAX,IAY,IAZ,IBCNV,IBCOV,IJCNV,IJCOV

#if !defined(vec_pot)
        parameter    (IRHO=1,IVX=2,IVY=3,IVZ=4,IBX=5,IBY=6,IBZ=7,ITMP=8
     .               ,IJX=9,IJY=10,IJZ=11,IAX=12,IAY=13,IAZ=14
     .               ,IJCNV=1,IJCOV=2,IBCNV=3,IBCOV=4)
#else
        parameter    (IRHO=1,IVX=2,IVY=3,IVZ=4,IAX=5,IAY=6,IAZ=7,ITMP=8
     .               ,IBX=9,IBY=10,IBZ=11,IJX=12,IJY=13,IJZ=14
     .               ,IJCNV=1,IJCOV=2,IBCNV=3,IBCOV=4)
#endif

#if defined(samrai)
#if !defined(vec_pot)
        real(8),pointer,dimension(:,:,:):: rho,px,py,pz,bx,by,bz,tmp
#else
        real(8),pointer,dimension(:,:,:):: rho,px,py,pz,ax,ay,az,tmp
#endif
#endif
        real(8)    :: dlambda,rshear,vparflow,vperflow

        character(5)  :: equil
        character(20) :: equ_file,prt_file
        real(8)       :: eq_params(6)

      end module equilibrium

c module grid_aliases
c ######################################################################
      module grid_aliases

        use grid

        real(8),pointer,dimension(:) :: xx,yy,zz,dxh,dyh,dzh,dx,dy,dz

        integer    :: igx,igy,igz,nx,ny,nz

        real(8)    :: xim,yim,zim,xip,yip,zip
     .               ,xjm,yjm,zjm,xjp,yjp,zjp
     .               ,xkm,ykm,zkm,xkp,ykp,zkp

        real(8)    :: x0,y0,z0,xh,yh,zh

        real(8)    :: jacip,jacim,jacjp,jacjm,jackp,jackm
     .               ,jacp,jacm,jach,jac0

        real(8)    :: gsub(3,3),gsuper(3,3),jac

        real(8)    :: nabla_v(3,3),hessian(3,3,3)
     .               ,cov_tnsr(3,3),cnv_tnsr(3,3)

        logical    :: cartesian

      end module grid_aliases

c module auxiliaryVariables
c ######################################################################
      module auxiliaryVariables

        use variables

        use equilibrium
#if defined(samrai)
#ifdef defined(vec_pot)

           real(8),pointer,dimension(:,:,:,:) :: acnv
           real(8),pointer,dimension(:,:,:,:) :: acov
           real(8),pointer,dimension(:,:,:,:) :: vlap_a
           real(8),pointer,dimension(:,:,:,:) :: bcnv
           real(8),pointer,dimension(:,:,:,:) :: bcov
#else
           real(8),pointer,dimension(:,:,:,:) :: bcnv
           real(8),pointer,dimension(:,:,:,:) :: bcov
#endif

           real(8),pointer,dimension(:,:,:,:) :: pcnv
           real(8),pointer,dimension(:,:,:,:) :: vecnv
           real(8),pointer,dimension(:,:,:,:) :: vdummy
           real(8),pointer,dimension(:,:,:,:) :: vcnv
           real(8),pointer,dimension(:,:,:,:) :: vcov
           real(8),pointer,dimension(:,:,:,:) :: grd

           real(8),pointer,dimension(:,:,:,:) :: b_n
           real(8),pointer,dimension(:,:,:,:) :: dv_dt

c          scalars           
           real(8),pointer,dimension(:,:,:,:) :: p_n

           real(8),pointer,dimension(:,:,:) :: eeta
           real(8),pointer,dimension(:,:,:) :: nuu
           real(8),pointer,dimension(:,:,:) :: divrgV
#else

#if defined(vec_pot)
        real(8),target,allocatable,dimension(:,:,:,:) :: acov,acnv
     .                                                  ,vlap_a
        real(8),pointer,dimension(:,:,:,:) :: bcnv,bcov
#else
        real(8),target,allocatable,dimension(:,:,:,:) :: bcnv,bcov
#endif
cc        real(8),target,allocatable,dimension(:,:,:,:) :: bcnv,bcov

        real(8),target,allocatable,dimension(:,:,:) :: eeta,nuu,divrgV

        real(8),target,allocatable,dimension(:,:,:,:) ::vcnv,vcov
     .         ,pcnv,vecnv,vdummy,grd

c SI operator
        real(8),target,allocatable,dimension(:,:,:,:) :: b_n,p_n,dv_dt
c SI operator

#endif

        real(8),pointer,dimension(:,:,:,:) :: jcnv,jcov
cc        real(8),target,allocatable,dimension(:,:,:,:) :: jcnv,jcov

        logical :: alt_eom

        logical :: nc_eom_jxb=.false.,nc_eom_gp=.false.,nc_eom_v=.false.


        real(8) :: k_si=0d0

      contains

c     allocAuxVariables
c     ###################################################################
      subroutine allocAuxVariables

c     -------------------------------------------------------------------
c     Allocates auxiliary variable storage.
c     -------------------------------------------------------------------

c     Begin program

c     Allocate storage

        if (.not.associated(gv%aux)) then
#if defined(vec_pot)
          call allocateAuxStruct(0,4,gv%aux)
#else
          call allocateAuxStruct(0,2,gv%aux)
#endif
        endif

        gv%aux%vec_list(IJCNV)%cnv = .true.
        jcnv => gv%aux%vec_list(IJCNV)%vec

        gv%aux%vec_list(IJCOV)%cnv = .false.
        jcov => gv%aux%vec_list(IJCOV)%vec

#if defined(vec_pot)
        gv%aux%vec_list(IBCNV)%cnv = .true.
        bcnv => gv%aux%vec_list(IBCNV)%vec

        gv%aux%vec_list(IBCOV)%cnv = .false.
        bcov => gv%aux%vec_list(IBCOV)%vec
#endif

c     End program

      end subroutine allocAuxVariables

c     defineAuxBCs
c     ###################################################################
cc      subroutine defineAuxBCs(neq,bbcs)
      subroutine defineAuxBCs

c     -------------------------------------------------------------------
c     Defines BCs for auxiliary quantities
c     -------------------------------------------------------------------

c     Call variables

cc      integer    :: neq,bbcs(6,neq)

c     Local variables

      integer    :: ieq,bcsq(6)

c     Begin program

cc#if defined(vec_pot)
cc      !Magnetic field
cc      bcsq = bcond
cc      where (bcsq == DEF) bcsq = DIR
cc      gv%aux%vec_list(IBCNV)%bconds(:,1) = bcsq
cc
cccc      if (equil == 'spnch' .or. equil == 'rfp2') then
cccc        bcsq = bcond
cccc        where (bcsq == DEF) bcsq = -EXT !On covariant components
cccc        gv%aux%vec_list(IBCNV)%bconds(:,2) = bcsq
cccc
cccc        bcsq = bcond
cccc        where (bcsq == DEF) bcsq = -EXT !On covariant components
cccc        gv%aux%vec_list(IBCNV)%bconds(:,3) = bcsq
cccc      else
cc        bcsq = bcond
cc        where (bcsq == DEF) bcsq = -NEU !On covariant components
cc        gv%aux%vec_list(IBCNV)%bconds(:,2) = bcsq
cc
cc        bcsq = bcond
cc        where (bcsq == DEF) bcsq = -NEU !On covariant components
cc        gv%aux%vec_list(IBCNV)%bconds(:,3) = bcsq
cccc      endif
cc
cc      gv%aux%vec_list(IBCOV)%bconds = gv%aux%vec_list(IBCNV)%bconds
cc
cc      !Current
cccc      if (equil == 'spnch' .or. equil == 'rfp2') then
cccc        bcsq = bcond
cccc        where (bcsq == DEF) bcsq = DIR
cccc        bcnd(:,1) = bcsq
cccc
cccc        bcsq = bcond
cccc        where (bcsq == DEF) bcsq = EXT
cccccc        where (bcsq == DEF) bcsq = -EQU  !On covariant components
cccc        bcnd(:,2) = bcsq
cccc
cccc        bcsq = bcond
cccc        where (bcsq == DEF) bcsq = EXT
cccccc        where (bcsq == DEF) bcsq = -EQU  !On covariant components
cccc        bcnd(:,3) = bcsq
cccc      else
cc        gv%aux%vec_list(IJCNV)%bconds(:,1) = bbcs(:,IAX)
cc        gv%aux%vec_list(IJCNV)%bconds(:,2) = bbcs(:,IAY)
cc        gv%aux%vec_list(IJCNV)%bconds(:,3) = bbcs(:,IAZ)
cccc      endif
cc
cc      gv%aux%vec_list(IJCOV)%bconds = gv%aux%vec_list(IJCNV)%bconds
cc#else
cc      !Current
cc      gv%aux%vec_list(IJCNV)%bconds = bbcs(:,IBX:IBZ)
cc      where (gv%aux%vec_list(IJCNV)%bconds == -NEU)
cccc        gv%aux%vec_list(IJCNV)%bconds = -DIR  !Use covariant components for tangential BCs
cc        gv%aux%vec_list(IJCNV)%bconds = -EQU  !Use covariant components for tangential BCs
cc      end where
cc
cc      gv%aux%vec_list(IJCOV)%bconds = gv%aux%vec_list(IJCNV)%bconds
cc#endif

#if defined(vec_pot)
      !Magnetic field
      bcsq = bcond
      where (bcsq == DEF) bcsq = DIR
      gv%aux%vec_list(IBCNV)%bconds(:,1) = bcsq

cc      if (equil == 'spnch' .or. equil == 'rfp2') then
cc        bcsq = bcond
cc        where (bcsq == DEF) bcsq = -EXT !On covariant components
cc        gv%aux%vec_list(IBCNV)%bconds(:,2) = bcsq
cc
cc        bcsq = bcond
cc        where (bcsq == DEF) bcsq = -EXT !On covariant components
cc        gv%aux%vec_list(IBCNV)%bconds(:,3) = bcsq
cc      else
        bcsq = bcond
        where (bcsq == DEF) bcsq = -NEU !On covariant components
        gv%aux%vec_list(IBCNV)%bconds(:,2) = bcsq

        bcsq = bcond
        where (bcsq == DEF) bcsq = -NEU !On covariant components
        gv%aux%vec_list(IBCNV)%bconds(:,3) = bcsq
cc      endif

      gv%aux%vec_list(IBCOV)%bconds = gv%aux%vec_list(IBCNV)%bconds

      !Current
cc      if (equil == 'spnch' .or. equil == 'rfp2') then
cc        bcsq = bcond
cc        where (bcsq == DEF) bcsq = DIR
cc        bcnd(:,1) = bcsq
cc
cc        bcsq = bcond
cc        where (bcsq == DEF) bcsq = EXT
cccc        where (bcsq == DEF) bcsq = -EQU  !On covariant components
cc        bcnd(:,2) = bcsq
cc
cc        bcsq = bcond
cc        where (bcsq == DEF) bcsq = EXT
cccc        where (bcsq == DEF) bcsq = -EQU  !On covariant components
cc        bcnd(:,3) = bcsq
cc      else
        gv%aux%vec_list(IJCNV)%bconds(:,1) = u_0%array_var(IAX)%bconds
        gv%aux%vec_list(IJCNV)%bconds(:,2) = u_0%array_var(IAY)%bconds
        gv%aux%vec_list(IJCNV)%bconds(:,3) = u_0%array_var(IAZ)%bconds
cc      endif

      gv%aux%vec_list(IJCOV)%bconds = gv%aux%vec_list(IJCNV)%bconds
#else
      !Current
      gv%aux%vec_list(IJCNV)%bconds(:,1) = u_0%array_var(IBX)%bconds
      gv%aux%vec_list(IJCNV)%bconds(:,2) = u_0%array_var(IBY)%bconds
      gv%aux%vec_list(IJCNV)%bconds(:,3) = u_0%array_var(IBZ)%bconds

      where (gv%aux%vec_list(IJCNV)%bconds == -NEU)
        gv%aux%vec_list(IJCNV)%bconds = -EQU  !Use covariant components for tangential BCs
      end where

      gv%aux%vec_list(IJCOV)%bconds = gv%aux%vec_list(IJCNV)%bconds
#endif

c     End program

      end subroutine defineAuxBCs

      end module auxiliaryVariables
