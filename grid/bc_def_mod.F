
c module bc_def
c #####################################################################
      module bc_def

c     -----------------------------------------------------------------
c     This module defines the types of BCs possible:
c     Topological BCs (set topology of computational domain):
c       PER: Periodic
c       SP : Singular point (as in cylindrical, toroidal)
c       SYM: Symmetry (> 0 -> homogeneous Neumann
c                     ,< 0 -> homogeneous Dirichlet)
c     Standard BCs:
c       EQU: Inhomogeneous dirichlet (generally imposed by equilibrium)
c       DIR: Homogeneous dirichlet
c       NEU: Neumann (homogeneous/inhomogeneous)
c       SYM: Symmetry (odd,even depend on sign)
c       EXT: By extrapolation
c       DEF: Default (does nothing)
c     The actual integers chosen indicate a particular order that the BC
c     scheduler follows to impose the boundary conditions.
c     -----------------------------------------------------------------

        use grid_structures

        implicit none

        integer(4) :: PER,DIR,NEU,SP,EQU,DEF,EXT,SYM,IFC
        parameter (DEF=0,PER=1,EQU=2,EXT=3,NEU=4,SYM=5,DIR=6,SP=7,IFC=8)

        integer(4),parameter :: BCLIM=8

        integer(4) :: bcond(6)

        logical    :: asm_dir(3)=.false.,asm=.false.,par_bc=.false.

        real(8)    :: SP_flsv=1d-5
        
      contains

c     bcSP
c     #################################################################
      function bcSP(ibc) result(sing_bc)

        integer(4),optional :: ibc
        logical    :: sing_bc

        integer(4) :: iibc

        if (PRESENT(ibc)) then
          iibc = ibc
        else
          iibc = 1
        endif

        sing_bc = (bcond(iibc) == SP)

      end function bcSP

c     isSP
c     #################################################################
      function isSP(i,j,k,igx,igy,igz) result(sing_point)

        integer(4) :: i,j,k,igx,igy,igz
        logical    :: sing_point

        sing_point = (i+grid_params%ilo(igx)-1==1).and.bcSP()

      end function isSP

c     isSYM
c     #################################################################
      function isSYM(i,igr,dim,loc)

        integer(4) :: i,igr,dim,loc
        logical    :: isSYM

        integer(4) :: ibc,ilog,itst

        ibc = (1+loc)+2*(dim-1)

        isSYM=.false.

        select case(dim)
        case(1)
          ilog = grid_params%ilo(igr)
          select case(loc)
          case (0)
            itst = 1
          case (1)
            itst = grid_params%nxgl(igr)
          end select
        case(2)
          ilog = grid_params%jlo(igr)
          select case(loc)
          case (0)
            itst = 1
          case (1)
            itst = grid_params%nygl(igr)
          end select
        case(3)
          ilog = grid_params%klo(igr)
          select case(loc)
          case (0)
            itst = 1
          case (1)
            itst = grid_params%nzgl(igr)
          end select
        end select

        if (i+ilog-1 == itst .and. bcond(ibc) == SYM)  isSYM=.true.

      end function isSYM

c     alt__eom
c     #################################################################
      function alt__eom()

        logical :: alt__eom

        alt__eom=.false.
        if (bcSP()) alt__eom=.true.

      end function alt__eom

      end module bc_def
