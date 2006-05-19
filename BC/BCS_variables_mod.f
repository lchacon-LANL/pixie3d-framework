c module BCS_variables
c####################################################################
      module BCS_variables

        use grid

        use icond

        use constants

        use variables

        use error

        use generalPurposeFunctions

        use oned_int

        integer(4) :: nnvar,imax,imin,jmax,jmin,kmax,kmin
        integer(4) :: nxbc,nybc,nzbc,igxbc,igybc,igzbc

        real(8),allocatable,dimension(:,:) :: rhs

        real(8),allocatable,dimension(:,:,:,:) :: v_cov,v_cnv,v0

        logical :: symm=.false.,bc_debug=.false.

      end module BCS_variables
