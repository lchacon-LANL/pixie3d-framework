c get_var_names
c######################################################################
      subroutine get_var_names(patch_var,depVarLabels,
     &  auxScalarLabels,auxVectorLabels)

c----------------------------------------------------------------------
c     Initializes MG and creates grid
c----------------------------------------------------------------------

      use parameters

      use grid

      use variables

      use timeStepping

      use newtongm

c      use constants

      use iosetup

c      use generalpurposefunctions

      implicit none

c Call variables

      type(patch), TARGET :: patch_var
      CHARACTER depVarLabels(21,patch_var%u_0%nvar)
      CHARACTER auxScalarLabels(21,patch_var%aux%nvar)
      CHARACTER auxVectorLabels(21,patch_var%aux%nvec)
      integer i, j

c Copy the character arrays

      depVarLabels(:,:) = CHAR(0)
      do i = 1, patch_var%u_0%nvar
         do j = 1, MIN(20,LEN_TRIM(patch_var%u_0%array_var(i)%descr))
            depVarLabels(j,i) = patch_var%u_0%array_var(i)%descr(j:j)
         enddo
      enddo
      auxScalarLabels(:,:) = CHAR(0)
      do i = 1, patch_var%aux%nvar
         do j = 1, MIN(20,LEN_TRIM(patch_var%aux%var_list(i)%descr))
            auxScalarLabels(j,i) = patch_var%aux%var_list(i)%descr(j:j)
         enddo
      enddo
      auxVectorLabels(:,:) = CHAR(0)
      do i = 1, patch_var%aux%nvec
         do j = 1, MIN(20,LEN_TRIM(patch_var%aux%vec_list(i)%descr))
            auxVectorLabels(j,i) = patch_var%aux%vec_list(i)%descr(j:j)
         enddo
      enddo


c End program

      end subroutine get_var_names


