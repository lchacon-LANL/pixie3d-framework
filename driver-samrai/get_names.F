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
      integer i, j, k
      logical used
      CHARACTER*(*) name_NULL
      PARAMETER(name_NULL='NULL')

c Copy the character arrays

      depVarLabels(:,:) = CHAR(0)
      do i = 1, patch_var%u_0%nvar
         do j = 1, MIN(20,LEN_TRIM(patch_var%u_0%array_var(i)%descr))
            depVarLabels(j,i) = patch_var%u_0%array_var(i)%descr(j:j)
         enddo
      enddo
      auxScalarLabels(:,:) = CHAR(0)
      do i = 1, patch_var%aux%nvar
c        Determine if this variable is used
         used = .false.
         do j = 1, patch_var%nbc_grp
            do k = 1, patch_var%bc_grp(j)%nbc_seq
               if ( i==-patch_var%bc_grp(j)%bc_seq(k,1) .AND. 
     &               patch_var%bc_grp(j)%bc_seq(k,2)==0 ) then
                  used = .true.
               endif
            enddo
         enddo
         if ( used ) then
           do j = 1, MIN(20,LEN_TRIM(patch_var%aux%var_list(i)%descr))
             auxScalarLabels(j,i) = patch_var%aux%var_list(i)%descr(j:j)
           enddo
         else
           do j = 1, MIN(20,LEN_TRIM(name_NULL))
             auxScalarLabels(j,i) = name_NULL(j:j)
           enddo
         endif
      enddo
      auxVectorLabels(:,:) = CHAR(0)
      do i = 1, patch_var%aux%nvec
         used = .false.
         do j = 1, patch_var%nbc_grp
            do k = 1, patch_var%bc_grp(j)%nbc_seq
               if ( i==-patch_var%bc_grp(j)%bc_seq(k,1) .AND. 
     &               patch_var%bc_grp(j)%bc_seq(k,2)==1 ) then
                  used = .true.
               endif
            enddo
         enddo
         if ( used ) then
           do j = 1, MIN(20,LEN_TRIM(patch_var%aux%vec_list(i)%descr))
             auxVectorLabels(j,i) = patch_var%aux%vec_list(i)%descr(j:j)
           enddo
         else
           do j = 1, MIN(20,LEN_TRIM(name_NULL))
             auxVectorLabels(j,i) = name_NULL(j:j)
           enddo
         endif
      enddo


c End program

      end subroutine get_var_names


