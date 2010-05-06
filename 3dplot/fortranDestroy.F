c fortranDestroy
c ######################################################################
      subroutine fortranDestroy

      use variables

      use iosetup

      use grid

      use graphics

      implicit none

c Call variables

c Local variables

      integer :: ierr,igr

c Begin program

      deallocate(zeros,vzeros,ones)

      call deallocateStructures

      call deallocateDerivedType(u_ic)

      call deallocateApplicationVariables   !External

      call deallocateAuxStruct(gv%aux)

      call destroyGrid(gv%gparams)

cc      deallocate(gv%bc_seq,STAT=ierr)
      do igr=1,size(gv%bc_grp)
        deallocate(gv%bc_grp(igr)%bc_seq,STAT=ierr)
      enddo
      deallocate(gv%bc_grp,STAT=ierr)

      deallocate(gv,STAT=ierr)

c End program

      end subroutine fortranDestroy
