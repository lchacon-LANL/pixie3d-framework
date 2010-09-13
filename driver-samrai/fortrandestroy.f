c fortrandestroy
c ######################################################################
      subroutine fortrandestroy

      use variables

      use iosetup

      implicit none

c Call variables

c Local variables

c Begin program

      call deallocateStructures

      call destroyGrid(gv%gparams)

cc      deallocate(gv%bc_seq,STAT=ierr)
      if (associated(gv%bc_grp)) then
        do igr=1,size(gv%bc_grp)
          deallocate(gv%bc_grp(igr)%bc_seq,STAT=ierr)
        enddo
        deallocate(gv%bc_grp,STAT=ierr)
      endif

      deallocate(gv,STAT=ierr)

c End program

      end subroutine fortrandestroy
