      module slatec_splines

      integer(4) :: kx,ky,kz,nnx,nny,nnz,dim,flg,order

      real(8),allocatable,dimension(:) :: sx,sy,sz
      real(8), dimension(:),allocatable:: tx,ty,tz,work
      real(8), dimension(:,:,:),allocatable:: bcoef

      real(8)    :: dbvalu,db2val,db3val
      external      dbvalu,db2val,db3val

      end module slatec_splines
