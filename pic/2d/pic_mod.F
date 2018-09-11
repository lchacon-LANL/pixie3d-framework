#include "pic_defs.h"

#if !defined(plotting)
c pic_mover
c #####################################################################    
      module pic_mover

      use pic_map

      use pic_moments

      use pic_diag

      !Particle orbit iteration counters
      integer(8) :: pst_counter=0    !Total number of particle substeps
     .             ,pit_counter=0    !Total number of picard iterations
     .             ,psh_counter=0    !Total number of particle pushes

      integer(8) :: pit_max=0           !Maximum picard iterations per substep
     .             ,pst_max=0           !Maximum number of particle substeps per timestep
     .             ,ntot_pcles=0        !Total number of particles

      real(8) :: pit_avg=0           !Average picard iterations per substep
     .          ,pst_avg=0           !Average number of particle substeps per timestep

      !Diagnostics
      real(8) :: push_timer=0d0
      integer,private :: clock_rate,clock_start,clock_stop

      contains

c     pic_mover_2d
c     ##################################################################
      subroutine pic_mover_2d(dt,A,phi,DtA,B0,tfv,tfx,cm,pm,cm_exch)
        
      implicit none
        
c     Call variables

      real(8),intent(inout) :: A(0:,0:,:),phi(0:,0:),DtA(0:,0:,:)
      real(8),intent(inout) :: B0(3),dt
      real(8),intent(in),optional :: tfv,tfx
      real(8),intent(out),optional :: cm(3),pm(3),cm_exch(3)

c     Local variables

      integer :: ipc,isp,i0(_Npg),j0(_Npg),i(_Npg),j(_Npg)
      integer :: im,it,k,psh_sum
     .          ,pits,pit_sum
     .          ,psts,pst_sum
      integer(8) :: pitmax,pstmax
      integer :: numpcles,numspcs,nxg,nyg,nzg
      real(_prec) :: d_x (_Npg,3)
      real(_prec) :: v_np(_Npg,3)
      real(_prec) :: x_nh(_Npg,3)

      real(_prec) :: dt_sub(_Npg),dtp(_Npg)
      real(_prec) :: v0(_Npg,3),dx(_Npg,3)
      real(_prec) :: xp(_Npg),yp(_Npg),xfp(_Npg),yfp(_Npg)
     .              ,v2(_Npg),w_pcle(_Npg),m_sp,q_sp

      real(_prec),allocatable,dimension(:,:) :: garr_n  
     $                                         ,garr_vx 
     $                                         ,garr_vy 
     $                                         ,garr_vz 
     $                                         ,garr_vx2
     $                                         ,garr_vy2
     $                                         ,garr_vz2

      real(_prec) :: My2(_Npg,1,3),Mx2(_Npg,3,1),My1(_Npg,1,2)
      real(_prec) :: Mx1(_Npg,2,1),Mxy2(_Npg,3,3)
      real(_prec) :: M_w_prod(_Npg,3,3),dum1(_Npg,3),dum2(_Npg,3)

      real(_prec) :: rhon,rho,jip,ji,ecc,cml(3),pml(3),cmex(3)
     .              ,lacc(9),gacc(9)

      real(_prec) :: dtv,dtx
      real(8) :: ihx,ihy,ihz

      type(particle) :: pcle_nu

c     Begin program

      allocate(garr_n  (0:size(A,1)-1,0:size(A,2)-1)
     $        ,garr_vx (0:size(A,1)-1,0:size(A,2)-1)
     $        ,garr_vy (0:size(A,1)-1,0:size(A,2)-1)
     $        ,garr_vz (0:size(A,1)-1,0:size(A,2)-1)
     $        ,garr_vx2(0:size(A,1)-1,0:size(A,2)-1)
     $        ,garr_vy2(0:size(A,1)-1,0:size(A,2)-1)
     $        ,garr_vz2(0:size(A,1)-1,0:size(A,2)-1))

c     Particle pusher timer - C. Mikida, 05/16

      call system_clock(count_rate=clock_rate) 
      call system_clock(count=clock_start)

c      call start_collection()

      if(present(tfv)) then
         dtv = dt*tfv
      else
         dtv = dt
      end if

      if(present(tfx)) then
         dtx = dt*tfx
      else
         dtx = dt
      end if

      cml = 0d0 ; pml = 0d0 ; cmex = 0d0

      call get_global_grid_size(nxg,nyg,nzg)
      call get_inverse_h(ihx,ihy,ihz)

c     Move particles 
      
      psh_sum = 0
      pst_sum = 0
      pit_sum = 0
      pitmax  = 0
      pstmax  = 0

      numspcs = size(spcs)

      ntot_pcles = 0

      do isp=1, numspcs
         numpcles = size(spcs(isp)%pcles)

         ntot_pcles = ntot_pcles + (numpcles - nip)

         garr_n = 0d0

         garr_vx = 0d0
         garr_vy = 0d0
         garr_vz = 0d0

         garr_vx2= 0d0
         garr_vy2= 0d0
         garr_vz2= 0d0

         m_sp = spcs(isp)%m
         q_sp = spcs(isp)%q

!     note that there is some issue with valgrind and libgomp.
!$omp parallel DEFAULT(SHARED) private(ipc,it,v0,v2,dtp,dt_sub
!$omp.              ,d_x,dx,xp,yp,xfp,yfp,x_nh,dum1,dum2
!$omp.              ,i,j,k,i0,j0,Mx1,Mx2,My1,My2,Mxy2,w_pcle
!$omp.              ,psts,M_w_prod,pits,pcle_nu)
!$         call set_omp_thread_id()
!$omp DO SCHEDULE (DYNAMIC,4)
!$omp.reduction(+:garr_n,garr_vx,garr_vy,garr_vz
!$omp.           ,garr_vx2,garr_vy2,garr_vz2,pml
#if defined(PIC_CM_DIAG)        
!$omp.           ,cmex,cml
#endif
!$omp.           ,pst_sum,psh_sum,pit_sum)
!$omp.reduction(max:pitmax,pstmax)
         do ipc=nip+1,numpcles

           pcle_nu%ijk_n = spcs(isp)%pcles(ipc)%ijk_n
           pcle_nu%x_n   = spcs(isp)%pcles(ipc)%x_n
           pcle_nu%v_n   = spcs(isp)%pcles(ipc)%v_n
           pcle_nu%w_n   = spcs(isp)%pcles(ipc)%w_n

           it = 0

           ! Use local registers for performance (C. Mikida)
           w_pcle = spcs(isp)%pcles(ipc)%w_n

#if defined(PIC_CM_DIAG)
           d_x = 0d0
#endif
           dtp = dtx

           do while(maxval(dtp)>0)
             it = it+1

             !Perform pcle substep
             call ACC_CN_incell_offset(spcs(isp)%q_m
     .                                ,spcs(isp)%pcles(ipc),pcle_nu
     $                                ,dt,dtp,dt_sub
     $                                ,DtA,B0,A,phi,dx
     $                                ,psts,pits)

             !Update pcle counters
             pst_sum = pst_sum + 1      !Particle substep counter

             pitmax  = max(pits,pitmax) !Maximum number of Picard iterations per substep

             pit_sum = pit_sum + pits   !Total number of picard iterations
             psh_sum = psh_sum + psts   !Total number of particle pushes (including ALL picard iterations)

             !Gather moments
#if defined(PIC_CM_DIAG)
             !accumulate particle displacement (for diagnose canonical mom conservation)
             d_x = d_x + dx
#endif

             call xform_pcle_idx(pcle_nu%ijk_n,nxg,i,j)

             x_nh = pcle_nu%x_n + dx*0.5d0
             xp = x_nh(:,1)*ihx
             yp = x_nh(:,2)*ihy 
             
             call shape1_offset(xp,Mx1(:,1:2,1))
             call shape1_offset(yp,My1(:,1,1:2))               

             xp = pcle_nu%x_n(:,1)*ihx
             yp = pcle_nu%x_n(:,2)*ihy

             xfp = (pcle_nu%x_n(:,1)+dx(:,1))*ihx
             yfp = (pcle_nu%x_n(:,2)+dx(:,2))*ihy

             call shape2_offset(xp ,dum1)
             call shape2_offset(xfp,dum2)
             Mx2(:,1:3,1) = (dum2-dum1)*0.5d0 !careful

             call shape2_offset(yp ,dum1)
             call shape2_offset(yfp,dum2)
             My2(:,1,1:3) = (dum2-dum1)*0.5d0

             !VECTOR SIMD
             do k=1,_Npg
#ifdef __GNU
                Mxy2(k,1:3,1:3)= opt_vecvec_3(Mx2(k,:,:),My2(k,:,:))
#else
                Mxy2(k,1:3,1:3)= matmul      (Mx2(k,:,:),My2(k,:,:))
#endif
             end do

             call shape2_offset(xp ,dum1)
             call shape2_offset(xfp,dum2)
             Mx2(:,1:3,1) = (dum1+dum2)*0.5d0 !careful

             call shape2_offset(yp ,dum1)
             call shape2_offset(yfp,dum2)
             My2(:,1,1:3) = (dum1+dum2)*0.5d0

             !Do NOT vectorize: memory race issues
             do k=1, _Npg
#ifdef __GNU
                Mxy2(k,1:3,1:3)= Mxy2(k,1:3,1:3)
     $                          +opt_vecvec_3(Mx2(k,:,:),My2(k,:,:))
#else
                Mxy2(k,1:3,1:3)= Mxy2(k,1:3,1:3)
     $                          +matmul(Mx2(k,:,:),My2(k,:,:))
#endif

#ifdef __GNU
                !accumulate vx (linear in x, 2nd in y)
                garr_vx(i(k)-1:i(k)  ,j(k)-1:j(k)+1) = 
     .          garr_vx(i(k)-1:i(k)  ,j(k)-1:j(k)+1) + 
     .          dx(k,1)*opt_vecvec_23(Mx1(k,:,:),My2(k,:,:))*w_pcle(k)

                !accumulate vy (linear in y, 2nd in x)
                garr_vy(i(k)-1:i(k)+1,j(k)-1:j(k)) = 
     .          garr_vy(i(k)-1:i(k)+1,j(k)-1:j(k)) + 
     .          dx(k,2)*opt_vecvec_32(Mx2(k,:,:),My1(k,:,:))*w_pcle(k)
#else                
                !accumulate vx (linear in x, 2nd in y)
                garr_vx(i(k)-1:i(k)  ,j(k)-1:j(k)+1) = 
     .          garr_vx(i(k)-1:i(k)  ,j(k)-1:j(k)+1) + 
     .          dx(k,1)*matmul(Mx1(k,:,:),My2(k,:,:))*w_pcle(k)

                !accumulate vy (linear in y, 2nd in x)
                garr_vy(i(k)-1:i(k)+1,j(k)-1:j(k)) = 
     .          garr_vy(i(k)-1:i(k)+1,j(k)-1:j(k)) +
     .          dx(k,2)*matmul(Mx2(k,:,:),My1(k,:,:))*w_pcle(k)
#endif
                !accumulate vz (2nd in x,y)
                garr_vz(i(k)-1:i(k)+1,j(k)-1:j(k)+1) =
     .          garr_vz(i(k)-1:i(k)+1,j(k)-1:j(k)+1) +
     .          dx(k,3)*Mxy2(k,:,:)*w_pcle(k)

             end do
                  
             !Prepare for next substep
             pcle_nu%ijk_n = spcs(isp)%pcles(ipc)%ijk_np
             pcle_nu%x_n   = spcs(isp)%pcles(ipc)%x_np
             pcle_nu%v_n   = spcs(isp)%pcles(ipc)%v_np

             dtp = dtp-dt_sub

           end do               !sub-stepping

           !Maximum number of substeps per timestep
           pstmax  = max(it,pstmax)
             
           !2nd order (for checking charge conservation) careful
           xp = pcle_nu%x_n(:,1)*ihx
           yp = pcle_nu%x_n(:,2)*ihy

           call shape1_offset(xp,Mx1(:,1:2,1))
           call shape1_offset(yp,My1(:,1,1:2))               
           call shape2_offset(xp,Mx2(:,1:3,1))
           call shape2_offset(yp,My2(:,1,1:3))

           !VECTOR SIMD
           do k=1,_Npg
#ifdef __GNU
              Mxy2(k,1:3,1:3)= opt_vecvec_3(Mx2(k,:,:),My2(k,:,:))
#else
              Mxy2(k,1:3,1:3)= matmul(Mx2(k,:,:),My2(k,:,:))
#endif
              M_w_prod(k,:,:) = Mxy2(k,:,:)*w_pcle(k)
           enddo

           !Do NOT put SIMD statement, reduction!!!!
           v0 = pcle_nu%v_n
           call xform_pcle_idx(pcle_nu%ijk_n,nxg,i0,j0)

           !Do NOT vectorize: memory race issues
           do k=1,_Npg
              !accumulate n (quadratic in x,y)
              garr_n(i0(k)-1:i0(k)+1,j0(k)-1:j0(k)+1) = 
     .        garr_n(i0(k)-1:i0(k)+1,j0(k)-1:j0(k)+1) + M_w_prod(k,:,:)

              !accumulate v^2 (2nd in x,y, cell-centered)
              v2(k) = v0(k,1)*v0(k,1)
              garr_vx2(i0(k)-1:i0(k)+1,j0(k)-1:j0(k)+1) = 
     .        garr_vx2(i0(k)-1:i0(k)+1,j0(k)-1:j0(k)+1) + 
     .                             v2(k)*M_w_prod(k,:,:)

              v2(k) = v0(k,2)*v0(k,2)
              garr_vy2(i0(k)-1:i0(k)+1,j0(k)-1:j0(k)+1) = 
     .        garr_vy2(i0(k)-1:i0(k)+1,j0(k)-1:j0(k)+1) + 
     .                             v2(k)*M_w_prod(k,:,:)

              v2(k) = v0(k,3)*v0(k,3)
              garr_vz2(i0(k)-1:i0(k)+1,j0(k)-1:j0(k)+1) = 
     .        garr_vz2(i0(k)-1:i0(k)+1,j0(k)-1:j0(k)+1) + 
     .                             v2(k)*M_w_prod(k,:,:)
           end do 

           !Diagnostics
           !Do NOT vectorize: memory race issues
           do k=1,_Npg
              !Gather momenta
              pml = pml + m_sp*v0(k,:)*w_pcle(k)

#if defined(PIC_CM_DIAG)
              !Gather canonical momenta
              cml = cml + m_sp*v0(k,:)
     $             + cmp(spcs(isp)%q,i0(k),j0(k),xp(k),yp(k),A(:,:,:))
     .      + 0.5*dt*cmp(spcs(isp)%q,i0(k),j0(k),xp(k),yp(k),DtA(:,:,:))

              !Canonical momenta exchange
              cmex(2) = cmex(2)+q_sp*B0(3)*d_x(k,1)
     .                         -q_sp*B0(1)*d_x(k,3)
              cmex(3) = cmex(3)+q_sp*B0(1)*d_x(k,2)
     .                         -q_sp*B0(2)*d_x(k,1)
#endif
           end do 

         end do                 !particles
!$omp end do 
!$omp end parallel 

         !orbit-averaged
         moms_g(:,:,1,isp,Ipx) = garr_vx/dtx
         moms_g(:,:,1,isp,Ipy) = garr_vy/dtx
         moms_g(:,:,1,isp,Ipz) = garr_vz/dtx

         !non-orbit-averaged
         moms_g(:,:,1,isp,Ins ) = garr_n !careful        
         moms_g(:,:,1,isp,IS11) = garr_vx2
         moms_g(:,:,1,isp,IS22) = garr_vy2
         moms_g(:,:,1,isp,IS33) = garr_vz2

      end do                    !species

      deallocate(garr_n  
     $          ,garr_vx 
     $          ,garr_vy 
     $          ,garr_vz 
     $          ,garr_vx2
     $          ,garr_vy2
     $          ,garr_vz2)

c     Particle pusher timing end - C. Mikida, 05/16
      
c      call stop_collection()
      call system_clock(count=clock_stop)

c     Append to particle pusher timer/counters

      push_timer= push_timer + real(clock_stop - clock_start)
     .                        /real(clock_rate)

c     Synchronize moments

      call wrap_moments_nodes(moms_g(:,:,1,:,:),wrap_density=.false.)
      moms_g(:,:,0,:,:) = moms_g(:,:,1,:,:)
      moms_g(:,:,2,:,:) = moms_g(:,:,1,:,:)

c     Scatter moments to local mesh

      do im=0,size(moms,5)-1
       call find_local(gpicdef,1,moms_g(:,:,:,:,im),moms(:,:,:,:,im))
      enddo

c     Apply BCs to moments

      call applyBC_moms(gpicdef)

c     Output info

#if defined(petsc)
      lacc(1:3) = cml
      lacc(4:6) = pml
      lacc(7:9) = cmex
      call MPI_AllReduce(lacc,gacc,size(lacc),MPI_DOUBLE_PRECISION
     $                  ,MPI_SUM,MPI_COMM_WORLD,mpierr)
      cml = gacc(1:3)
      pml = gacc(4:6)
      cmex= gacc(7:9)
#endif

      if (PRESENT(cm)) cm = cml
      if (PRESENT(pm)) pm = pml
      if (PRESENT(cm_exch)) cm_exch = cmex

      cm_ex = cmex

c     Substep counter

      pit_counter = pit_counter + pit_sum  !Total number of picard iterations
      psh_counter = psh_counter + psh_sum  !Total number of particle pushes
      pst_counter = pst_counter + pst_sum  !Total number of particle substeps

      pit_max = pitmax                     !Maximum number of picard iterations per substep
      pst_max = pstmax                     !Maximum number of particle substeps per timestep

      pit_avg = 1d0*pit_sum/pst_sum        !Average number of picard iterations per substep
      pst_avg = 1d0*pst_sum/ntot_pcles     !Average number of substeps per particle

      ntot_pcles = ntot_pcles*_Npg         !Account for pcle grouping for external averages

c     End program

      end subroutine pic_mover_2d

c     ACC_CN_incell_offset
c     ##################################################################
      subroutine ACC_CN_incell_offset (q_m,pcle,pcle_nu,dt,dtp,dt_sub
     $                                ,DtA,B0,A,phi,dx,substeps,it)

      implicit none
      
c     Call variables

      type(particle), intent(inout) :: pcle,pcle_nu

      real(8),intent(in) :: dt,q_m
      real(8),intent(in) :: A(0:,0:,:),phi(0:,0:),DtA(0:,0:,:),B0(3)

      real(_prec),intent(inout) :: dtp(_Npg)
      real(_prec),intent(out)   :: dt_sub(_Npg),dx(_Npg,3)

      integer,intent(out):: substeps,it

c     Local variables

      integer :: i,j,k,nxg,nyg,nzg
      integer :: i0(_Npg),j0(_Npg),inp_loc(_Npg),jnp_loc(_Npg)

      real(8) :: ihx,ihy,ihz,hx,hy,hz

      real(_prec) :: rx(_Npg),ry(_Npg),rz(_Npg)

      real(_prec), parameter :: eps  = 1d-12 !Buffer for dt computation in pcle Xing
     .                         ,eps2 = 1d-8  !Buffer for pcle Xing

      real(_prec) :: ptola,alpha=0d0,dt_coeff=0.2

      real(_prec) :: am(_Npg,_Dim)

      real(_prec) :: x_k(_Npg),x_nh(_Npg),y_nh(_Npg),d_err(_Npg)
      real(_prec) :: dx_k(_Npg,3),vh_k(_Npg,3),dx_err
      
      real(_prec) :: ac0(_Npg,3),v0curv(_Npg,3),vh(_Npg,3),vhes(_Npg,3)
     $              ,vhcurv(_Npg,3)

      real(_prec) :: Omegac(_Npg),omb(_Npg),qm(_Npg)
     $              ,Bmax(_Npg),vdotB(_Npg)

      real(_prec) :: Ep(_Npg,3),Bp(_Npg,3)
      real(_prec) :: xnp(_Npg),ynp(_Npg)
      real(_prec) :: dxl(_Npg,2),dxr(_Npg,2)

      real(_prec) :: DtA_1 (_Npg,2,3),DtA_2 (_Npg,3,2),DtA_3 (_Npg,3,3)
      real(_prec) :: A_1   (_Npg,2,3),A_2   (_Npg,3,2),A_3   (_Npg,3,3)
      real(_prec) :: E_1   (_Npg,2,3),E_2   (_Npg,3,2),E_3   (_Npg,3,3)

      real(_prec) :: phi_loc(_Npg,3,3),xnp_loc(_Npg,3),vnp_loc(_Npg,3)
      real(_prec) :: x0_loc(_Npg,3),v0_loc(_Npg,3)
      real(_prec) :: x0p(_Npg),y0p(_Npg)
      real(_prec) :: Sx2_o(_Npg,3),Sy2_o(_Npg,3)
      real(_prec) :: dteq0(_Npg)

c     Begin program

      qm = q_m

      dteq0 = 0d0

      pcle%x_np = pcle_nu%x_n
      pcle%v_np = pcle_nu%v_n
      pcle%ijk_np = pcle_nu%ijk_n
      pcle%w_np = pcle_nu%w_n

      call get_global_grid_size(nxg,nyg,nzg)
      call get_inverse_h(ihx,ihy,ihz)
      call get_h(hx,hy,hz)

c     Initialize local registers

      call xform_pcle_idx(pcle_nu%ijk_n,nxg,i0,j0)

      x0_loc = pcle_nu%x_n
      v0_loc = pcle_nu%v_n

c     Also define some insets of arrays for later use (C. Mikida)

      !VECTOR SIMD
      do i=1,_Npg
         DtA_1(i,:,:) = DtA(i0(i)-1:i0(i)  , j0(i)-1:j0(i)+1, 1)
         DtA_2(i,:,:) = DtA(i0(i)-1:i0(i)+1, j0(i)-1:j0(i)  , 2)
         DtA_3(i,:,:) = DtA(i0(i)-1:i0(i)+1, j0(i)-1:j0(i)+1, 3)

         A_1(i,:,:)   = A(i0(i)-1:i0(i)  , j0(i)-1:j0(i)+1, 1)
         A_2(i,:,:)   = A(i0(i)-1:i0(i)+1, j0(i)-1:j0(i)  , 2)
         A_3(i,:,:)   = A(i0(i)-1:i0(i)+1, j0(i)-1:j0(i)+1, 3)

         phi_loc(i,:,:) = phi(i0(i)-1:i0(i)+1,j0(i)-1:j0(i)+1)

         E_1(i,:,:) =-DtA_1(i,:,:)
     .               -(phi_loc(i,2:3,:)-phi_loc(i,1:2,:))*ihx
         E_2(i,:,:) =-DtA_2(i,:,:)
     .               -(phi_loc(i,:,2:3)-phi_loc(i,:,1:2))*ihy
         E_3(i,:,:) =-DtA_3(i,:,:)
      end do

c     Estimate sub-time-step 

      dt_sub = dtp

      call estimate_dtp(dt_sub)

c     Initialize Picard

      call initialize_it(dt_sub,dx,vh)

      it=0

      substeps = 0

      do
        it = it+1

        vh_k = vh
        dx_k = dx

        substeps = substeps + _Npg

        !Perform Boris push
        call Boris_solve(it,x0_loc,v0_loc,dx,vh,dt_sub)

        !Check convergence (very expensive!!!)
        if (mod(it,2) /= 0) cycle  !Check everyother iteration

#if defined(L1_NORM)
        dx_err = 0d0
        do i=1,_Npg
           dx_err = max(dx_err, abs(dx(i,1)-dx_k(i,1))
     $                         +abs(dx(i,2)-dx_k(i,2))
     $                         +abs(dx(i,3)-dx_k(i,3))
     $                         +abs(vh(i,1)-vh_k(i,1))
     $                         +abs(vh(i,2)-vh_k(i,2))
     $                         +abs(vh(i,3)-vh_k(i,3)))
        enddo
#else
c$$$        dx_err = (sum(abs(dx-dx_k))
c$$$     $           +sum(abs(vh-vh_k)))/(2*_Npg)
        dx_err = max(maxval(abs(dx-dx_k)),maxval(abs(vh-vh_k)))
#endif

#if _prec == 8
        if (dx_err < 1d-13 .or. it>500) exit
#else
        if (dx_err < 1d-6  .or. it>500) exit
#endif
      end do             

c     Update particle quantities

      vnp_loc = v0_loc + 2d0*(vh - v0_loc) !2d0*vh - v0
      xnp_loc = x0_loc + dx

      inp_loc = i0
      jnp_loc = j0

c     Cell crossing and particle b.c.

      call cell_crossing(dx,xnp_loc,vnp_loc,inp_loc,jnp_loc)

c     Back to the global registers now:

      pcle%ijk_np = inp_loc + nxg*(jnp_loc - 1)

      pcle%v_np = vnp_loc
      pcle%x_np = xnp_loc

c     End program

      contains

c     initialize_it
c     ######################################################################
      subroutine initialize_it(dts,dx,vh)

      implicit none

c     Call variables

      real(_prec) :: dts(_Npg),dx(_Npg,3),vh(_Npg,3)

c     Local variables

      real(_prec) :: dx_pn,dv_pn(_Npg),vha(_Npg),ama(_Npg),am(_Npg,_Dim)
      real(_prec) :: vhcurv_inv(_Npg,3)
cc      real(_prec) :: Chrtof2(_Npg,3,3,3),ainert(_Npg,3)

c     Begin program

c     Find vh

      vh = (v0_loc+pcle%v_np)*0.5d0 

c     Find dx

      if (no_map) then
        vhcurv = vh
      else
        vhcurv = find_vcurv(i0,j0,x0_loc,vh)

c WARNING: CHRISTOFFEL SYMBOL NOT IMPLEMENTED IN 3D.

c$$$           vhcurv = find_vcurv(i0,j0,x0_loc,vh,gamma=Chrtof2)
c$$$           ainert(:,:) = 0d0
c$$$
c$$$           do i=1,3
c$$$              do j=1,3
c$$$                !VECTOR SIMD
c$$$                 do ii=1,3
c$$$                    ainert(ii,:)=ainert(ii,:)
c$$$     .                      +vhcurv(:,j)*Chrtof2(ii,i,j,:)*vhcurv(i,:)
c$$$                 end do
c$$$              end do
c$$$           end do
      endif

      vhcurv_inv = 1d0/(vhcurv + 1d-44)

      dxl(:,1) =   -x0_loc(:,1)
      dxr(:,1) = hx-x0_loc(:,1)
            
      dxl(:,2) =   -x0_loc(:,2)
      dxr(:,2) = hy-x0_loc(:,2)

      dts = min(dts
     .            ,max((dxl(:,1)-eps2)*vhcurv_inv(:,1)
     .                ,(dxr(:,1)+eps2)*vhcurv_inv(:,1))
     .            ,max((dxl(:,2)-eps2)*vhcurv_inv(:,2)
     .                ,(dxr(:,2)+eps2)*vhcurv_inv(:,2)))

      dx(:,1) = vhcurv(:,1)*dts
      dx(:,2) = vhcurv(:,2)*dts
      dx(:,3) = vhcurv(:,3)*dts

c     Guess for dx error

c$$$      if (no_map) then
c$$$        call scatter_fields_nomap(x0_loc,dteq0,dx,A_1,A_2,A_3
c$$$     $                           ,dtA_1,dtA_2,dtA_3,Ep,Bp)
c$$$      else
c$$$        call scatter_fields_curv (x0_loc,dteq0,dx,A_1,A_2,A_3
c$$$     $                           ,dtA_1,dtA_2,dtA_3,Ep,Bp)
c$$$      endif
c$$$
c$$$      am(:,1) = qm*(Ep(:,1) + vh(:,2)*Bp(:,3) - vh(:,3)*Bp(:,2))
c$$$      am(:,2) = qm*(Ep(:,2) + vh(:,3)*Bp(:,1) - vh(:,1)*Bp(:,3))
c$$$      am(:,3) = qm*(Ep(:,3) + vh(:,1)*Bp(:,2) - vh(:,2)*Bp(:,1))
c$$$
c$$$      !VECTOR SIMD
c$$$      do i=1, _Npg
c$$$#if defined(L1_NORM)
c$$$        ama(i) = sum(abs(am(i,:))) !L1 norm
c$$$#else
c$$$        ama(i) = sqrt(dot_product(am(i,:),am(i,:))) !L2 norm
c$$$#endif
c$$$      end do
c$$$
c$$$      dv_pn = dts*ama
c$$$
c$$$#if defined(L1_NORM)
c$$$      vha = sum(abs(vh),2) !L1 norm
c$$$#else
c$$$      vha = sqrt(vh(:,1)*vh(:,1) + vh(:,2)*vh(:,2) + vh(:,3)*vh(:,3)) !L2 norm
c$$$#endif
c$$$
c$$$      dv_pn = max(dts*vha + 0.5*dts**2*ama,dv_pn)
c$$$      dx_pn = maxval(dv_pn)

      end subroutine initialize_it

c     estimate_dtp
c     ######################################################################
      subroutine estimate_dtp(dt_sub)

      implicit none
      
c     Call variables

        real(_prec) :: dt_sub(_Npg)

c     Logical variables

        integer :: i,j,k
        real(_prec) :: dts(_Npg)    !Test timestep to estimate qtys
     .                ,xtoll=1d-3   !Tolerance for timestep estimator
        real(_prec) :: iDenom(_Npg),at(_Npg,3),accl(_Npg,3),anup(_Npg)
        real(_prec) :: qmdth(_Npg)

c     Begin program

        dts = 1d-8     !Test timestep to estimate qtys

c     Scatter fields

        dx = 0d0

        if (no_map) then
          call scatter_fields_nomap(x0_loc,dteq0,dx,A_1,A_2,A_3
     $                             ,dtA_1,dtA_2,dtA_3,Ep,Bp)
        else
          call scatter_fields_curv (x0_loc,dteq0,dx,A_1,A_2,A_3
     $                             ,dtA_1,dtA_2,dtA_3,Ep,Bp)
        endif

c     Cyclotron freq.

#if _Dim>1
#if defined(L1_NORM)
        Bmax = sum(abs(Bp),2)     !L1 norm
#else
        Bmax = sqrt(sum(Bp**2,2)) !L2 norm
#endif
        Omegac = qm * Bmax
        Omegac = abs(Omegac) + 1d-44  !careful
        omb = 1d0/Omegac
        dt_sub = min(dt_sub,dt_coeff*omb)
#endif

c     Initial acceleration

        vh = v0_loc

        !VECTOR SIMD
        do i=1,_Npg
           ac0(i,1)= qm(i)*(Ep(i,1) + vh(i,2)*Bp(i,3) - vh(i,3)*Bp(i,2))
           ac0(i,2)= qm(i)*(Ep(i,2) + vh(i,3)*Bp(i,1) - vh(i,1)*Bp(i,3))
           ac0(i,3)= qm(i)*(Ep(i,3) + vh(i,1)*Bp(i,2) - vh(i,2)*Bp(i,1))
        enddo

c     Test for dts (presumably no cell-crossing with a small timestep)

        if (no_map) then
           vhcurv = vh
        else
           vhcurv = find_vcurv(i0,j0,x0_loc,vh)
        endif
        v0curv = vhcurv

        dx(:,1) = dts*vhcurv(:,1)
        dx(:,2) = dts*vhcurv(:,2)
        dx(:,3) = dts*vhcurv(:,3)

        if (no_map) then
          call scatter_fields_nomap(x0_loc,dts,dx,A_1,A_2,A_3
     $                             ,dtA_1,dtA_2,dtA_3,Ep,Bp) !mid-point
        else
          call scatter_fields_curv (x0_loc,dts,dx,A_1,A_2,A_3
     $                             ,dtA_1,dtA_2,dtA_3,Ep,Bp) !mid-point
        endif
        
        !Velocity update
        qmdth = qm*dts*0.5d0
        vhes(:,1) = v0_loc(:,1) + qmdth*Ep(:,1)
        vhes(:,2) = v0_loc(:,2) + qmdth*Ep(:,2)
        vhes(:,3) = v0_loc(:,3) + qmdth*Ep(:,3)

        !VECTOR SIMD
        do i=1,_Npg
           iDenom(i) = (1d0
     $                 + qmdth(i)*qmdth(i)*dot_product(Bp(i,:),Bp(i,:)))
           iDenom(i) = 1d0/iDenom(i)
        end do

        vdotB = (vhes(:,1)*Bp(:,1)+vhes(:,2)*Bp(:,2)+vhes(:,3)*Bp(:,3))

        vdotB = qmdth*vdotB

        vh(:,1) = (vhes(:,1)+qmdth*(vhes(:,2)*Bp(:,3)-vhes(:,3)*Bp(:,2)
     .                            + vdotB*Bp(:,1)))*iDenom
        vh(:,2) = (vhes(:,2)+qmdth*(vhes(:,3)*Bp(:,1)-vhes(:,1)*Bp(:,3) 
     .                            + vdotB*Bp(:,2)))*iDenom
        vh(:,3) = (vhes(:,3)+qmdth*(vhes(:,1)*Bp(:,2)-vhes(:,2)*Bp(:,1) 
     .                            + vdotB*Bp(:,3)))*iDenom

        at(:,1) = qm*(Ep(:,1) + vh(:,2)*Bp(:,3) - vh(:,3)*Bp(:,2))
        at(:,2) = qm*(Ep(:,2) + vh(:,3)*Bp(:,1) - vh(:,1)*Bp(:,3))
        at(:,3) = qm*(Ep(:,3) + vh(:,1)*Bp(:,2) - vh(:,2)*Bp(:,1))

        !VECTOR SIMD
        do i=1,_Npg 
          accl(i,:) = (at(i,:) - ac0(i,:))*2d0 !/dts added below in dt_sub
#if defined(L1_NORM)
          anup(i) = sum(abs(accl(i,:)))                    !L1 norm
#else
          anup(i) = sqrt(dot_product(accl(i,:),accl(i,:))) !L2 norm
#endif
        end do

        dt_sub = min(dt_sub,(6d0*xtoll*dts/(anup+1d-30*dts))
     $                      **(0.333333333333333d0))

      end subroutine estimate_dtp

c     Boris_solve
c     ######################################################################
      subroutine Boris_solve(it,x0,v0,dx,vh,dt_sub) 

      implicit none
      
c     Call variables

        integer :: it
        real(_prec) :: x0(_Npg,3)
        real(_prec) :: v0(_Npg,3)
        real(_prec) :: dx(_Npg,3)
        real(_prec) :: vh(_Npg,3),dt_sub(_Npg)

c     Logical variables

        real(_prec) :: vhes(_Npg,3),vhcurv(_Npg,3),vhcurv_inv(_Npg,3)
        real(_prec) :: qmdth(_Npg),iDenom(_Npg),dtb(_Npg,3)
        real(_prec) :: xp_vhcurv(_Npg,3)

c     Begin program

        if (no_map) then
          call scatter_fields_nomap(x0,dt_sub,dx,A_1,A_2,A_3
     $                           ,dtA_1,dtA_2,dtA_3,Ep,Bp)
        else
          call scatter_fields_curv (x0,dt_sub,dx,A_1,A_2,A_3
     $                           ,dtA_1,dtA_2,dtA_3,Ep,Bp)
        endif

        !Velocity update
        qmdth = qm*dt_sub*0.5d0
        vhes(:,1) = v0(:,1) + qmdth*Ep(:,1)
        vhes(:,2) = v0(:,2) + qmdth*Ep(:,2)
        vhes(:,3) = v0(:,3) + qmdth*Ep(:,3)

        !VECTOR SIMD
        do i=1,_Npg
           iDenom(i) = (1d0
     $                 + qmdth(i)*qmdth(i)*dot_product(Bp(i,:),Bp(i,:)))
           iDenom(i) = 1d0/iDenom(i)
        end do

        vdotB = (vhes(:,1)*Bp(:,1)+vhes(:,2)*Bp(:,2)+vhes(:,3)*Bp(:,3))

        vdotB = qmdth*vdotB

        vh(:,1) = (vhes(:,1)+qmdth*(vhes(:,2)*Bp(:,3)-vhes(:,3)*Bp(:,2)
     .                            + vdotB*Bp(:,1)))*iDenom
        vh(:,2) = (vhes(:,2)+qmdth*(vhes(:,3)*Bp(:,1)-vhes(:,1)*Bp(:,3) 
     .                            + vdotB*Bp(:,2)))*iDenom
        vh(:,3) = (vhes(:,3)+qmdth*(vhes(:,1)*Bp(:,2)-vhes(:,2)*Bp(:,1) 
     .                            + vdotB*Bp(:,3)))*iDenom

        !Time step estimate for cell crossing
        if (no_map) then
          vhcurv = vh
        else
          xp_vhcurv = x0 + 0.5*dx
          vhcurv = find_vcurv(i0,j0,xp_vhcurv,vh)
        endif

        dx(:,1) = dt_sub*vhcurv(:,1)
        dx(:,2) = dt_sub*vhcurv(:,2)
        dx(:,3) = dt_sub*vhcurv(:,3)

        !VECTOR SIMD
        do i=1,_Npg
           !x-direction
           dtb(i,1) = dt_sub(i)
           if((dx(i,1)<dxl(i,1)-eps) .or. (dx(i,1)>dxr(i,1)+eps) ) then
             if(v0curv(i,1)*vhcurv(i,1)<0d0) then
               dtb(i,1) = 0.5d0*dt_sub(i)
               vh(i,:) = v0(i,:)
               vhcurv(i,:) = v0curv(i,:)
             else
               if(it<5) then  !This branch needs to be taken out of the pcle loop!
                 vhcurv_inv(i,1) = 1d0/vhcurv(i,1)
                 dtb(i,1) = min(dt_sub(i),
     .              max((dxl(i,1)-eps)*vhcurv_inv(i,1)
     .                 ,(dxr(i,1)+eps)*vhcurv_inv(i,1)))
               else if((dx(i,1)<dxl(i,1)-eps2)
     .             .or. dx(i,1)>dxr(i,1)+eps2) then
                 dtb(i,1) = 0.5d0*dtb(i,1)
               end if
             end if
           end if

           !y-direction
           dtb(i,2) = dt_sub(i)
           if((dx(i,2)<dxl(i,2)-eps) .or. dx(i,2)>dxr(i,2)+eps) then
             if(v0curv(i,2)*vhcurv(i,2)<0d0) then
               dtb(i,2) = 0.5d0*dt_sub(i)
               vh(i,:) = v0(i,:)
               vhcurv(i,:) = v0curv(i,:)
             else 
               if(it<5) then  !This branch needs to be taken out of the pcle loop!
                 vhcurv_inv(i,2) = 1d0/vhcurv(i,2)
                 dtb(i,2) = min(dt_sub(i),
     .              max((dxl(i,2)-eps)*vhcurv_inv(i,2)
     .                 ,(dxr(i,2)+eps)*vhcurv_inv(i,2)))
               else if(dx(i,2)<dxl(i,2)-eps2 
     .            .or. dx(i,2)>dxr(i,2)+eps2) then
                 dtb(i,2) = 0.5d0*dtb(i,2)
               end if
             end if
           end if

           dt_sub(i) = min(dt_sub(i),dtb(i,1),dtb(i,2))
 
        end do

        dx(:,1) = dt_sub*vhcurv(:,1)
        dx(:,2) = dt_sub*vhcurv(:,2)
        dx(:,3) = dt_sub*vhcurv(:,3)

c     End program

      end subroutine Boris_solve

c     scatter_fields_nomap
c     ######################################################################
      subroutine scatter_fields_nomap(x0,dtau,dx,A_1,A_2,A_3
     $                               ,dtA_1,dtA_2,dtA_3,Ep,Bp)

      implicit none

c     Call variables

      real(_prec) :: x0(_Npg,3)
      real(_prec) :: dx(_Npg,3),Ep(_Npg,3),Bp(_Npg,3)
      real(_prec) :: dtau(_Npg)
      real(_prec) :: DtA_1 (_Npg,2,3),DtA_2 (_Npg,3,2),DtA_3 (_Npg,3,3)
      real(_prec) :: A_1   (_Npg,2,3),A_2   (_Npg,3,2),A_3   (_Npg,3,3)

c     Local variables

      integer     :: i,j
      real(_prec) :: Exp,Eyp,Ezp,Bxp,Byp,Bzp
     .              ,xhp(_Npg,3),xnp(_Npg,3),x0p(_Npg,3)
     .              ,x_nh(_Npg,3),cov(3,3),cnv(3,3),xlog(_Npg,3)
     .              ,dum1(_Npg),dum2(_Npg,2),dum3(_Npg,3),dt_nuh(_Npg)
      real(_prec) :: Sx2_n(_Npg,3), Sy2_n(_Npg,3)
      real(_prec) :: Sx1(_Npg,2),Sx2(_Npg,3),DSx(_Npg,3)
      real(_prec) :: Sy1(_Npg,2),Sy2(_Npg,3),DSy(_Npg,3)
      real(_prec) :: Sxd(_Npg,3),Syd(_Npg,3)
      real(_prec) :: Sx2_o(_Npg,3),Sy2_o(_Npg,3)

      real(_prec) :: Ax_nuh(_Npg,2,3),Ay_nuh(_Npg,3,2),Az_nuh(_Npg,3,3)

c     Begin program

      x0p(:,1) = x0(:,1)*ihx
      x0p(:,2) = x0(:,2)*ihy
      x0p(:,3) = x0(:,3)*ihz

      x_nh = x0 + dx*0.5d0
      xhp(:,1) = x_nh(:,1)*ihx
      xhp(:,2) = x_nh(:,2)*ihy      
      xhp(:,3) = x_nh(:,3)*ihz

      xnp  = x0 + dx
      xnp(:,1) = xnp(:,1)*ihx
      xnp(:,2) = xnp(:,2)*ihy
      xnp(:,3) = xnp(:,3)*ihz

      call shape2_offset(x0p(:,1),Sx2_o)
      call shape2_offset(x0p(:,2),Sy2_o)
      call shape2_offset(xnp(:,1),Sx2_n)
      call shape2_offset(xnp(:,2),Sy2_n)
      call shape1_offset(xhp(:,1),Sx1)
      call shape1_offset(xhp(:,2),Sy1)

      Sx2 = (Sx2_n+Sx2_o)*0.5d0
      Sy2 = (Sy2_n+Sy2_o)*0.5d0

      call Dx_shape2_offset(xhp(:,1),DSx)
      call Dy_shape2_offset(xhp(:,2),DSy)

      Sxd = (Sx2_n-Sx2_o)*0.5d0
      Syd = (Sy2_n-Sy2_o)*0.5d0

c     Find E

      call opt_matvec_23_pcle(E_1,Sy2,dum2)
      call dot_product_2_pcle(dum2,Sx1,Ep(:,1))
     
      call opt_matvec_2_pcle(Sx2,E_2,dum2)
      call dot_product_2_pcle(dum2,Sy1,Ep(:,2))

      call opt_matvec_tr_pcle(Sx2,E_3,dum3)
      call dot_product_3_pcle(dum3,Sy2,dum1)
      call opt_matvec_tr_pcle(Sxd,E_3,dum3)
      Ep(:,3) = dum1
      call dot_product_3_pcle(dum3,Syd,dum1)
      Ep(:,3) = Ep(:,3) + dum1

c     Find B

      !Prepare vector potential at nu+1/2
      dt_nuh = 0.5d0*(-2d0*dtp+dt+dtau)

      !VECTOR SIMD
      do i=1,_Npg
         Az_nuh(i,:,:) = A_3(i,:,:) + DtA_3(i,:,:)*dt_nuh(i)
      enddo

      call opt_matvec_tr_pcle(Sx2,Az_nuh,dum3)
      call dot_product_3_pcle(dum3,DSy,dum1)
      Bp(:,1) = dum1 + B0(1)

      call opt_matvec_tr_pcle(DSx,Az_nuh,dum3)
      call dot_product_3_pcle(dum3,Sy2,dum1)
      Bp(:,2) =-dum1 + B0(2)

      !VECTOR SIMD
      do i=1,_Npg
         Ay_nuh(i,:,:) = A_2(i,:,:) + DtA_2(i,:,:)*dt_nuh(i)
      enddo

      call opt_matvec_2_pcle (DSx,Ay_nuh,dum2)
      call dot_product_2_pcle(dum2,Sy1,dum1)
      Bp(:,3) = dum1

      !VECTOR SIMD
      do i=1,_Npg
         Ax_nuh(i,:,:) = A_1(i,:,:) + DtA_1(i,:,:)*dt_nuh1(i)
      enddo

      call opt_matvec_22_pcle(Sx1,Ax_nuh,dum3)
      call dot_product_3_pcle(dum3,DSy,dum1)
      Bp(:,3) = Bp(:,3) - dum1 + B0(3)

c     End program

      end subroutine scatter_fields_nomap

c     scatter_fields_curv
c     ######################################################################
      subroutine scatter_fields_curv(x0,dtau,dx,A_1,A_2,A_3
     $                              ,dtA_1,dtA_2,dtA_3,Ep,Bp)

      implicit none

c     Call variables

      real(_prec) :: x0(_Npg,3)
      real(_prec) :: dx(_Npg,3),Ep(_Npg,3),Bp(_Npg,3)
      real(_prec) :: dtau(_Npg)
      real(_prec) :: DtA_1 (_Npg,2,3),DtA_2 (_Npg,3,2),DtA_3 (_Npg,3,3)
      real(_prec) :: A_1   (_Npg,2,3),A_2   (_Npg,3,2),A_3   (_Npg,3,3)

c     Local variables

      integer     :: i
      real(_prec) :: xhp(_Npg,3),xnp(_Npg,3)
     .              ,x_nh(_Npg,3),cov(_Npg,3,3),cnv(_Npg,3,3)
     .              ,xlog(_Npg,3),x0p(_Npg,3)
     .              ,dum1(_Npg),dum2(_Npg,2),dum3(_Npg,3),dt_nuh(_Npg)
      real(_prec) :: Sx2_n(_Npg,3), Sy2_n(_Npg,3)
      real(_prec) :: Sx1(_Npg,2),Sx2(_Npg,3),DSx(_Npg,3)
      real(_prec) :: Sy1(_Npg,2),Sy2(_Npg,3),DSy(_Npg,3)
      real(_prec) :: Sxd(_Npg,3),Syd(_Npg,3)
      real(_prec) :: Sx2_o(_Npg,3),Sy2_o(_Npg,3)

      real(_prec) :: Ax_nuh(_Npg,2,3),Ay_nuh(_Npg,3,2),Az_nuh(_Npg,3,3)

c     Begin program

      x0p(:,1) = x0(:,1)*ihx
      x0p(:,2) = x0(:,2)*ihy
      x0p(:,3) = x0(:,3)*ihz

      x_nh = x0 + dx*0.5d0
      xhp(:,1) = x_nh(:,1)*ihx
      xhp(:,2) = x_nh(:,2)*ihy      
      xhp(:,3) = x_nh(:,3)*ihz

      xnp  = x0 + dx
      xnp(:,1) = xnp(:,1)*ihx
      xnp(:,2) = xnp(:,2)*ihy
      xnp(:,3) = xnp(:,3)*ihz

      call shape2_offset(x0p(:,1),Sx2_o)
      call shape2_offset(x0p(:,2),Sy2_o)
      call shape1_offset(xhp(:,1),Sx1)
      call shape1_offset(xhp(:,2),Sy1)
      call shape2_offset(xnp(:,1),Sx2_n)
      call shape2_offset(xnp(:,2),Sy2_n)
      Sx2 = (Sx2_n+Sx2_o)*0.5d0
      Sy2 = (Sy2_n+Sy2_o)*0.5d0
      call Dx_shape2_offset(xhp(:,1),DSx)
      call Dy_shape2_offset(xhp(:,2),DSy)

c     Find E

      Sxd = (Sx2_n-Sx2_o)*0.5d0
      Syd = (Sy2_n-Sy2_o)*0.5d0

      call opt_matvec_23_pcle(E_1,Sy2,dum2)
      call dot_product_2_pcle(dum2,Sx1,Ep(:,1))
     
      call opt_matvec_2_pcle(Sx2,E_2,dum2)
      call dot_product_2_pcle(dum2,Sy1,Ep(:,2))

      call opt_matvec_tr_pcle(Sx2,E_3,dum3)
      call dot_product_3_pcle(dum3,Sy2,dum1)
      call opt_matvec_tr_pcle(Sxd,E_3,dum3)
      Ep(:,3) = dum1
      call dot_product_3_pcle(dum3,Syd,dum1)
      Ep(:,3) = Ep(:,3) + dum1

c     Find B

      !Prepare vector potential at nu+1/2
      dt_nuh = 0.5d0*(-2d0*dtp+dt+dtau)

      !VECTOR SIMD
      do i=1,_Npg
         Az_nuh(i,:,:) = A_3(i,:,:) + DtA_3(i,:,:)*dt_nuh(i)
      enddo

      call opt_matvec_tr_pcle(Sx2,Az_nuh,dum3)
      call dot_product_3_pcle(dum3,DSy,dum1)
      Bp(:,1) = dum1

      call opt_matvec_tr_pcle(DSx,Az_nuh,dum3)
      call dot_product_3_pcle(dum3,Sy2,dum1)
      Bp(:,2) =-dum1

      !VECTOR SIMD
      do i=1,_Npg
         Ay_nuh(i,:,:) = A_2(i,:,:) + DtA_2(i,:,:)*dt_nuh(i)
      enddo

      call opt_matvec_2_pcle (DSx,Ay_nuh,dum2)
      call dot_product_2_pcle(dum2,Sy1,dum1)
      Bp(:,3) = dum1

      !VECTOR SIMD
      do i=1,_Npg
         Ax_nuh(i,:,:) = A_1(i,:,:) + DtA_1(i,:,:)*dt_nuh(i)
      enddo

      call opt_matvec_22_pcle(Sx1,Ax_nuh,dum3)
      call dot_product_3_pcle(dum3,DSy,dum1)
      Bp(:,3) = Bp(:,3) - dum1

c     Transform to Cartesian

      xlog(:,1) = x_nh(:,1)*ihx
      xlog(:,2) = x_nh(:,2)*ihy
      xlog(:,3) = 1d0
      call transf_vector_Bspl(i0,j0,xlog,cov,cnv)
      !VECTOR SIMD
      do i=1, _Npg
#ifdef __GNU
        Ep(i,:) = opt_matvec_tr(Ep(i,:),cov(i,:,:))
        Bp(i,:) = opt_matvec_tr(Bp(i,:),cnv(i,:,:))
#else        
        Ep(i,:) = matmul(Ep(i,:),cov(i,:,:))
        Bp(i,:) = matmul(Bp(i,:),cnv(i,:,:))
#endif
      end do

      Bp(:,1) = Bp(:,1) + B0(1)
      Bp(:,2) = Bp(:,2) + B0(2)
      Bp(:,3) = Bp(:,3) + B0(3)

c     End program

      end subroutine scatter_fields_curv

      end subroutine ACC_CN_incell_offset

      end module pic_mover

c pic_explicit_mover
c ######################################################################    
      module pic_explicit_mover

        use pic_io

      contains

ccc     Explicit_subcyc
ccc     ##################################################################
cc      subroutine Explicit_subcyc(pcle,i0,x0,v0,dt,E,phi,theta_in)
cc      implicit none
cc      
ccc     Call variables
cc
cc      type(particle), intent(inout) :: pcle
cc      real(8),intent(in) :: dt,E(0:),phi(0:)
cc      real(8),intent(in) :: x0(1),v0(1)
cc      integer,intent(in) :: i0
cc      real(8),intent(in),optional :: theta_in
cc
ccc     Local variables
cc
cc      integer :: it
cc      real(8) :: vnp(1),a0(1),am(1),x_k(1),anp(1)
cc      real(8) :: dx_pn,dx_pk
cc      real(8) :: beta,tt,dtt
cc      real(8) :: theta
cc
ccc     Begin program
ccc     Optional keyword
cc
cc      if(present(theta_in)) then
cc         theta = theta_in
cc      else
cc         theta = 5d-1
cc      end if
cc
cc
cc!     mid-point acceleration
cc      pcle%x_np = x0 + 0.5*dt*v0
cc      call par_bc_pcles(pcle%x_np(1),0)             
cc      pcle%i_np = locate(xx,pcle%x_np(1)) !careful         
cc      pcle%E_n = fld_on_pcleh(E,pcle)  !careful
cc      am = pcle%q_m*  pcle%E_n
cc
cc      pcle%v_np = v0 + dt*am
cc
cc      pcle%x_nh = pcle%x_np
cc      pcle%i_nh = pcle%i_np
cc      pcle%v_nh = 0.5*(v0+pcle%v_np)
cc
cc      pcle%x_np = pcle%x_np + 0.5*dt*pcle%v_np
cc      call par_bc_pcles(pcle%x_np(1),0)             
cc      pcle%i_np = locate(xx,pcle%x_np(1)) !careful         
cc
ccc     End program
cc
cc      end subroutine Explicit_subcyc
cc
ccc     move_pcles_ex_subcyc
ccc     ##################################################################
cc        subroutine move_pcles_ex_subcyc(nx,spcs,dt,E,phi,js)
cc        
cc        implicit none
cc        
ccc     Call variables
cc        integer :: nx
cc        type(species) :: spcs(:)
cc        real(8) :: dt,E(0:),phi(0:)
cc        real(8), intent(out) :: js(0:nx+1)
cc
ccc     Local variables
cc        type(particle),pointer,dimension(:) :: pcles
cc        integer :: ipc,isp,nt_sub,isub,i0
cc        real(8) :: h,ih,dt_sub,x0(1),v0(1)
cc        
cc        real(8) :: garr(0:nx)
ccc     Begin program
cc        h = xx(2) - xx(1)
cc        ih = 1d0/h
cc        js = 0d0
cc
ccc     Move particles 
cc        do isp=1,size(spcs)
cc           pcles => spcs(isp)%pcles
cc           garr = 0d0
cc           nt_sub = nt_suba(isp)
cc           dt_sub = dt/nt_sub
cc           do ipc=1,size(pcles)
cc              i0 = pcles(ipc)%i_n
cc              x0 = pcles(ipc)%x_n
cc              v0 = pcles(ipc)%v_n
cc
cc              call Explicit_subcyc(
cc     $             pcles(ipc),i0,x0,v0,dt_sub,E,phi)
cc              call collect_j(pcles(ipc),i0,x0,v0,dt_sub)
cc              
cc              do isub = 1,nt_sub-1,1
cc                 i0 = pcles(ipc)%i_np
cc                 x0 = pcles(ipc)%x_np
cc                 v0 = pcles(ipc)%v_np                
cc                 call Explicit_subcyc(
cc     $                pcles(ipc) ,i0,x0,v0,dt_sub,E,phi)
cc                 call collect_j(pcles(ipc),i0,x0,v0,dt_sub)
cc              end do
cc
cc           end do
cc           garr = garr*ih/dt    !careful, to be dimensionally correct
cc!     BCs
cc           garr(0)  = garr(0) + garr(nx)
cc           garr(nx) = garr(0 )
cc           js(0:nx) = js(0:nx) + spcs(isp)%q*garr(0:nx)
cc         
cc        end do
cc
cc        CONTAINS
cc              
cc
cc
cc        subroutine collect_j(pcle,i0,x0,v0,factor)
cc        type(particle), intent(in) :: pcle
cc        integer :: i0
cc        real(8) :: x0(1),v0(1)
cc        integer :: i
cc        real(8) ::vh,xp,factor
cc
cc        vh = (v0(1) + pcle%v_np(1)) * 0.5d0 * factor
cc        i  = pcle%i_nh
cc        xp = pcle%x_nh(1)
cc        call garrp(i,xp,vh)
cc
cc        end subroutine collect_j
cc
cc      subroutine garrp(i,xp,vh)
ccc     Call variables
cc      integer :: i
cc      real(8) :: xp,vh
ccc     Local variables
cc      integer :: is,ig,ileft
cc      real(8) :: xxh0,xxh
cc      
cc
cc      xxh0 = (xx(i) + xx(i+1))*0.5d0 
cc      xxh = xxh0 - (sp_order/2)*h
cc      do is = 1, sp_order
cc         ig  = cellh_sp_index(i,is)
cc         garr(ig) = garr(ig) + pshapej((xxh-xp)*ih)*vh
cc         xxh = xxh + h
cc      end do
cc      
cc      end subroutine garrp
cc
cc
ccc     End program
cc
cc      end subroutine move_pcles_ex_subcyc

cc     LF_v_udpate
cc     ##################################################################
c      subroutine LF_v_update (pcle,i0,j0,x0,v0,dt,DtA,B0,A,phi,dx,ipc)
c
c      implicit none
c      
cc     Call variables
c
c      type(particle), intent(inout) :: pcle
c      real(8),intent(in) :: A(0:,0:,:),phi(0:,0:)
c      real(8),intent(in) :: dt,DtA(0:,0:,:),B0(3)
c      real(8),intent(out):: dx(_Dim)
c      real(8),intent(in) :: x0(_Dim),v0(_Dim)
c      integer,intent(in) :: i0,ipc,j0
c
cc     Local variables
c      real(8) :: Ep(3),Bp(3),vh(3),vhes(3),vnp(3),xnp(3),vdotB, qmdth
c     $     ,iDenom
c      integer :: inp,jnp
c
cc     Begin program
c
cc     from last time step
c
c      pcle%x_np = x0
c      pcle%v_np = v0 
c
c      pcle%i_np = i0    
c      pcle%j_np = j0
c
cc     Estimate sub-time-step 
c
c      dx = 0d0
c      
c      call scatter_fields_ex(x0,dx,Ep,Bp)
c
cc     Velocity update
c
c      qmdth = pcle%q_m*dt*0.5d0
c      vhes = v0 + qmdth*Ep
c
c      iDenom = 1d0/(1d0 + qmdth*qmdth*dot_product(Bp,Bp))
c      vdotB = dot_product(vhes,Bp)
c      vh(1) = (vhes(1) + qmdth*(vhes(2)*Bp(3)-vhes(3)*Bp(2)
c     .     + qmdth*vdotB*Bp(1)))*iDenom
c      vh(2) = (vhes(2) + qmdth*(vhes(3)*Bp(1)-vhes(1)*Bp(3) 
c     .     + qmdth*vdotB*Bp(2)))*iDenom
c      vh(3) = (vhes(3) + qmdth*(vhes(1)*Bp(2)-vhes(2)*Bp(1) 
c     .     + qmdth*vdotB*Bp(3)))*iDenom
c
c      vnp = v0 + 2d0*(vh - v0) !2d0*vh - v0
c
cc$$$      dx  = dt_sub*vh
cc$$$
cc$$$!     Position update (2 steps)
cc$$$      xnp = x0 + dx*0.5d0
cc$$$      inp = i0
cc$$$      jnp = j0
cc$$$      call cell_crossing(dx,xnp,vnp,inp,jnp)
cc$$$
cc$$$      pcle%x_nh = xnp
cc$$$      pcle%i_nh = inp
cc$$$      pcle%j_nh = jnp
cc$$$
cc$$$!     2nd half step
cc$$$      xnp = xnp + dx*0.5d0
cc$$$      call cell_crossing(dx,xnp,vnp,inp,jnp)
cc$$$
cc$$$      pcle%x_np = xnp
cc$$$      pcle%i_np = inp
cc$$$      pcle%j_np = jnp      
cc$$$
c      pcle%v_np = vnp
c
cc     End program
c
c      contains
c
cc     scatter_fields_ex
cc     ######################################################################
c      subroutine scatter_fields_ex(x0,dx,Ep,Bp)
c
c      implicit none
c
cc     Call variables
c
c      real(8) :: x0(2),dx(2),Ep(3),Bp(3)
c
cc     Local variables
c
c      real(8) :: Exp,Eyp,Ezp,Bxp,Byp,Bzp ,x0p,y0p,xhp,yhp,xnp,ynp,x_nh
c     $     ,y_nh ,Sy1(2),Sy2(3),DSy(3),Sxd(3),Syd(3),Sx1(2),Sx2(3)
c     $     ,DSx(3),cov(3,3),cnv(3,3),Ax_nuh(2,3),Az_nuh(3,3),Ay_nuh(3,2)
c     .     ,xlog(3)
c
c      real(8) :: xp,xp2,E0,E1,E2,Bz1,Bz2,By1,By2,Ey0,Ey1,Ey2,Ez0,Ez1,Ez2
c
cc     Begin program
c
c      x0p = x0(1)/hx
c      y0p = x0(2)/hy
c
c      x_nh = x0(1) + dx(1)*0.5d0
c      xhp = x_nh/hx
c
c      y_nh = x0(2) + dx(2)*0.5d0
c      yhp = y_nh/hy      
c
c      xnp = (x0(1) + dx(1))/hx
c      ynp = (x0(2) + dx(2))/hy
c
c      Sx1 = shape1_offset(xhp)
c      Sy1 = shape1_offset(yhp)
c      Sx2 = (shape2_offset(xnp)+shape2_offset(x0p))*0.5d0
c      Sy2 = (shape2_offset(ynp)+shape2_offset(y0p))*0.5d0
c      DSx = Dx_shape2_offset(xhp)
c      DSy = Dy_shape2_offset(yhp)
c
c      !Find E
c      Exp   =-dot_product(MATMUL(DSx,phi(i0-1:i0+1,j0-1:j0+1)),Sy2) 
c     $       -dot_product(MATMUL(DtA(i0-1:i0,j0-1:j0+1,1),Sy2),Sx1)
c
c      Eyp = -dot_product(MATMUL(Sx2,DtA(i0-1:i0+1,j0-1:j0,2)),Sy1)
c     .      -dot_product(MATMUL(Sx2,phi(i0-1:i0+1,j0-1:j0+1)),DSy)
c
c      Ezp = -dot_product(MATMUL(Sx2,DtA(i0-1:i0+1,j0-1:j0+1,3)),Sy2)
c
c      Sxd = (shape2_offset(xnp)-shape2_offset(x0p))*0.5d0
c      Syd = (shape2_offset(ynp)-shape2_offset(y0p))*0.5d0
c      Ezp = Ezp-dot_product(MATMUL(Sxd,DtA(i0-1:i0+1,j0-1:j0+1,3)),Syd)
c
c      Ep = (/Exp,Eyp,Ezp/)
c
c      !Find B
c      Az_nuh =  A(i0-1:i0+1,j0-1:j0+1,3)
c      Ay_nuh =  A(i0-1:i0+1,j0-1:j0  ,2)
c      Ax_nuh =  A(i0-1:i0  ,j0-1:j0+1,1)
c
c      Bxp = dot_product(MATMUL(Sx2,Az_nuh),DSy)
c      Byp =-dot_product(MATMUL(DSx,Az_nuh),Sy2)
c      Bzp = dot_product(MATMUL(DSx,Ay_nuh),Sy1)
c     $     -dot_product(MATMUL(Sx1,Ax_nuh),DSy)
c
c      Bp = (/Bxp,Byp,Bzp/)
c
c      !Transform to Cartesian
c      if (no_map) then
c        Bp = Bp + B0
c      else
c        xlog(1) = pxx(i0)       !????
c        xlog(2) = pyy(j0)       !????
c        xlog(3) = pzz(1)      
c        cov = cov_vector_anal(xlog)
c        cnv = cnv_vector_anal(xlog)
ccc        xlog(1) = (x_nh)/hx
ccc        xlog(2) = (y_nh)/hy
ccc        xlog(3) = 1d0
ccc        call transf_vector_Bspl(i0,j0,xlog,cov,cnv=cnv)
c        Ep = matmul(Ep,cov)
c        Bp = matmul(Bp,cnv) + B0
c      end if
c
cc     End program
c
c      end subroutine scatter_fields_ex      
c
c      end subroutine LF_v_update
c
cc     LF_x_update
cc     ##################################################################
c      subroutine LF_x_update (pcle,i0,j0,x0,v0,dt,dtp,dt_sub
c     $     ,dx,ipc)
c
c      implicit none
c      
cc     Call variables
c
c      type(particle), intent(inout) :: pcle
c      real(8),intent(in) :: dt,dtp
c      real(8),intent(out):: dt_sub,dx(_Dim)
c      real(8),intent(in) :: x0(_Dim),v0(_Dim)
c      integer,intent(in) :: i0,j0,ipc
c
cc     Local variables
c
c      integer :: in(3),in0(3),nsize(3),dim
c      real(8), parameter ::eps= 1d-14
c      real(8) :: dxl(2),dxr(2),xp(3)
c
cc     Begin program
c
cc     from last time step
c
c      pcle%x_np = x0
c      pcle%i_np = i0    
c      pcle%j_np = j0
c
cc     Estimate sub-time-step 
c
c      dt_sub = dtp
c
c      dxl(1) = -x0(1)
c      dxr(1) = hx-x0(1)      
c      dxl(2) = -x0(2)
c      dxr(2) = hy-x0(2)
c
c      dt_sub = min(dt_sub,max((dxl(1)-eps)/v0(1),(dxr(1)+eps)/v0(1)))
c      dt_sub = min(dt_sub,max((dxl(2)-eps)/v0(2),(dxr(2)+eps)/v0(2)))
c
c      dx  = dt_sub*v0
c      pcle%x_nh = x0 + dx*0.5d0
c      pcle%x_np = x0 + dx
c      pcle%i_nh = i0                                                 
c      pcle%i_np = i0                                    
c
c      pcle%j_nh = j0
c      pcle%j_np = j0 
c
cc     Cell crossing and particle b.c.
c
c      call cell_crossing(dx,pcle%x_np,pcle%v_np,pcle%i_np,pcle%j_np)
c
cc     End program
c
c      end subroutine LF_x_update

cc     init_vp_exp
cc     ##################################################################
c      subroutine init_vp_exp(dt,Ex,Ey,Ez,By,Bz,mments)
c
c      implicit none
c      
cc     Call variables
c
c      real(8), intent(in)  :: dt,Ex(0:), Ey(0:), Ez(0:), By(0:), Bz(0:)
c      real(8), intent(out) :: mments(0:,0:,:,0:)
c
cc     Local variables 
c
c      type(particle),pointer,dimension(:) :: pcles
c      integer :: ipc,isp, i0,inp,inh
c      real(8) :: xp,E0,E1,E2,Exp,Bz1,Bz2,Bzp,By1,By2,Byp,Bxp
c     .     ,xp2,Ey0,Ey1,Ey2,Eyp,Ez0,Ez1,Ez2,Ezp,vh(3),vhes(3)
c     .     ,qmdth,Denom,vdotB,x0(3),dx(3),xnp(3),v0(3),vnp(3)
c
c
c      integer :: id_mp, nth_mp,id
c
c      real(8) :: garr_n  (0:nxg+1,0:nyg+1)
c     $          ,garr_vx (0:nxg+1,0:nyg+1)
c     $          ,garr_vy (0:nxg+1,0:nyg+1)
c     $          ,garr_vz (0:nxg+1,0:nyg+1)
c     $          ,garr_vx2(0:nxg+1,0:nyg+1)
c     $          ,garr_vy2(0:nxg+1,0:nyg+1)
c     $          ,garr_vz2(0:nxg+1,0:nyg+1)
c
c      real(8) :: My2(1,3),Mx2(3,1),My1(1,2),Mx1(2,1),Mxy2(3,3)
c      real(8) :: dtp,dt_sub,yp,v2
c      integer :: j0,i,j,it
c      real(8) :: dth
c
c      dth = dt*0.5d0
c      nth_mp = 1
c      id_mp  = 1
c      mments = 0d0
c
c!$omp parallel private(nth_mp)
c!$      nth_mp = omp_get_num_threads()
cc      id_mp  = omp_get_thread_num()
c!$omp end parallel                                                             
c
cc     Move particles 
c        do isp=1, size(spcs)
c           pcles => spcs(isp)%pcles
c           garr_n  = 0d0
c           garr_vx = 0d0
c           garr_vy = 0d0
c           garr_vz = 0d0
c           garr_vx2= 0d0
c           garr_vy2= 0d0
c           garr_vz2= 0d0
c!$omp parallel DEFAULT(SHARED) private(ipc,i0,inp
c!$omp.     ,xp,E0,E1,E2,Exp,Bz1,Bz2,Bzp,By1,By2,Byp,Bxp
c!$omp.     ,xp2,Ey0,Ey1,Ey2,Eyp,Ez0,Ez1,Ez2,Ezp,vh,vhes
c!$omp.     ,qmdth,Denom,vdotB,x0,dx,xnp,v0,vnp,id_mp)
c
c           id_mp  = 1
c!$         id_mp  = omp_get_thread_num() + id_mp  !careful
c
c!$           call set_omp_thread_id()
c
c!$omp DO 
c           do ipc=nip+1,size(pcles)
cc     Boris's push
c              v0 = pcles(ipc)%v_n
c              x0 = pcles(ipc)%x_n
c              
c
c              xp = x0(1)/hx
c              xp2 = xp - 0.5d0
c              i0 = pcles(ipc)%i_n
c
c              E0 = Ex(i0-1)
c              E1 = Ex(i0)
c              E2 = Ex(i0+1)
c              Exp = (E0*(1d0-xp) + E1*xp) 
cc$$$              Exp = E0*(0.5d0*(0.5d0-xp2)*(0.5d0-xp2)) +
cc$$$     .              E1*(0.75d0-xp2*xp2) +
cc$$$     .              E2*(0.5d0*(0.5d0+xp2)*(0.5d0+xp2))
c
c
c              Bz1 = Bz(i0-1) 
c              Bz2 = Bz(i0)   
c              Bzp = Bz1*(1d0-xp) + Bz2*xp
c
c              By1 = By(i0-1) 
c              By2 = By(i0)   
c              Byp = By1*(1d0-xp) + By2*xp
c
c
c
c              Ey0 = Ey(i0-1)
c              Ey1 = Ey(i0)
c              Ey2 = Ey(i0+1)
c
c              Eyp = Ey0*(0.5d0*(0.5d0-xp2)*(0.5d0-xp2)) +
c     .              Ey1*(0.75d0-xp2*xp2) +
c     .              Ey2*(0.5d0*(0.5d0+xp2)*(0.5d0+xp2))
c
c
c              Ez0 = Ez(i0-1)
c              Ez1 = Ez(i0)
c              Ez2 = Ez(i0+1)
c
c              Ezp = Ez0*(0.5d0*(0.5d0-xp2)*(0.5d0-xp2)) +
c     .              Ez1*(0.75d0-xp2*xp2) +
c     .              Ez2*(0.5d0*(0.5d0+xp2)*(0.5d0+xp2))
cc              Ezp = Ez0*(1d0-xp) + Ez1*xp
c
c              Bxp = 0d0
c
c              qmdth = pcles(ipc)%q_m*dth*0.25d0 !dt/4
c
c
c              vhes(1) = v0(1) + qmdth*Exp
c              vhes(2) = v0(2) + qmdth*Eyp
c              vhes(3) = v0(3) + qmdth*Ezp
c              
c              Denom = 1d0/(1d0 + qmdth*qmdth*(Bxp*Bxp+Byp*Byp+Bzp*Bzp))
c              
c              vdotB = vhes(1)*Bxp + vhes(2)*Byp + vhes(3)*Bzp
c              vh(1) = (vhes(1) + qmdth*(vhes(2)*Bzp-vhes(3)*Byp
c     .             + qmdth*vdotB*Bxp))*Denom
c              vh(2) = (vhes(2) + qmdth*(vhes(3)*Bxp-vhes(1)*Bzp 
c     .             + qmdth*vdotB*Byp))*Denom
c              vh(3) = (vhes(3) + qmdth*(vhes(1)*Byp-vhes(2)*Bxp 
c     .             + qmdth*vdotB*Bzp))*Denom
c              
cc     update velocity
cc              vnp = 2d0*vh - v0
c              vnp = v0 + 2d0*(vh - v0)
c             
c      pcles(ipc)%v_np    = vnp
c
cc     update position
c              i0 = pcles(ipc)%i_n
c              j0 = pcles(ipc)%j_n
c              x0 = pcles(ipc)%x_n
c              v0 = pcles(ipc)%v_n   
c              it = 0
c              dtp = dth
c              do while(dtp>0)
c                 call LF_x_update(pcles(ipc),i0,j0,x0,vnp,dth
c     $                ,dtp,dt_sub,dx,ipc)
c               xp = pcles(ipc)%x_nh(1)/hx 
c               yp = pcles(ipc)%x_nh(2)/hy  
c               
c               Mx1(1:2,1) = shape1_offset(xp)
c               My1(1,1:2) = shape1_offset(yp)               
cc               Mx2(1:3,1) = shape2_offset(xp)
cc               My2(1,1:3) = shape2_offset(yp)
c
c               Mx2(1:3,1) = (-shape2_offset(x0(1)/hx)+
c     $              shape2_offset((x0(1)+dx(1))/hx))*0.5d0 !careful
c               My2(1,1:3) = (-shape2_offset(x0(2)/hy)+
c     $              shape2_offset((x0(2)+dx(2))/hy))*0.5d0
c               Mxy2(1:3,1:3)= matmul(Mx2,My2)
c
c               Mx2(1:3,1) = (shape2_offset(x0(1)/hx)+
c     $              shape2_offset((x0(1)+dx(1))/hx))*0.5d0 !careful
c               My2(1,1:3) = (shape2_offset(x0(2)/hy)+
c     $              shape2_offset((x0(2)+dx(2))/hy))*0.5d0
c
c               Mxy2(1:3,1:3)= Mxy2(1:3,1:3) + matmul(Mx2,My2)
c
c               i  = pcles(ipc)%i_nh
c               j  = pcles(ipc)%j_nh ! j0
c               
c!     accumulate vx (linear in x, 2nd in y)
c               garr_vx(i-1:i  ,j-1:j+1) = 
c     .         garr_vx(i-1:i  ,j-1:j+1) + dx(1)*matmul(Mx1,My2)
c     .                                     *pcles(ipc)%w
c!     accumulate vy (linear in y, 2nd in x)
c               garr_vy(i-1:i+1,j-1:j  ) = 
c     .         garr_vy(i-1:i+1,j-1:j  ) + dx(2)*matmul(Mx2,My1)               
c     .                                     *pcles(ipc)%w
c
!!     accumulate vz (2nd in x,y)
c               garr_vz(i-1:i+1,j-1:j+1) =
c     .         garr_vz(i-1:i+1,j-1:j+1) + dx(3)*Mxy2                      
c     .                                     *pcles(ipc)%w
c
c
c               it = it+1
c
c               i0 = pcles(ipc)%i_np
c               j0 = pcles(ipc)%j_np
c               x0 = pcles(ipc)%x_np
c               v0 = pcles(ipc)%v_np
c               dtp = dtp-dt_sub
c
c            end do              !sub-stepping
c!     update particle position to 1/2
c            pcles(ipc)%x_n     = x0 !careful
c            pcles(ipc)%i_n     = i0
c
c!     2nd order (for checking charge conservation) careful
c            xp = x0(1)/hx
c            yp = x0(2)/hy
c
c            Mx1(1:2,1) = shape1_offset(xp)
c            My1(1,1:2) = shape1_offset(yp)               
c            Mx2(1:3,1) = shape2_offset(xp)
c            My2(1,1:3) = shape2_offset(yp)
c
c            Mxy2(1:3,1:3)= matmul(Mx2,My2)
c
c!     accumulate n (quadratic in x,y)
c            garr_n  (i0-1:i0+1,j0-1:j0+1) = 
c     .      garr_n  (i0-1:i0+1,j0-1:j0+1) +    Mxy2*pcles(ipc)%w
c
c!     accumulate v^2 (2nd in x,y, cell-centered)
c            v2 = pcles(ipc)%v_np(1)*pcles(ipc)%v_np(1)
c            garr_vx2(i0-1:i0+1,j0-1:j0+1) = 
c     .      garr_vx2(i0-1:i0+1,j0-1:j0+1) + v2*Mxy2*pcles(ipc)%w
c
c            v2 = pcles(ipc)%v_np(2)*pcles(ipc)%v_np(2)
c            garr_vy2(i0-1:i0+1,j0-1:j0+1) = 
c     .      garr_vy2(i0-1:i0+1,j0-1:j0+1) + v2*Mxy2*pcles(ipc)%w   
c
c            v2 = pcles(ipc)%v_np(3)*pcles(ipc)%v_np(3)
c            garr_vz2(i0-1:i0+1,j0-1:j0+1) = 
c     .      garr_vz2(i0-1:i0+1,j0-1:j0+1) + v2*Mxy2*pcles(ipc)%w
c              
c           end do
c!$omp end do 
c!$omp end parallel 
c
cc     collect moments
c
c!     orbit-average
c         mments(:,:,isp,Ipx) = garr_vx/dth
c         mments(:,:,isp,Ipy) = garr_vy/dth
c         mments(:,:,isp,Ipz) = garr_vz/dth
c!     non-orbit-average
c         mments(:,:,isp,Ins)  = garr_n !careful        
c         mments(:,:,isp,IS11) = garr_vx2
c         mments(:,:,isp,IS22) = garr_vy2
c         mments(:,:,isp,IS33) = garr_vz2
c        
c        end do !species
c
c        call wrap_moments_nodes(mments)
c
cc     End program
c
c      end subroutine init_vp_exp

c     update_xp
c     ##################################################################
      subroutine update_xp

c     ------------------------------------------------------------------
c     Save old positions and velocities
c     ------------------------------------------------------------------
      
      implicit none

c     Call variables

c     Local variables

      integer :: isp,ipc,numpcles,numspcs

c     Begin program
     
      numspcs = size(spcs)
      do isp=1,numspcs
        numpcles = size(spcs(isp)%pcles)
!$OMP PARALLEL DO PRIVATE(ipc)
        do ipc=nip+1,numpcles
          spcs(isp)%pcles(ipc)%x_n = spcs(isp)%pcles(ipc)%x_np
          spcs(isp)%pcles(ipc)%ijk_n = spcs(isp)%pcles(ipc)%ijk_np              
        end do
!$OMP END PARALLEL DO
      end do

c     End program

      end subroutine update_xp

      end module pic_explicit_mover

c pic
c ######################################################################    
      module pic 

      use pic_moments

      use pic_mass_matrix

      use pic_setup

      use pic_explicit_mover

      use pic_sort

      use pic_mover

      integer, parameter :: NDim = _Dim

      contains

c     pic_update_tnp
c     ##################################################################
      subroutine pic_update_tnp(A,DtA,test,init)

c     ------------------------------------------------------------------
c     Update particle position and velocities in time, compute inflow
c     np filling.
c     ------------------------------------------------------------------

      implicit none

c     Call variables

      real(8), intent(in) :: A(0:,0:,1:),DtA(0:,0:,1:)

      character(*) :: test

      logical ::  init

c     Local variables

      integer :: isp,ipc,Np,npc1,ierr !count the # of particle in the first cell
      integer :: numspcs,nxg,nyg,nzg,nng,npl
      real(8) :: cm(3),ihx,ihy,ihz
      real(_prec) :: v2
      real(_prec) :: xp(_Npg), yp(_Npg)
      integer :: i0(_Npg), j0(_Npg)
      integer :: i, i_n(_Npg),j_n(_Npg)

      integer,pointer,dimension(:,:,:,:) :: Npcle_cell

c     Begin program

c     Update moments

      momsn = moms

cc      if (init) return
      if (.not.init) then

c     Update particle info

      call get_inverse_h(ihx,ihy,ihz)
      call get_global_grid_size(nxg,nyg,nzg)
      nng = nxg*nyg*nzg

      numspcs = size(spcs)

      select case(test)
      case('slow_shock')
        npc1 = 0

        if (particle_sort) then
          allocate(Npcle_cell(nxg,nyg,nzg,numspcs))
          call sort_particles(Npart_cell=Npcle_cell)
          isp = 2 !Electrons
          npc1 = Npcle_cell(1,1,1,isp)+Npcle_cell(2,1,1,isp) !Sum first two cells
          deallocate(Npcle_cell)
        else
          do isp=1,numspcs
            Np = size(spcs(isp)%pcles)
!$OMP PARALLEL DEFAULT(SHARED) private(ipc,i,i_n,j_n)
!$OMP DO 
!$OMP.REDUCTION (+:npc1)
            do ipc=nip+1,Np
              !Save old positions and velocities
              spcs(isp)%pcles(ipc)%x_n = spcs(isp)%pcles(ipc)%x_np
              spcs(isp)%pcles(ipc)%v_n = spcs(isp)%pcles(ipc)%v_np
              spcs(isp)%pcles(ipc)%w_n = spcs(isp)%pcles(ipc)%w_np
              spcs(isp)%pcles(ipc)%ijk_n = spcs(isp)%pcles(ipc)%ijk_np    

              call xform_pcle_idx(spcs(isp)%pcles(ipc)%ijk_n
     .                           ,nxg,i_n,j_n)
              do i=1, _Npg
                if((i_n(i).eq.1.or.i_n(i).eq.2).and.(isp.eq.2)) then
                  npc1 = npc1 + 1
                end if
              end do
            end do
!$OMP END DO
!$OMP END PARALLEL 
          end do
        endif

        !Calculate inflow np filling
        if(nip>0) then
          !keep the first cell particle number constant
          npl = (npc1/_Npg-npc1_0(1))
           
          do while(npl<0) 
            nip = nip - 1
            npl = npl + 1
          end do

          if (my_rank_pic == 0.and.(.not.tst_flg)) then
            write (*,*) 'Particle buffer number',nip
          endif

          cm = 0d0
          do isp=1,numspcs
            v2 = 0            
!$OMP PARALLEL DEFAULT(SHARED) private(ipc,xp,yp,i0,j0,i)
!$OMP DO 
!$OMP.REDUCTION (+:v2,cm)
            do ipc = nip+1,nip+npa_int  !What is npa_int for????
              xp = spcs(isp)%pcles(ipc)%x_np(:,1)*ihx 
              yp = spcs(isp)%pcles(ipc)%x_np(:,2)*ihy
              call xform_pcle_idx(spcs(isp)%pcles(ipc)%ijk_np,nxg,i0,j0)
              do i=1, _Npg
                 v2 = v2
     .              +( spcs(isp)%pcles(ipc)%v_np(i,1)
     .                *spcs(isp)%pcles(ipc)%v_np(i,1)
     .              +  spcs(isp)%pcles(ipc)%v_np(i,2)
     .                *spcs(isp)%pcles(ipc)%v_np(i,2)
     .              +  spcs(isp)%pcles(ipc)%v_np(i,3)
     .                *spcs(isp)%pcles(ipc)%v_np(i,3)
     .               )*spcs(isp)%pcles(ipc)%w_np(i)
#ifdef PIC_CM_DIAG
                 cm = cm + cmp(spcs(isp)%q,i0(i),j0(i),xp(i),yp(i),A)
     .            + 0.5*dt*cmp(spcs(isp)%q,i0(i),j0(i),xp(i),yp(i),DtA)
                 cm(1)= cm(1)+spcs(isp)%m*spcs(isp)%pcles(ipc)%v_np(i,1)
                 cm(2)= cm(2)+spcs(isp)%m*spcs(isp)%pcles(ipc)%v_np(i,2)
                 cm(3)= cm(3)+spcs(isp)%m*spcs(isp)%pcles(ipc)%v_np(i,3)
#endif
              end do
            end do
!$OMP END DO
!$OMP END PARALLEL 
            m_np(Ien) = m_np(Ien) + 0.5d0*spcs(isp)%m*v2
          end do
        end if

        m_np(Icmx:Icmz) = m_np(Icmx:Icmz) + cm
 
        if(nip < 0) then
          print *, "no more inactive particles, quit."
          stop
        end if

      case default

        if (particle_sort) then
          call sort_particles
        else
          do isp=1,numspcs
            Np = size(spcs(isp)%pcles)
!$omp parallel DEFAULT(SHARED) private(ipc)
            do ipc=nip+1,Np
              spcs(isp)%pcles(ipc)%x_n = spcs(isp)%pcles(ipc)%x_np
              spcs(isp)%pcles(ipc)%v_n = spcs(isp)%pcles(ipc)%v_np
              spcs(isp)%pcles(ipc)%ijk_n = spcs(isp)%pcles(ipc)%ijk_np
              spcs(isp)%pcles(ipc)%w_n = spcs(isp)%pcles(ipc)%w_np
            end do
!$omp end parallel 
          end do
        endif

      end select

      endif

c     Dump particles for restart

      if ((TS_chk_dump().or.itime==0).and.(.not.tst_flg)) then
        ierr = pic_dump_pcles(time)
      endif

c     End program

      end subroutine pic_update_tnp

c$$$c     find_cb
c$$$c     ##################################################################
c$$$      function find_cb(xx,i0,in) result(xcb)
c$$$
c$$$      implicit none
c$$$
c$$$c     Call variables
c$$$
c$$$      real(8) :: xx(0:)
c$$$      integer, intent(in)  :: i0,in
c$$$      real(8) :: xcb
c$$$      
c$$$c     Begin program
c$$$      xcb = 0.5d0*(xx(i0) + xx(in))
c$$$        
c$$$c     End program
c$$$
c$$$      end function find_cb
c$$$
c$$$c     f0v
c$$$c     ##################################################################
c$$$      function f0v(v) result(f0)
c$$$
c$$$      implicit none
c$$$
c$$$c     Call variables
c$$$
c$$$      real(8) , intent(in)  :: v
c$$$      real(8) :: f0
c$$$      
c$$$c     Local variables
c$$$
c$$$      real(8) , parameter :: isq2pi  = 0.3989422804014327
c$$$
c$$$c     Begin program
c$$$
c$$$      f0 = exp(-v*v*0.5)*isq2pi
c$$$
c$$$c     End program
c$$$
c$$$      end function f0v

      end module pic
#endif
