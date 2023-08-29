c     module EM_utils
c     ##################################################################
      module EM_utils

      use gm_init

      INTEGER :: MAXSPLIT
CGMM_OMP!$OMP THREADPRIVATE(MAXSPLIT)

C     Constants that define array sizes at compilation time.

      INTEGER,PARAMETER :: MITERMAX=30000
C         maximum number of iterations for the EM algorithm

      INTEGER,PARAMETER :: DENMAX=175
C          maximum value of the A term in exp(-A) used when calculating
C          the density of a point. Above this value exp(-A) is equated
C          to zero.

      CONTAINS

c     CAPART
c     ##################################################################
      SUBROUTINE CAPART(NIND,NATT,NG,W,IDT,XCC)
C     This subroutine determines  the partition of entities,
C     from the posterior probabilities W, into NG groups
      implicit double precision (a-h,o-z)
      DIMENSION W(NIND,NG),IDT(NIND),XCC(NG)

      XCC=0.0

      if (NIND == 0) return

      DO I=1,NIND
        MAX=1

        DO K=2,NG
          IF (W(I,K).GT.W(I,MAX)) MAX=K
        ENDDO
        
        XCC(MAX)=XCC(MAX)+W(I,MAX)

        IDT(I)=MAX

        CKOUT= SUM(W(I,1:NG))

        IF (CKOUT.LT.0.0001) IDT(I)=0
      END DO
      
      END SUBROUTINE CAPART
      
C     ENFORCE_MOMS_CONSERVATION
C     ######################################################################
      SUBROUTINE ENFORCE_MOMS_CONSERVATION(NG,NATT,TS,XMUS,XVARS
     $                  ,T,XMU,XVAR)

C     Input: TS = (should be 1): sample's proportion
C	     XMUS = sample's mean
C	     XVARS = sample's covariance
C     Output: Conserved quantities T,XMU,XVAR
      IMPLICIT NONE

C     Call variables

      INTEGER, INTENT(IN) :: NG,NATT
      REAL(8), INTENT(IN) :: TS,XMUS(NATT),XVARS(NATT,NATT)
      REAL(8), INTENT(INOUT) ::T(NG),XMU(NG,NATT)
     $          ,XVAR(NG,NATT,NATT)

C     Local variables

      REAL(8) :: XMU_D(NATT),XVAR_D(NATT,NATT)
      INTEGER :: I,J,K
      REAL(8) :: DIFF,TOL

C     Conserve proportion (mass) - zeroth-order mom

      IF ( NG .EQ. 1 ) THEN
        T(NG) = 1D0
      ELSE
        T(NG) = 1D0-SUM(T(1:NG-1))
      ENDIF

C     Conserve momentum (mean) - first-order mom

      XMU_D(1:NATT) = XMUS(1:NATT)
      DO K=1,NG
        XMU_D(1:NATT) = XMU_D(1:NATT) - T(K)*XMU(K,1:NATT)
      END DO
      DO K=1,NG
        XMU(K,1:NATT) = XMU(K,1:NATT) + XMU_D(1:NATT)
      END DO
      
C     Conserve covar - second-order mom

      XVAR_D = XVARS
      DO K=1,NG
        DO J=1,NATT      
          DO I=1,J
            XVAR_D(I,J)=XVAR_D(I,J)-XVAR(K,I,J)*T(K)
            XVAR_D(I,J)=XVAR_D(I,J)-XMU(K,I)*XMU(K,J)*T(K)
          END DO
        END DO
      END DO      
      DO J=1,NATT
        DO I=1,J
          XVAR_D(I,J)=XVAR_D(I,J)+XMUS(I)*XMUS(J)
          XVAR_D(J,I)=XVAR_D(I,J)
        END DO
      END DO

      DO K=1,NG
       XVAR(K,1:NATT,1:NATT)=XVAR(K,1:NATT,1:NATT)+XVAR_D(1:NATT,1:NATT)
      END DO

      END SUBROUTINE ENFORCE_MOMS_CONSERVATION

c     std_EM_step
c     ###############################################################
      subroutine std_EM_step(NIND,NATT,NG,X,WX,W,GUM,DV,T,XMU,XVAR
     .                      ,V,WSUM,WTOT)

      implicit none

c     Call variables

      integer :: NIND,NATT,NG
      real(8) :: X(NIND,NATT),WX(NIND),GUM(NIND),W(NIND,NG),DV(NG)
     .          ,T(NG),XMU(NG,NATT),XVAR(NG,NATT,NATT),V(NG,NATT,NATT)
     .          ,WSUM(NG),WTOT
      
c     Local variables
      
      real(8) :: AL(NIND,NG)
      real(8) :: p_ik,T_OLD(NG),LWSUM,CN
      integer :: IER,NULL,JJ,KK,K

c     Begin program

      CN=1d0/(2d0*PI)**(NATT*0.5d0)
      
c     Calculate multivariate density of each point for every group

      DO JJ=1,NIND              !particle
        DO KK=1,NG              !group
          !Calculate component density
          CALL Mahalanobis(NATT,X(JJ,:),XMU(KK,:),V(KK,:,:),p_ik)
          !Compute particle probability per Gaussian
          p_ik=(-0.5d0)*MIN(p_ik,1d0*DENMAX)
          AL(JJ,KK)=EXP(p_ik)/SQRT(DV(KK))*CN
          !Calculate mixture density
          W(JJ,KK)=T(KK)*AL(JJ,KK) 
        END DO                  !group
        GUM(JJ) = sum(W(JJ,1:NG))
      END DO                    !particle
      
C     MSTEP Compute new estimate of mixing proportion (T) for each group.
C     Calculate sum of posterior probabilities in each group (=effective
C     number of points assigned to cluster K)
      DO K=1,NG
        LWSUM=0d0
        !Responsibility (of component K generating particle JJ)
        DO JJ=1,NIND            !particle
          W(JJ,K)= WX(JJ)*W(JJ,K)/GUM(JJ)
        END DO

        WSUM(K)=PSUM(W(:,K))
        T(K)=WSUM(K)/WTOT

        call update_gaussian_moments(NIND,NATT,K,X,W(:,K),DV(K)
     .                              ,XMU(K,:),XVAR(K,:,:)
     .                              ,V(K,:,:),WSUM(K))
          
      END DO                    !Gaussians

      !Recompute mixing coefficients to avoid round-off errors
      if(NG>1) then
        T(NG) = 1d0-sum(T(1:NG-1)) 
      else
        T(NG) = 1d0
      end if

C     Check conservation

      if (DIAG_IO) then
        CALL CHECK_CONSERVATION(NG,NATT,TS,XMUS,XVARS
     $            ,T(1:NG),XMU(1:NG,1:NATT),XVAR(1:NG,1:NATT,1:NATT))
        write (ML_IO,*)
      end if

      end subroutine std_EM_step

c     CW_AEM_step
c     ###############################################################
      subroutine CW_AEM_step(NIND,NATT,NG,X,WX,W,AL,GUM,DV,T,XMU,XVAR
     .                      ,V,WSUM,WTOT,ALIVE_GRP)
      ! This subroutine do adaptive EM step in a component-wise manner
      ! [Celeux2001] Truong 10/16/2020: We now know that this can be
      ! expensive and can be improved with simultaneous approach given
      ! the use of GS--Kmeans to avoid choosing too large value for the
      ! initial number of components.

      implicit none
      
c     Call variables

      INTEGER :: NIND,NATT,NG,ALIVE_GRP(:)
      REAL(8) :: X(NIND,NATT),WX(NIND),GUM(NIND),W(NIND,NG),DV(NG)
     .     ,T(NG),XMU(NG,NATT),XVAR(NG,NATT,NATT),V(NG,NATT,NATT)
     .     ,WSUM(NG),WTOT,AL(NIND,NG)

c     Local variables

      real(8) :: NPARBY2
      real(8) :: p_ik,T_OLD(NG),TK,ISQDV,CN,LWSUM,EIG_SPRD
      real(8) :: semiAL(NIND,NG),AL_JJ(NG),XMUK(NATT),VK(NATT,NATT)

      integer :: K,KK,JJ,NG_OLD,NG_Killed,K_tmp,NINDT

      logical :: KILL_G
      
c     Begin program

      NINDT = ipsum((/NIND/))

      CN=1d0/(2d0*PI)**(NATT*0.5d0)
      
      NPARBY2 = (NATT+NATT*(NATT+1)/2d0)*0.5d0 !T/2 (of restart paper)
      
      NG_OLD = NG
      NG_Killed = 0
      WSUM = 0d0

      ! Make a copy of Gaussians' proportions at the current
      ! iteration before updating it for merging purposes
      T_OLD(1:NG) = T(1:NG)
      
      !Cycle over gaussians
      K = 1
      DO WHILE(K.LE.NG_OLD)
        LWSUM = 0d0

        if (NIND > 0) then
          DO JJ=1,NIND          !particle
            W(JJ,K) = WX(JJ)*AL(JJ,K)/(GUM(JJ)+1e-15)
            LWSUM = LWSUM + W(JJ,K)
          END DO                !particle
        endif
#if defined(petsc)
        call MPI_Allreduce(lwsum,WSUM(K),1,MPI_DOUBLE_PRECISION
     .                    ,MPI_SUM,MPI_COMM_WORLD,mpierr)
#else
        WSUM(K) = LWSUM
#endif

C	Mixing proportion calc. by Figuireido's 2002 paper

        T(K) = max(WSUM(K)-NPARBY2,0d0)/(WTOT-NG_OLD*NPARBY2)

	!Weed out too-small Gaussians and/or update moments
	!The following IF-statement is consistent with algorithm
	!listed in Figueiredo's 2002 paper
c$$$        KILL_G = (T(K) .LT. 5D-3)

        !An easier condition for killing Gaussian:
        !Kill if it cannot support at least 2*NATT pcles, or if too elongated
        !(very inaccurate Cholesky factorization). Elongation is 
        !estimated by ratio of arithmetic over geometric eigenvalue average.

        EIG_SPRD = trace(NATT,XVAR(K,:,:))/NATT
     $           /(det(XVAR(K,:,:)))**(1d0/NATT)
        
        KILL_G = (T(K)*NINDT < 2**NATT).or.(EIG_SPRD > 2d1)

        if(KILL_G .and. NG > 1) then ! kill grp K 

          T(K) = 0D0
          XMU(K,1:NATT) = 0D0
          XVAR(K,1:NATT,1:NATT) = 0D0
          NG = NG - 1
          NG_Killed = NG_killed+1
          ! need to remove the probabilities of dead component from GUM
          DO JJ=1,NIND
            GUM(JJ)= GUM(JJ) - AL(JJ,K)!+1d-14 ! LC 08/11/23: broke Killian's test case
          END DO

        else

          !Eliminate killed Gaussians from list
          if(NG_Killed>0) then
            K_tmp = K-NG_killed
            T(K_tmp) = T(K) ; T(K) = 0d0          ! Truong's fix (01/31/2022)
            WSUM(K_tmp) = WSUM(K) ; WSUM(K) = 0d0
          else
            K_tmp = K
          end if
          
          ! Truong's fix (01/31/2022)
          T(1:NG_OLD) = T(1:NG_OLD)/SUM(T(1:NG_OLD)) ! Follow Figueiredo&Jain2002

          call update_gaussian_moments(NIND,NATT,K,X,W(:,K),DV(K_tmp)
     .                                ,XMU(K_tmp,:),XVAR(K_tmp,:,:)
     .                                ,V(K_tmp,:,:),WSUM(K_tmp))
           
C        Compute particle prob. semiAL(:,K) of grp K 
C	 after update its moments (means & cov)

          XMUK = XMU(K_tmp,:)
          VK   = V(K_tmp,:,:)
          ISQDV= CN/SQRT(DV(K_tmp))
          TK   = T(K_tmp)
          DO JJ=1,NIND              !particle
             !Calculate component density
             CALL Mahalanobis(NATT,X(JJ,:),XMUK,VK,p_ik)

             !Compute particle probability for Gaussian K
             p_ik=(-0.5d0)*MIN(p_ik,1d0*DENMAX)
             semiAL(JJ,K_tmp)=EXP(p_ik)*ISQDV
             !Component-wise update 
             GUM(JJ) = GUM(JJ) - AL(JJ,K) 
             AL(JJ,K_tmp)=TK*semiAL(JJ,K_tmp)
             GUM(JJ) = GUM(JJ) + AL(JJ,K_tmp)
          END DO                     !particle

        end if

        K = K + 1               !Next gaussian          

      END DO                    !Gaussians - end-do-while(K)

      !Recompute mixing coefficients to avoid round-off errors
      if(NG>1) then
         T(1:NG) = T(1:NG)/SUM(T(1:NG))         
         T(NG) = 1d0-sum(T(1:NG-1)) 
      else
         T(NG) = 1d0
      end if

      ALIVE_GRP(1:NG) = 1
      ALIVE_GRP(NG+1:) = 0
      
      end subroutine CW_AEM_step

c     SIM_AEM_step
c     ###############################################################
      subroutine SIM_AEM_step(NIND,NATT,NG,X,WX,W,AL,GUM,DV,T,XMU,XVAR
     .                      ,V,WSUM,WTOT,ALIVE_GRP)
      ! This subroutine do adaptive EM step in a simultaneous manner,
      ! reordering the Gaussian list so that finite Gaussians are contiguous.

      ! (TN -- 02/08/2021): while this works great for the simultaneous EM
      ! approach, it is not compatible with the current set up for AA,
      ! called in AMEM() subroutine. For now, we made a copy of this
      ! subroutine, named it SIM_AEM_step_2(), which does the same thing
      ! but without Gaussian list reordering.

      implicit none
      
c     Call variables

      INTEGER :: NIND,NATT,NG,ALIVE_GRP(:)
      REAL(8) :: X(NIND,NATT),WX(NIND),GUM(NIND),W(NIND,NG),DV(NG)
     .     ,T(NG),XMU(NG,NATT),XVAR(NG,NATT,NATT),V(NG,NATT,NATT)
     .     ,WSUM(NG),WTOT,AL(NIND,NG)

c     Local variables

      real(8) :: NPARBY2
      real(8) :: p_ik,T_OLD(NG),TK,ISQDV,CN,LWSUM,EIG_SPRD
      real(8) :: semiAL(NIND,NG),AL_JJ(NG),XMUK(NATT),VK(NATT,NATT)

      integer :: K,KK,JJ,NG_OLD,NG_Killed,K_tmp,NINDT

      logical :: KILL_G
      
c     Begin program

      NINDT = IPSUM((/NIND/))
c$$$#if defined(petsc)
c$$$      call MPI_Allreduce(NIND,NINDT,1,MPI_INTEGER
c$$$     .                  ,MPI_SUM,MPI_COMM_WORLD,mpierr)
c$$$#else
c$$$      NINDT = NIND
c$$$#endif

      CN=1d0/(2d0*PI)**(NATT*0.5d0)
      
      NPARBY2 = (NATT+NATT*(NATT+1)/2d0)*0.5d0 !T/2 (of restart paper)
      
      NG_OLD = NG
      NG_Killed = 0

      ! Make a copy of Gaussians' proportions at the current
      ! iteration before updating it for merging purposes
      T_OLD(1:NG) = T(1:NG)
      
      !Cycle over gaussians
      K = 1
      DO WHILE(K.LE.NG_OLD)
        LWSUM = 0d0

        if (NIND > 0) then
          DO JJ=1,NIND          !particle
            W(JJ,K) = WX(JJ)*AL(JJ,K)/GUM(JJ)
            LWSUM = LWSUM + W(JJ,K)
          END DO                !particle
        endif
#if defined(petsc)
        call MPI_Allreduce(lwsum,WSUM(K),1,MPI_DOUBLE_PRECISION
     .                    ,MPI_SUM,MPI_COMM_WORLD,mpierr)
#else
        WSUM(K) = LWSUM
#endif

C	Mixing proportion calc. by Figuireido's 2002 paper

        T(K) = max(WSUM(K)-NPARBY2,0d0)/(WTOT-NG_OLD*NPARBY2)

	!Weed out too-small Gaussians and/or update moments
	!The following IF-statement is consistent with algorithm
	!listed in Figueiredo's 2002 paper
c        if(T(K) .LT. 5D-6 .and. NG > 1) then	! kill grp K if T(K)==0

        !An easier condition for killing Gaussian:
        !Kill if it cannot support at least 2*NATT pcles, or if too elongated
        !(very inaccurate Cholesky factorization). Elongation is 
        !estimated by ratio of arithmetic over geometric eigenvalue average.

        EIG_SPRD = trace(NATT,XVAR(K,:,:))/NATT
     $           /(det(XVAR(K,:,:)))**(1d0/NATT)
        
        KILL_G = (T(K)*NINDT < 2**NATT).or.(EIG_SPRD > 2d1)

        if(KILL_G .and. NG > 1) then !Zero out Gaussian qtys

          T(K) = 0D0
          XMU(K,1:NATT) = 0D0
          XVAR(K,1:NATT,1:NATT) = 0D0
          WSUM(K) = 0d0 ; DV(K) = 0d0 ; V(K,:,:) = 0d0

          NG = NG - 1
          NG_Killed = NG_killed+1
        
        else                         !Recompute Gaussian qtys

          K_tmp = K-NG_killed
          T(K_tmp) = T(K)
          WSUM(K_tmp) = WSUM(K)

          call update_gaussian_moments(NIND,NATT,K,X,W(:,K),DV(K_tmp)
     .                                ,XMU(K_tmp,:),XVAR(K_tmp,:,:)
     .                                ,V(K_tmp,:,:),WSUM(K_tmp))

        end if

        K = K + 1               !Next gaussian          

      END DO                    !Gaussians - end-do-while(K)

      !Recompute mixing coefficients to avoid round-off errors
      if(NG>1) then
        T(1:NG) = T(1:NG)/SUM(T(1:NG))
        !T(NG) = 1d0-sum(T(1:NG-1)) 
      else
        T(NG) = 1d0
      end if

C     We do not need to update AL & GUM for the Picard solution
C     Since we don't know if Picard or AA solutions which one are
C     better (02/01/2021).
cC     Update GUM for computing XLOGL (to compare with the new one)
cC     not needed with component-wise update (already done during the iterations)
c      DO JJ=1,NIND              !particle
c        DO KK=1,NG             !group
c          AL(JJ,KK)=T(KK)*semiAL(JJ,KK)
c        END DO                  !group
c        !Calculate mixture density
c        GUM(JJ) = SUM(AL(JJ,1:NG))
c      END DO                     !particle

      ALIVE_GRP(1:NG) = 1
      ALIVE_GRP(NG+1:) = 0
      
      end subroutine SIM_AEM_step

c     SIM_AEM_step_2
c     ###############################################################
      subroutine SIM_AEM_step_2(NIND,NATT,NG_ORG,NG,X,WX,W,AL,GUM,DV
     .                      ,T,XMU,XVAR,V,WSUM,WTOT,ALIVE_GRP)
      ! This subroutine do adaptive EM step in a simultaneous manner
      ! without Gaussian list reordering. 
      ! This implementation is compatible with AA in AMEM() subroutine

      ! NG_ORG = starting/original number of components
      ! NG = current number of components @ current iteration

      implicit none

c     Call variables

      INTEGER :: NIND,NATT,NG_ORG,NG,ALIVE_GRP(:)
      REAL(8) :: X(NIND,NATT),WX(NIND),GUM(NIND),W(NIND,NG_ORG)
     .     ,DV(NG_ORG),T(NG_ORG),XMU(NG_ORG,NATT),XVAR(NG_ORG,NATT,NATT)
     .     ,V(NG_ORG,NATT,NATT),WSUM(NG_ORG),WTOT,AL(NIND,NG_ORG)

c     Local variables

      real(8) :: NPARBY2
      real(8) :: p_ik,T_OLD(NG_ORG),TK,ISQDV,CN,LWSUM,EIG_SPRD
      real(8) :: semiAL(NIND,NG_ORG),AL_JJ(NG_ORG)
     ,          ,XMUK(NATT),VK(NATT,NATT)

      integer :: K,KK,JJ,NG_OLD,K_tmp,NINDT

      logical :: KILL_G
      
c     Begin program

      NINDT = IPSUM((/NIND/))
c$$$#if defined(petsc)
c$$$      call MPI_Allreduce(NIND,NINDT,1,MPI_INTEGER
c$$$     .                  ,MPI_SUM,MPI_COMM_WORLD,mpierr)
c$$$#else
c$$$      NINDT = NIND
c$$$#endif

      CN=1d0/(2d0*PI)**(NATT*0.5d0)

      NPARBY2 = (NATT+NATT*(NATT+1)/2d0)*0.5d0 !T/2 (of restart paper)

      WSUM(1:NG_ORG) = 0d0
      NG_OLD = NG   ! NG=current number of alive grps

      ! Make a copy of Gaussians' proportions at the current
      ! iteration before updating it for merging purposes
      T_OLD(1:NG_ORG) = T(1:NG_ORG)

      !Cycle over gaussians
      DO K=1,NG_ORG
        IF (ALIVE_GRP(K) .EQ. 1) THEN
          LWSUM = 0d0

          if (NIND > 0) then
            DO JJ=1,NIND          !particle
              W(JJ,K) = WX(JJ)*AL(JJ,K)/GUM(JJ)
              LWSUM = LWSUM + W(JJ,K)
            END DO                !particle
          endif
#if defined(petsc)
          call MPI_Allreduce(lwsum,WSUM(K),1,MPI_DOUBLE_PRECISION
     .                    ,MPI_SUM,MPI_COMM_WORLD,mpierr)
#else
          WSUM(K) = LWSUM
#endif

          !Calculate grp's proportion using MML criterian
          T(K)=max(WSUM(K)-NPARBY2,0d0)/(WTOT-NG_OLD*NPARBY2)

          !Kill Gaussians 

          EIG_SPRD = trace(NATT,XVAR(K,:,:))/NATT
     $           /(det(XVAR(K,:,:)))**(1d0/NATT)
          KILL_G = (T(K)*NINDT < 2**NATT).or.(EIG_SPRD > 2d1)

          if (KILL_G .and. NG>1) then !Zero out corresponding qtys
            T(K) = 0D0  ; XMU(K,:) = 0D0   ; XVAR(K,:,:) = 0D0
            NG = NG-1   ; ALIVE_GRP(K) = 0
            DV(K) = 0d0 ; V(K,:,:) = 0d0   ; WSUM(K) = 0d0
          else                        !Recompute Gaussian qtys
            call update_gaussian_moments(NIND,NATT,K,X,W(:,K),DV(K)
     .                                ,XMU(K,:),XVAR(K,:,:)
     .                                ,V(K,:,:),WSUM(K))
          end if

        END IF   ! in ALIVE_GRP(K) == 1
      END DO     ! in K=1,NG_ORG

      ! Renormalization of Gaussian weights to have sum of 1
      T(1:NG_ORG) = T(1:NG_ORG)/SUM(T(1:NG_ORG))

      end subroutine SIM_AEM_step_2

c     update_gaussian_moments
c     ###############################################################
      subroutine update_gaussian_moments(NIND,NATT,K,X,W,DV,XMU,XVAR
     .                                  ,V,WSUM)
C     update XMU(K,:), XVAR(K,:,:) & compute inv & det of XVAR(K,:,:)

      implicit none

      integer :: NIND,NATT,K
      real(8) :: X(NIND,NATT),WX(NIND),GUM(NIND),W(NIND),DV
     .          ,XMU(NATT),XVAR(NATT,NATT),V(NATT,NATT),WSUM

      integer :: IER,NULL
      
      !Compute new estimates of group means and covariance matrices for group K
      CALL LCALK(NIND,NATT,X,W,WSUM,XMU,XVAR)

      !Obtain inverse and determinant of each estimated covariance matrix
      CALL GDETK(NATT,K,XVAR,V,DV,IER,NULL)
      
      end subroutine update_gaussian_moments

c$$$c     loglikelihood
c$$$c     ###############################################################
c$$$      subroutine loglikelihood(it,NATT,T,GUM,WX,XLOGL,MLENG)
c$$$
c$$$      implicit none
c$$$
c$$$c     Call variables
c$$$
c$$$      integer :: it,NATT
c$$$      real(8) :: T(:),GUM(:),WX(:),XLOGL,MLENG
c$$$      
c$$$c     Local variables
c$$$
c$$$      integer :: JJ,NG,NIND
c$$$      real(8) :: PENALIZED_TERM,LXLOG,tmpsum,NPARBY2
c$$$
c$$$c     Begin program
c$$$
c$$$      NPARBY2 = (NATT+NATT*(NATT+1)/2d0)*0.5d0 !T/2 (of restart paper)
c$$$      
c$$$      NG = size(T)
c$$$
c$$$      NIND = size(GUM)
c$$$      
c$$$      LXLOG = 0d0
c$$$      DO JJ = 1,NIND
c$$$        LXLOG = LXLOG+LOG(GUM(JJ))*WX(JJ)
c$$$      END DO
c$$$
c$$$#if defined(petsc)
c$$$      call MPI_Allreduce(LXLOG,XLOGL,1
c$$$     .                  ,MPI_DOUBLE_PRECISION
c$$$     .                  ,MPI_SUM,MPI_COMM_WORLD,mpierr)
c$$$#else
c$$$      XLOGL = LXLOG
c$$$#endif
c$$$      
c$$$      tmpsum = sum(LOG(T(1:NG)))
c$$$
c$$$C     --------------------------------------------------------
c$$$C     Formula for MLENG in Figueiredo's matlab code
c$$$C     This formula is not consistent with the paper
c$$$c      MLENG = -XLOGL + NPARBY2*tmpsum + 
c$$$c     .       (NPARBY2+0.5D0)*NG*LOG(NINDT*1d0)
c$$$
c$$$C     This formula is consistent with Figueiredo's paper 
c$$$      MLENG = -XLOGL+ NPARBY2*tmpsum + 0.5d0*NG
c$$$     .       *LOG(NINDT/12d0)+NG*(NPARBY2+1d0)*0.5d0
c$$$C     --------------------------------------------------------
c$$$
c$$$C    This formula is consistent with our GMM restart paper
c$$$      PENALIZED_TERM = -NPARBY2*tmpsum
c$$$     .                 -(NG*NPARBY2+0.5*NG-0.5)*LOG(1d0*NINDT) 
c$$$
c$$$      XLOGL = XLOGL + PENALIZED_TERM
c$$$
c$$$      if (DIAG_IO) write(ML_IO,*) it,XLOGL,MLENG,MINVAL(T(1:NG))
c$$$
c$$$      end subroutine loglikelihood

c     COMPUTE_PLLH
c     ###############################################################
      SUBROUTINE COMPUTE_PLLH(NATT,T,GUM,WX,ALIVE_GRP,XLOGL,MLENG)
c     Given GUM(1:NIND) = probability of X(J,:) in the mixture, we compute
c     the XLOGL, the (PENALIZED) LOGLIKELIHOOD value, of GMM,
c     PLLH(\theta) = XLOGL

      IMPLICIT NONE

      !	Call variables
      INTEGER, INTENT(IN) :: NATT
      INTEGER, INTENT(IN) :: ALIVE_GRP(:)
      REAL(8), INTENT(IN) :: WX(:),GUM(:),T(:)
      REAL(8), INTENT(INOUT) :: XLOGL,MLENG

      !	Local variables
      INTEGER :: JJ,K,NG_ORG,NG,NIND,NINDT
      REAL(8) :: NPARBY2,DVAL,X_WGHTS,LXLOGL,LOGT

      ! Begin program

      NG_ORG = size(ALIVE_GRP)
      NG = count(ALIVE_GRP == 1)
      NIND = size(GUM)
      XLOGL = 0d0

      NINDT = IPSUM((/NIND/))
      
      NPARBY2 = 0.5D0*(NATT+NATT*(NATT+1)*0.5)

      ! Compute standard LLH (taking into account particles' weights)
      XLOGL = PSUM(WX*LOG(GUM))

      LOGT = sum(LOG(T),mask=(ALIVE_GRP==1))

      ! Take into account the penalized term in adaptive-EM approach:
      DVAL = (NG*NPARBY2+0.5*NG-0.5)  !Figuereido & Jain minus one to account for weight norm.
c$$$      DVAL = (NG*NPARBY2+0.5*NG)  !Figuereido & Jain
      XLOGL = XLOGL -NPARBY2*LOGT - DVAL*LOG(1.0D0*NINDT)
    
      MLENG = -XLOGL + NPARBY2*LOGT + 0.5d0*NG
     .       *LOG(NINDT/12d0)+NG*(NPARBY2+1d0)*0.5d0

      END SUBROUTINE COMPUTE_PLLH

c     split_gaussian
c     ###############################################################
      function split_gaussian(isplit,TOL,NIND,NATT,NG,X,WX,W,GUM,DV,T
     .                       ,XMU,XVAR,V,WSUM,XLOGL,MLENG)
     .         result(split)

      implicit none

c     Call variables
      
      logical :: split

      integer :: isplit

      INTEGER :: NIND,NATT,NG
      REAL(8) :: X(NIND,NATT),WX(NIND),GUM(NIND),W(NIND,NG),DV(NG)
     .     ,T(NG),XMU(NG,NATT),XVAR(NG,NATT,NATT),V(NG,NATT,NATT)
     .     ,WSUM(NG),WTOT,XLOGL,MLENG,TOL

c     Local variables
      
      real(8) :: SDIFF,SALIM,Np,Nd,KOK,K0
      real(8),save :: SXLOGL(MITERMAX),SMLENG(MITERMAX)
      integer :: kspl

      !check for split-loop convergence
      IF(ISPLIT>1)then
        SDIFF=SMLENG(ISPLIT)-SMLENG(ISPLIT-1)
cc        SALIM=TOL*XLOGL(ISPLIT)
        SALIM=TOL*0.1d0
      ELSE
        SDIFF = 1d0
        SALIM = 0d0        
      END IF

      IF(ISPLIT>MAXSPLIT .OR. ABS(SDIFF)<ABS(SALIM)) THEN
        SPLIT = .false.
      ELSE
        SPLIT = .true.

        SXLOGL(ISPLIT) = XLOGL
        SMLENG(ISPLIT) = MLENG
       
        !find split candidate
        call find_split_cand(NG,NATT,XMU(1:NG,1:NATT)
     .                      ,XVAR(1:NG,1:NATT,1:NATT),T(1:NG),Kspl)
        
        !split the candidate
        KOK = KURTOSIS(NIND,NATT,X(1:NIND,1:NATT),W(1:NIND,Kspl)
     $       ,XMU(Kspl,1:NATT),WSUM(Kspl),V(Kspl,1:NATT,1:NATT))
        Np = WSUM(Kspl)
        Nd = DBLE(NATT)
        K0 = EKurt(Np,Nd)

        if (STD_IO) print *
c$$$        print *,"To split component #",kspl
c$$$        print *,"  Kurtosis=",KOK," K0=",K0
        Call splitK(Np,Kspl,NG,NATT,XMU(1:NG+1,1:NATT),XVAR(1:NG+1
     $       ,1:NATT,1:NATT),T(1:NG+1),DV(1:NG+1),V(1:NG+1,1:NATT
     $       ,1:NATT))

        ISPLIT=ISPLIT+1
      END IF

      if (STD_IO) print *,"End SPLIT-loop#",ISPLIT-1
      if (STD_IO) print *

      end function split_gaussian

      end module EM_utils
