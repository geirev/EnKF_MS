module m_enkf
contains
subroutine enkf(mem,nrens,obs,nroceanobs,nratmosobs,mode_analysis,truncation,covmodel,dx,rh,Rexact,rd,&
               &lrandrot,lsymsqrt,&
               &inflate, infmult,&
               &local,robs, obstreshold)

   use mod_dimensions
   use mod_state
   use mod_observation
   use m_obs_pert
   use m_getD
   implicit none
   integer, intent(in) :: nrens
   integer, intent(in) :: nroceanobs
   integer, intent(in) :: nratmosobs
   integer, intent(in) :: mode_analysis

   type(state),          intent(inout) :: mem(nrens)
   type(observation),    intent(in)    :: obs(nroceanobs+nratmosobs)
   logical, intent(in) :: Rexact

   logical, intent(in) :: lrandrot
   logical, intent(in) :: lsymsqrt

   integer, intent(in) :: local
   real,    intent(in) :: robs          ! influence radii for the measurements
   real,    intent(in) :: obstreshold

   integer, intent(in) :: inflate
   real, intent(in)    :: infmult

   real,    intent(in) :: truncation
   character(len=8), intent(in) :: covmodel
   real,    intent(in) :: dx
   type(substate),   intent(in) :: rh
   real,    intent(in) :: rd


   real, allocatable :: R(:,:)
   real, allocatable :: E(:,:)
   real, allocatable :: D(:,:)
   real, allocatable :: S(:,:)
   real, allocatable :: meanS(:)
   real, allocatable :: innovation(:)
   real, allocatable :: scaling(:)

   integer m,i,j
   logical :: lupdate_randrot=.true.

   integer  nrobs
   nrobs=nroceanobs+nratmosobs
   if (mode_analysis == 0) return

! Local analysis variables
!   integer l,icall
!   logical lobs(nrobs)           ! which measurements are active
!   integer nobs
!   logical local_rot
!   real corr(nrobs),stdA,aveA,stdS
!   real, allocatable, dimension(:,:) :: subS,subE,subD, subR
!   real, allocatable, dimension(:,:)   :: submem
!   real, allocatable :: subinnovation(:)


! End Local analysis variables

   allocate(E(nrobs,nrens))
   allocate(D(nrobs,nrens))
   allocate(S(nrobs,nrens))
   allocate(meanS(nrobs))
   allocate(scaling(nrobs))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Observe ensemble to construct the matrix S=HA
   do j =1, nrens
      do m =1, nroceanobs
         S(m,j)      =  mem(j)%ocean(obs(m)%pos)
      enddo
      i=nroceanobs
      do m =1,nratmosobs
         S(i+m,j)      =  mem(j)%atmos(obs(i+m)%pos)
      enddo
   enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Construct observation perturbations E
   call obs_pert(E(1:nroceanobs,1:nrens),nrens,nroceanobs,.true.,dx,rh%ocean,covmodel,obs(1:nroceanobs)%pos)
   i=nroceanobs
   call obs_pert(E(i+1:i+nratmosobs,1:nrens),nrens,nratmosobs,.true.,dx,rh%atmos,covmodel,obs(i+1:i+nratmosobs)%pos)

! Introduce correct variances
   do j=1,nrens
      do m=1,nrobs
         E(m,j)=sqrt(obs(m)%var)*E(m,j)
      enddo
   enddo


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Construct ensemble of measurements D=d+E
   do j=1,nrens
      do m=1,nrobs
         D(m,j)=obs(m)%d+E(m,j)
      enddo
   enddo


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Compute innovation D'=D-HA
   D=D-S
   write(*,*)'D:'
   write(*,'(10F10.4)')S(1:10,1:10)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Compute mean(HA)
   meanS=0.0
   do j=1,nrens
   do m=1,nrobs
      meanS(m)=meanS(m)+S(m,j)
   enddo
   enddo
   meanS=(1.0/float(nrens))*meanS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Compute HA'=HA-mean(HA)
   do j=1,nrens
      S(:,j)=S(:,j)-meanS(:)
   enddo
   write(*,*)'S:'
   write(*,'(10F10.4)')S(1:10,1:10)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   allocate(R(nrobs,nrobs))
   R=0.0

   if (Rexact) then
      print '(a,a)','       enkf: Exact R using covariance model: ',trim(covmodel)
      select case (trim(covmodel))
      case ('diagonal')
         do m=1,nrobs
            R(m,m)=obs(m)%var
         enddo
      case ('gaussian')
         do i=1,nrobs
         do j=1,nrobs
            R(i,j)=sqrt(obs(i)%var)*sqrt(obs(j)%var)*exp(-real(obs(i)%pos-obs(j)%pos)**2/rd**2)
         enddo
         enddo
      case default
         print '(a,a)','       enkf: covmodel is invalid : ',trim(covmodel)
      end select
   else
      print '(a,a)','       enkf: lowrank R using covariance model: ',trim(covmodel)
      R=matmul(E,transpose(E))/float(nrens-1)
   endif
   write(*,*)'R:'
   write(*,'(10F10.4)')R(1:10,1:10)


   allocate(innovation(nrobs))
   do m =1, nrobs
      innovation(m)      = obs(m)%d- meanS(m)
   enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Scaling of matrices
   print '(a)','       enkf: scale matrices'
   do m=1,nrobs
      scaling(m)=1./sqrt(R(m,m))
      S(m,:)=scaling(m)*S(m,:)
      E(m,:)=scaling(m)*E(m,:)
      D(m,:)=scaling(m)*D(m,:)
      innovation(m)=scaling(m)*innovation(m)
   enddo

   do j=1,nrobs
   do i=1,nrobs
      R(i,j)=scaling(i)*R(i,j)*scaling(j)
   enddo
   enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Computing analysis
   if (local==0) then
      print '(a,i2)','       enkf: calling global analysis with mode: ',mode_analysis
      call analysis(mem, R, E, S, D, innovation, ndim, nrens, nrobs, .true., truncation, mode_analysis, &
                      lrandrot, lupdate_randrot, lsymsqrt, inflate, infmult, 1)

!LOC   else
!LOC      print '(a,i2)','       enkf: calling local analysis with mode: ',mode_analysis,local
!LOC      icall=0
!LOC      do i=1,nx
!LOC         nobs=0
!LOC         lobs(:)=.false.
!LOC
!LOC
!LOC
!LOC! Distace based localization
!LOC         if (local == 1) then
!LOC            do m=1,nrobs
!LOC               if (real(abs(obs(m)%pos-i)) < robs) then
!LOC                  lobs(m)=.true.
!LOC                  nobs=nobs+1
!LOC               endif
!LOC            enddo
!LOC         endif
!LOC
!LOC! Adaptive localization
!LOC         if(local == 2) then
!LOC            stdA=0.0
!LOC            aveA=0.0
!LOC            do l=1,nrens
!LOC               aveA=aveA+mem(i,l)
!LOC               stdA=stdA+mem(i,l)*mem(i,l)
!LOC            enddo
!LOC            aveA=aveA/real(nrens)
!LOC            stdA=stdA/real(nrens)
!LOC            stdA=stdA-aveA**2
!LOC            stdA=sqrt(stdA)
!LOC
!LOC            do m=1,nrobs
!LOC               corr(m)=0.0
!LOC               stdS=0.0
!LOC               do l=1,nrens
!LOC                  stdS=stdS+S(m,l)*S(m,l)
!LOC                  corr(m)=corr(m)+S(m,l)*mem(i,l)
!LOC               enddo
!LOC               stdS=sqrt(stdS/real(nrens))
!LOC               corr(m)=corr(m)/real(nrens)
!LOC               corr(m)=corr(m)/(stdS*stdA)
!LOC            enddo
!LOC
!LOC            do m=1,nrobs
!LOC               if (abs(corr(m)) > obstreshold) then
!LOC                  lobs(m)=.true.
!LOC                  nobs=nobs+1
!LOC               endif
!LOC            enddo
!LOC         endif
!LOC
!LOC
!LOC         if (nobs > 0) then
!LOC            allocate(subD(nobs,nrens))
!LOC            allocate(subE(nobs,nrens))
!LOC            allocate(subS(nobs,nrens))
!LOC            allocate(subR(nobs,nobs))
!LOC            call getD(D,subD,nrobs,nrens,lobs,nobs) ! the innovations to use
!LOC            call getD(E,subE,nrobs,nrens,lobs,nobs) ! the observation errors to use
!LOC            call getD(S,subS,nrobs,nrens,lobs,nobs) ! the HA' to use
!LOC            subR=matmul(subE,transpose(subE))/float(nrens)
!LOC            allocate(subinnovation(nobs))
!LOC            allocate(submem(1,nrens))
!LOC            l=0
!LOC            do m =1, nrobs
!LOC               if (lobs(m)) then
!LOC                  l=l+1
!LOC                  subinnovation(l)      = obs(m)- meanS(m)
!LOC               endif
!LOC            enddo
!LOC
!LOC
!LOC            submem(1)=mem(i)
!LOC            icall=icall+1
!LOC            if (lupdate_randrot .and. icall==1) then
!LOC                local_rot=.true.
!LOC            else
!LOC                local_rot=.false.
!LOC            endif
!LOC            call analysis(submem, subR, subE, subS, subD, subinnovation, 1, nrens, nobs, .false., truncation, mode_analysis, &
!LOC                         lrandrot, local_rot, lsymsqrt, inflate, infmult)
!LOC            mem(i)=submem(1)
!LOC            deallocate(subD, subE, subS, subR, subinnovation, submem)
!LOC         endif
!LOC      enddo
   endif
   print '(a)','       enkf: done'

   deallocate(innovation)
   deallocate(R)


end subroutine
end module
