module m_assimilation_old

contains

subroutine assimilation_old(obs,win,winref,DA,tini,tfin,nrobs,iter,obsoloc,obsaloc,obsotimes,obsatimes)
   use mod_dimensions
   use mod_state
   use mod_observation
   use m_readinfile
   use m_enkfprep
   implicit none
   integer, intent(in) :: nrobs
   integer, intent(in) :: iter
   type(state),    intent(in) :: win(0:nrt,nrens)
   type(state),    intent(in) :: winref(0:nrt)
   type(observation),  intent(in):: obs(nrobs)
   integer, intent(in) :: tini
   integer, intent(in) :: tfin
   integer, intent(in) :: obsoloc(nro)
   integer, intent(in) :: obsaloc(nra)
   integer, intent(in) :: obsotimes(nrt)
   integer, intent(in) :: obsatimes(nrt)
   real,  intent(in) :: DA(nrobs,nrens)

   real, allocatable :: Y(:,:)
   real, allocatable :: S(:,:)
   real, allocatable :: E(:,:)
   real, allocatable :: D(:,:)
   real, allocatable :: R(:,:)
   real, allocatable :: meanS(:)
   real, allocatable :: innov(:)

   allocate(Y(nrobs,nrens))
   allocate(S(nrobs,nrens))
   allocate(E(nrobs,nrens))
   allocate(D(nrobs,nrens))
   allocate(R(nrobs,nrobs))
   allocate(meanS(nrobs))
   allocate(innov(nrobs))

   print '(tr5,a,i3,a)','main: iter=',iter,' -> Calling enkf preprep'
   ! Returns all matrices UNSCALED by sqrt(nrens-1) for use in analysis.F90
   call enkfprep(obs,Y,S,E,D,DA,meanS,R,innov,winref,win,nrobs,&
                 tini,tfin,obsoloc,obsaloc,obsotimes,obsatimes)
   D=D-Y

!   if (.not.local) then
   print '(tr5,a,2(i3,a))','main: iter=',iter,' -> Calling analysis with mode: ',mode_analysis
   call analysis(win(tini:tfin,:), R, E, S, D, innov, ndim*(nrw+1), nrens, nrobs, .true., truncation, mode_analysis, &
                           lrandrot, lupdate_randrot, lsymsqrt, inflate, infmult, 1)

!   else 
!      print '(a,i2)','       enkf: calling local analysis with mode: ',mode_analysis,local
!      icall=0
!      do i=1,nx
!         nobs=0
!         lobs(:)=.false.
!
!
!
!! Distace based localization
!         if (local == 1) then
!            do m=1,nrobs
!               if (real(abs(obspos(m)-i)) < robs) then
!                  lobs(m)=.true.
!                  nobs=nobs+1
!               endif
!            enddo
!         endif
!
!! Adaptive localization
!         if(local == 2) then
!            stdA=0.0
!            aveA=0.0
!            do l=1,nrens
!               aveA=aveA+mem(i,l)
!               stdA=stdA+mem(i,l)*mem(i,l)
!            enddo
!            aveA=aveA/real(nrens)
!            stdA=stdA/real(nrens)
!            stdA=stdA-aveA**2
!            stdA=sqrt(stdA)
!
!            do m=1,nrobs
!               corr(m)=0.0
!               stdS=0.0
!               do l=1,nrens
!                  stdS=stdS+S(m,l)*S(m,l)
!                  corr(m)=corr(m)+S(m,l)*mem(i,l)
!               enddo
!               stdS=sqrt(stdS/real(nrens))
!               corr(m)=corr(m)/real(nrens)
!               corr(m)=corr(m)/(stdS*stdA)
!            enddo
!
!            do m=1,nrobs
!               if (abs(corr(m)) > obstreshold) then
!                  lobs(m)=.true.
!                  nobs=nobs+1
!               endif
!            enddo
!         endif
!
!
!         if (nobs > 0) then
!            allocate(subD(nobs,nrens))
!            allocate(subE(nobs,nrens))
!            allocate(subS(nobs,nrens))
!            allocate(subR(nobs,nobs))
!            call getD(D,subD,nrobs,nrens,lobs,nobs) ! the innovations to use 
!            call getD(E,subE,nrobs,nrens,lobs,nobs) ! the observation errors to use
!            call getD(S,subS,nrobs,nrens,lobs,nobs) ! the HA' to use
!            subR=matmul(subE,transpose(subE))/float(nrens)
!            allocate(subinnovation(nobs))
!            allocate(submem(1,nrens))
!            l=0
!            do m =1, nrobs
!               if (lobs(m)) then
!                  l=l+1
!                  subinnovation(l)      = obs(m)- meanS(m)
!               endif
!            enddo
!
!
!            submem(1,:)=mem(i,:)
!            icall=icall+1
!            if (lupdate_randrot .and. icall==1) then
!                local_rot=.true.
!            else
!                local_rot=.false.
!            endif
!            call analysis(submem, subR, subE, subS, subD, subinnovation, 1, nrens, nobs, .false., truncation, mode_analysis, &
!                         lrandrot, local_rot, lsymsqrt, inflate, infmult, 1)
!            mem(i,:)=submem(1,:)
!            deallocate(subD, subE, subS, subR, subinnovation, submem)
!
!   endif
   deallocate(Y,S,E,D,R,meanS,innov)

end subroutine
end module
