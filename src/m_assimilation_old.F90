module m_assimilation_old

contains

subroutine assimilation_old(obs,win,winref,DA,tini,tfin,nrobs,iter,obsoloc,obsaloc,obsotimes,obsatimes)
   use mod_dimensions
   use mod_state
   use mod_observation
   use m_readinfile
   use m_enkfprep
   use m_getD
   use m_wtime
   implicit none
   integer, intent(in) :: nrobs
   integer, intent(in) :: iter
   type(state),    intent(inout) :: win(0:nrt,nrens)
   type(state),    intent(in) :: winref(0:nrt)
   type(observation),  intent(in):: obs(nrobs)
   integer, intent(in) :: tini
   integer, intent(in) :: tfin
   integer, intent(in) :: obsoloc(nro)
   integer, intent(in) :: obsaloc(nra)
   integer, intent(in) :: obsotimes(nrt)
   integer, intent(in) :: obsatimes(nrt)
   real,  intent(in) :: DA(nrobs,nrens)

   real start,finish,cpu0,cpu1
   logical, allocatable :: lobs(:,:)
   integer, allocatable :: nobs(:)

   real, allocatable :: dist(:,:) ! Distance between obs and grid point
   real, allocatable :: Y(:,:)
   real, allocatable :: S(:,:)
   real, allocatable :: E(:,:)
   real, allocatable :: D(:,:)
   real, allocatable :: R(:,:)
   real, allocatable :: meanS(:)
   real, allocatable :: innov(:)

   real, allocatable, dimension(:,:)   :: subS,subE,subD, subR, subdist
   real, allocatable, dimension(:)     :: subinnov,corrA,corrO
   real, allocatable, dimension(:,:,:)   :: subwin
   integer i,j,l,m
   real stdA,aveA,stdO,aveO,stdS
   real Einfl,truncdist
   real lcovb


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

   cpu0=wtime()
   call cpu_time(start)
   if (local==0) then
      print '(tr5,a,2(i3,a))','        assimilation_old: iter=',iter,' -> Calling analysis with mode: ',mode_analysis
      allocate(subwin(ndim,tini:tfin,nrens))
      do j=1,nrens
      do l=tini,tfin
         subwin(1:nx,l,j)=win(l,j)%Atmos(1:nx)
         subwin(nx+1:2*nx,l,j)=win(l,j)%Ocean(1:nx)
      enddo
      enddo

      call analysis(subwin, R, E, S, D, innov, ndim*(nrw+1), nrens, nrobs, .true., truncation, mode_analysis, &
                              lrandrot, lupdate_randrot, lsymsqrt, inflate, infmult, 1)

      do j=1,nrens
      do l=tini,tfin
         win(l,j)%Atmos(1:nx)=subwin(1:nx,l,j)
         win(l,j)%Ocean(1:nx)=subwin(nx+1:2*nx,l,j)
      enddo
      enddo
      deallocate(subwin)

   else
      print '(a,2i2,f10.2,f10.3)','        assimilation_old: calling local analysis with mode: '&
                                           ,mode_analysis,local,obs_radius,obs_truncation

      allocate(lobs(nx,nrobs),nobs(nx),dist(nx,nrobs))
      lobs(:,:)=.false.
      nobs(:)=0

! Distace based localization compute active measurements for each grid point
      if (local == 1) then
         truncdist=obs_radius
         do m=1,nrobs
            do i=1,nx
               dist(i,m)=real(abs(obs(m)%xloc-i))
               if (dist(i,m) < truncdist) then
                  lobs(i,m)=.true.
                  nobs(i)=nobs(i)+1
               endif
            enddo
            if (i==nx/2)  print *,'        number of local measurements for i=',i,nobs(i)
         enddo
      endif

! Adaptive localization compute active measurements for each grid point
      if (local == 2) then
         truncdist=1.0-obs_truncation
         allocate(corrA(nrobs), corrO(nrobs))
!!$OMP PARALLEL DO private(aveA,aveO,stdA,stdO,j,m,stdS,corrA,corrO)
         do i=1,nx
            stdA=0.0
            aveA=0.0
            stdO=0.0
            aveO=0.0
            do j=1,nrens
               aveA=aveA+win(tfin,j)%Atmos(i)
               aveO=aveO+win(tfin,j)%Ocean(i)
               stdA=stdA+win(tfin,j)%Atmos(i)**2
               stdO=stdO+win(tfin,j)%Ocean(i)**2
            enddo
            aveA=aveA/real(nrens)
            aveO=aveO/real(nrens)
            stdA=stdA/real(nrens)
            stdO=stdO/real(nrens)

            stdA=stdA-aveA**2
            stdO=stdO-aveO**2

            stdA=sqrt(stdA)
            stdO=sqrt(stdO)

            do m=1,nrobs
               corrA(m)=0.0
               corrO(m)=0.0
               stdS=0.0
               do j=1,nrens
                  stdS=stdS+S(m,j)**2
                  corrA(m)=corrA(m)+S(m,j)*win(tfin,j)%Atmos(i)
                  corrO(m)=corrO(m)+S(m,j)*win(tfin,j)%Ocean(i)
               enddo
               stdS=stdS/real(nrens)
               stdS=sqrt(stdS)

               corrA(m)=corrA(m)/real(nrens)
               corrO(m)=corrO(m)/real(nrens)
               corrA(m)=corrA(m)/(stdS*stdA+tiny(stdA))
               corrO(m)=corrO(m)/(stdS*stdO+tiny(stdA))
            enddo

            do m=1,nrobs
               if ((abs(corrA(m)) > obs_truncation) .or. (abs(corrO(m)) > obs_truncation)) then
                  nobs(i)=nobs(i)+1
                  lobs(i,m)=.true.
                  dist(i,m)=1.0-max(abs(corrA(m)),abs(corrO(m)))
               endif
            enddo
            if (i==nx/2)  print *,'       assimilation_old: number of local measurements for i=',i,nobs(i)
         enddo
!!$OMP END PARALLEL DO
      endif


      do i=1,nx
         if (nobs(i) > 0) then
            allocate(subD(nobs(i),nrens))
            allocate(subdist(nobs(i),1))
            allocate(subE(nobs(i),nrens))
            allocate(subS(nobs(i),nrens))
            allocate(subinnov(nobs(i)))
            allocate(subR(nobs(i),nobs(i)))
            call getD(D,subD,nrobs,nrens,lobs(i,:),nobs(i)) ! the innovations
            call getD(E,subE,nrobs,nrens,lobs(i,:),nobs(i)) ! the observation errors
            call getD(S,subS,nrobs,nrens,lobs(i,:),nobs(i)) ! the HA'

            if (lcovsmooth) then
               lcovb=truncdist/(2.0*log(obsdamping))
               call getD(dist(i,:),subdist,nrobs,1,lobs(i,:),nobs(i)) ! distances
               do m=1,nobs(i)
                  if (subdist(m,1) > real(truncdist)/2.0 ) then
                     Einfl=exp((subdist(m,1)-real(truncdist)/2.0)/lcovb)
                     subE(m,:)=subE(m,:)*Einfl
                  endif
               enddo
            endif

            allocate(subwin(2,tini:tfin,nrens))
            do j=1,nrens
            do l=tini,tfin
               subwin(1,l,j)=win(l,j)%Atmos(i)
               subwin(2,l,j)=win(l,j)%Ocean(i)
            enddo
            enddo

             call analysis(subwin, subR, subE, subS, subD, subinnov, 2*(tfin-tini+1), nrens, nobs(i),&
                           .false., truncation, mode_analysis, &
                           lrandrot, lupdate_randrot, lsymsqrt, inflate, infmult, 1)

            do j=1,nrens
            do l=tini,tfin
               win(l,j)%Atmos(i)=subwin(1,l,j)
               win(l,j)%Ocean(i)=subwin(2,l,j)
            enddo
            enddo

            deallocate(subD, subE, subS, subinnov, subR, subwin, subdist)
         endif
      enddo
      print *,'       assimilation_old: local analysis done'

      if (local== 2) deallocate(corrA, corrO)
      deallocate(lobs,nobs)
      if (allocated(dist)) deallocate(dist)
   endif

   deallocate(Y,S,E,D,R,meanS,innov)
   cpu1=wtime()
   call cpu_time(finish)
   print '(2(a,f6.2))','        assimilation_old: cpu time=',finish-start,', wall-clock time=',cpu1-cpu0

end subroutine
end module
