module m_enkfprep
contains
subroutine enkfprep(mem,obs,Y,S,E,D,DA,meanY,R,innov,winref,win,nrobs,lwin,tini,tfin,obsoloc,obsaloc,obsotimes,obsatimes)
   use mod_dimensions
   use mod_state
   use mod_observation
   use m_pseudo1D
   use m_random
   use m_readinfile
   use m_model
   implicit none
   type(state),    intent(inout) :: mem(nrens)
   type(state),    intent(in) :: win(0:nrw,nrens)
   type(state),    intent(in) :: winref(0:nrw)
   integer, intent(in) :: nrobs
   integer, intent(in) :: lwin
   integer, intent(in) :: tini
   integer, intent(in) :: tfin
   integer, intent(in) :: obsoloc(nro)
   integer, intent(in) :: obsaloc(nra)
   integer, intent(in) :: obsotimes(nrt)
   integer, intent(in) :: obsatimes(nrt)
   type(observation),  intent(out):: obs(nrobs)
   real,  intent(out):: Y(nrobs,nrens)
   real,  intent(out):: S(nrobs,nrens)
   real,  intent(out):: E(nrobs,nrens)
   real,  intent(out):: D(nrobs,nrens)
   real,  intent(in) :: DA(nrobs,nrens)
   real,  intent(out):: R(nrobs,nrobs)
   real,  intent(out):: meanY(nrobs)
   real,  intent(out):: innov(nrobs)

   integer m,i,k,mm,m1,m2,j

   real tmpo(nro)
   real tmpa(nra)
   real, allocatable :: EEfield(:,:)
   real, allocatable :: scaling(:)

   m=0
   do i=1,nrt
      if (tini < obsotimes(i) .and. obsotimes(i) <= tfin .and. nro > 0) then
         do k=1,nrw
            if (obsotimes(i) == (lwin-1)*nrw+k ) then
               call random(tmpo,nro)
               m1=m+1
               m2=m+nro
               do mm=1,nro
                  m=m+1
!                  obs(m)%xloc=obsoloc(mm)
!                  obs(m)%tloc=obsotimes(i)
!                  obs(m)%var=obsvar%ocean
!                  obs(m)%observed='ocean'
!                  obs(m)%d=winref(k)%ocean(obs(m)%xloc) + sqrt(obsvar%ocean)*tmpo(mm)
                  Y(m,:) = win(k,:)%ocean(obs(m)%xloc)
               enddo
            endif
         enddo
      endif

      if (tini < obsatimes(i) .and. obsatimes(i) <= tfin .and. nra > 0) then
         do k=1,nrw
            if (obsatimes(i) == (lwin-1)*nrw+k ) then
               call random(tmpa,nra)
               m1=m+1
               m2=m+nra
               do mm=1,nra
                  m=m+1
!                  obs(m)%xloc=obsaloc(mm)
!                  obs(m)%tloc=obsatimes(i)
!                  obs(m)%var=obsvar%atmos
!                  obs(m)%observed='atmos'
!                  obs(m)%d=winref(k)%atmos(obs(m)%xloc) + sqrt(obsvar%atmos)*tmpa(mm)
                  Y(m,:) = win(k,:)%atmos(obs(m)%xloc)
               enddo
            endif
         enddo
      endif
   enddo

!   do m=1,nrobs
!      print '(a,i5,f13.5)','obs(m)%d: ',m,obs(m)%d
!   enddo
!   stop

   m=0
   do i=1,nrt
      if (tini < obsotimes(i) .and. obsotimes(i) <= tfin .and. nro > 0) then
         do k=1,nrw
            if (obsotimes(i) == (lwin-1)*nrw+k ) then
               m1=m+1
               m2=m+nro
               do mm=1,nro
                  m=m+1
               enddo
               if (trim(covmodel) == 'diagonal') then
                  call random(E(m1:m2,1:nrens),nro*nrens)
                  E(m1:m2,:)=sqrt(obsvar%ocean)*E(m1:m2,:)
               elseif(trim(covmodel) == 'gaussian') then
                  allocate(EEfield(nx,nrens))
                  call pseudo1D(EEfield,nx,nrens,rd,dx,nx)
                  do j=1,nrens
                  do m=m1,m2
                     E(m,j)=EEfield(obsaloc(m),j)
                  enddo
                  enddo
                  deallocate(EEfield)
                  E(m1:m2,:)=sqrt(obsvar%ocean)*E(m1:m2,:)
               else
                  stop 'invalid covariance model'
               endif
            endif
         enddo
      endif

      if (tini < obsatimes(i) .and. obsatimes(i) <= tfin .and. nra > 0) then
         do k=1,nrw
            if (obsatimes(i) == (lwin-1)*nrw+k ) then
               m1=m+1
               m2=m+nra
               do mm=1,nra
                  m=m+1
               enddo
               if (trim(covmodel) == 'diagonal') then
                  call random(E(m1:m2,1:nrens),nra*nrens)
                  E(m1:m2,:)=sqrt(obsvar%atmos)*E(m1:m2,:)
               elseif(trim(covmodel) == 'gaussian') then
                  allocate(EEfield(nx,nrens))
                  call pseudo1D(EEfield,nx,nrens,rd,dx,nx)
                  do j=1,nrens
                  do m=m1,m2
                     E(m,j)=EEfield(obsaloc(m),j)
                  enddo
                  enddo
                  E(m1:m2,:)=sqrt(obsvar%atmos)*E(m1:m2,:)
                  deallocate(EEfield)
               else
                  stop 'invalid covariance model'
               endif
            endif
         enddo
      endif
   enddo

!   print '(a)','E'
!   print '(5f13.5)',E(1:5,1:5)
!   print '(5f13.5)',E(45:49,1:5)
!   do m=1,nrobs
!      print '(a,i5,f13.5)','E: ',m,E(m,10)
!   enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   do m=1,nrobs
!      print '(a,i5,tr1,a5,2i5,2f12.4)','measurements: ',m,obs(m)%observed,obs(m)%tloc,obs(m)%xloc,obs(m)%d,obs(m)%var
!   enddo

! ESMDA adjustment of measurement perturbations and error variance
   if (cmethod == "MDA") E=sqrt(real(nmda))*E

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Construct ensemble of measurements D=d+E
   D=DA+E

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Compute innovation D'=D-HA
!!   D=D-Y

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Compute mean(HA)
   meanY=0.0
   do j=1,nrens
   do m=1,nrobs
      meanY(m)=meanY(m)+Y(m,j)
   enddo
   enddo
   meanY=(1.0/float(nrens))*meanY

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Compute HA'=HA-mean(HA)
   do j=1,nrens
      S(:,j)=Y(:,j)-meanY(:)
   enddo

   R=matmul(E,transpose(E))/real(nrens-1)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! used in square root filters
   do m =1, nrobs
      innov(m)      = obs(m)%d- meanY(m)
   enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Scaling of matrices
   allocate(scaling(nrobs))
   do m=1,nrobs
      scaling(m)=1./sqrt(R(m,m))
      S(m,:)=scaling(m)*S(m,:)
      Y(m,:)=scaling(m)*Y(m,:)
      E(m,:)=scaling(m)*E(m,:)
      D(m,:)=scaling(m)*D(m,:)
      innov(m)=scaling(m)*innov(m)
   enddo
!   print '(a,10f10.4)','scaling: ',scaling(1:10)

   do j=1,nrobs
   do i=1,nrobs
      R(i,j)=scaling(i)*R(i,j)*scaling(j)
   enddo
   enddo
   deallocate(scaling)


end subroutine enkfprep
end module m_enkfprep
