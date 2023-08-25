module m_pertE
contains
subroutine pertE(E,nrobs,tini,tfin,obsoloc,obsaloc,obsotimes,obsatimes)
   use mod_dimensions
   use mod_state
   use mod_observation
   use m_pseudo1D
   use m_random
   use m_readinfile
   use m_model
   implicit none
   integer, intent(in) :: nrobs
   integer, intent(in) :: tini
   integer, intent(in) :: tfin
   integer, intent(in) :: obsoloc(nro)
   integer, intent(in) :: obsaloc(nra)
   integer, intent(in) :: obsotimes(nrt)
   integer, intent(in) :: obsatimes(nrt)
   real,    intent(out):: E(nrobs,nrens)

   integer m,i,k,m1,m2,j,mm

   real, allocatable :: EEfield(:,:)

   m=0
   do i=1,nrt
      if (tini < obsotimes(i) .and. obsotimes(i) <= tfin .and. nro > 0) then
         do k=tini+1,tfin
            if (obsotimes(i) ==  k ) then
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
         do k=tini+1,tfin
            if (obsatimes(i) == k ) then

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
!   stop

end subroutine pertE
end module m_pertE

