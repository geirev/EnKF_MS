module m_prepD
contains
subroutine prepD(obs,D,winref,nrobs,lwin,tini,tfin,obsoloc,obsaloc,obsotimes,obsatimes,iter,nro,nra,nrt,nrw,nrens,obsvar)
   use mod_dimensions
   use mod_state
   use mod_observation
   use m_pseudo1D
   use m_random
!   use m_readinfile
   use m_model
   implicit none
   integer, intent(in) :: nrens
   integer, intent(in) :: nrw
   integer, intent(in) :: nro
   integer, intent(in) :: nra
   integer, intent(in) :: nrt
   integer, intent(in) :: nrobs
   integer, intent(in) :: lwin
   integer, intent(in) :: tini
   integer, intent(in) :: tfin
   integer, intent(in) :: obsoloc(nro)
   integer, intent(in) :: obsaloc(nra)
   integer, intent(in) :: obsotimes(nrt)
   integer, intent(in) :: obsatimes(nrt)
   integer, intent(in) :: iter
   type(observation),  intent(out):: obs(nrobs)
   real,  intent(out):: D(nrobs,nrens)
   type(substate), intent(in) :: obsvar
   type(state), intent(in) :: winref(0:nrw)

   integer m,i,k,mm,m1,m2,j

   real tmpo(nro)
   real tmpa(nra)


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
                  obs(m)%xloc=obsoloc(mm)
                  obs(m)%tloc=obsotimes(i)
                  obs(m)%var=obsvar%ocean
                  obs(m)%observed='ocean'
                  obs(m)%d=winref(k)%ocean(obs(m)%xloc) + sqrt(obsvar%ocean)*tmpo(mm)
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
                  obs(m)%xloc=obsaloc(mm)
                  obs(m)%tloc=obsatimes(i)
                  obs(m)%var=obsvar%atmos
                  obs(m)%observed='atmos'
                  obs(m)%d=winref(k)%atmos(obs(m)%xloc) + sqrt(obsvar%atmos)*tmpa(mm)
               enddo
            endif
         enddo
      endif
   enddo

! Construct ensemble of measurements D=d+E
   do j=1,nrens
      do m=1,nrobs
         D(m,j)=obs(m)%d
      enddo
   enddo
!   do m=1,nrobs
!      print '(a,i5,f13.5)','obs(m)%d: ',m,obs(m)%d
!   enddo
!   stop


end subroutine prepD
end module m_prepD

