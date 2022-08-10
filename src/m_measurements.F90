module m_measurements
contains
subroutine measurements(ana,obs,obsvar,nroceanobs,nratmosobs,iobs,mkobs,time)
   use mod_dimensions
   use mod_state
   use mod_observation
   use m_pseudo1D
   use m_random
   implicit none
   type(state),    intent(in) :: ana
   integer, intent(in) :: nroceanobs
   integer, intent(in) :: nratmosobs
   type(observation),  intent(out):: obs(nroceanobs+nratmosobs)
   integer, intent(in) :: iobs
   logical, intent(in) :: mkobs
   real,    intent(in) :: time
   type(substate), intent(in) :: obsvar

   integer m,i,reclA
   real tt
   logical ex

   real tmpo(nroceanobs)
   real tmpa(nratmosobs)

   inquire(iolength=reclA)tt,i,m,obs
   inquire(file='obs.uf',exist=ex)

   if (mkobs) then
      call random(tmpo,nroceanobs)
      call random(tmpa,nratmosobs)

      do m=1,nroceanobs
         obs(m)%var=obsvar%ocean
         obs(m)%d= ana%ocean(obs(m)%pos) + sqrt(obs(m)%var)*tmpo(m)
      enddo
      i=nroceanobs
      do m=1,nratmosobs
         obs(i+m)%var=obsvar%atmos
         obs(i+m)%d= ana%atmos(obs(i+m)%pos) + sqrt(obs(i+m)%var)*tmpa(m)
      enddo

      open(10,file='obs.uf',form='unformatted',access='direct',recl=reclA)
         write(10,rec=iobs)time,nroceanobs,nratmosobs,obs
      close(10)

   else

      if (ex) then
         open(10,file='obs.uf',form='unformatted',access='direct',recl=reclA)
            read(10,rec=iobs,err=100)tt,i,m,obs
         close(10)
         if (tt /= time .or. i /= nroceanobs .or. m /= nratmosobs) then
            print *,'Problem reading obs.uf :',iobs,tt,time,i,nroceanobs,m,nratmosobs
         endif
      else
         print *,'file obs.uf does not exist'
         stop
      endif

   endif

   return
   100 stop 'Error reading obs.uf (record does not exist?)'
end subroutine measurements
end module m_measurements


