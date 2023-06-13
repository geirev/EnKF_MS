module m_prepY
contains
subroutine prepY(Y,win,nrobs,lwin,tini,tfin,obsoloc,obsaloc,obsotimes,obsatimes)
   use mod_dimensions
   use mod_state
   use mod_observation
   use m_pseudo1D
   use m_random
   use m_readinfile
   use m_model
   implicit none
   type(state),    intent(in) :: win(0:nrw,nrens)
   integer, intent(in) :: nrobs
   integer, intent(in) :: lwin
   integer, intent(in) :: tini
   integer, intent(in) :: tfin
   integer, intent(in) :: obsoloc(nro)
   integer, intent(in) :: obsaloc(nra)
   integer, intent(in) :: obsotimes(nrt)
   integer, intent(in) :: obsatimes(nrt)
   real,  intent(out):: Y(nrobs,nrens)

   integer m,i,k,mm,m1,m2



   m=0
   do i=1,nrt
      if (tini < obsotimes(i) .and. obsotimes(i) <= tfin .and. nro > 0) then
         do k=1,nrw
            if (obsotimes(i) == (lwin-1)*nrw+k ) then
               m1=m+1
               m2=m+nro
               do mm=1,nro
                  m=m+1
                  Y(m,:) = win(k,:)%ocean(obsoloc(mm))
               enddo
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
                  Y(m,:) = win(k,:)%atmos(obsaloc(mm))
               enddo
            endif
         enddo
      endif
   enddo

end subroutine prepY
end module m_prepY

