module m_obscount
contains
integer function  obscount(nrt,tini,tfin,obsotimes,obsatimes,nro,nra)
   integer, intent(in) :: nrt
   integer, intent(in) :: tini
   integer, intent(in) :: tfin
   integer, intent(in) :: obsotimes(nrt)
   integer, intent(in) :: obsatimes(nrt)
   integer, intent(in) :: nro
   integer, intent(in) :: nra
   integer m,i

   m=0
   do i=1,nrt
!     print *,'obscount: ',tini,obsotimes(i),obsatimes(i),tfin, nro, nra
      if (tini < obsotimes(i) .and. obsotimes(i) <= tfin) then
         m=m+nro
      endif
      if (tini < obsatimes(i) .and. obsatimes(i) <= tfin) then
         m=m+nra
      endif
      if (tfin < obsotimes(i) .and. tfin < obsatimes(i)) exit
   enddo
   obscount=m
end function
end module
