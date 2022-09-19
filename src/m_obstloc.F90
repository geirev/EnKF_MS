module m_obstloc
contains
subroutine obstloc(nrt,t0,obsdt,obstimes)
   implicit none
   integer, intent(in)  :: nrt
   integer, intent(in)  :: t0
   integer, intent(in)  :: obsdt
   integer, intent(out) :: obstimes(nrt)
   integer i,tt
   obstimes(:)=0
   do i=1,nrt
      tt=t0+(i-1)*obsdt
      if (tt > nrt) exit
      obstimes(i)=tt
      print *,'obs times: ',obstimes(i)
   enddo
end subroutine
end module
