module m_obspoints
contains
subroutine obspoints(fname,obstimes,nrt,obsloc,nr)
   implicit none
   character(len=*), intent(in) :: fname
   integer, intent(in) :: nrt
   integer, intent(in) :: nr
   integer, intent(in) :: obstimes(nrt)
   integer, intent(in) :: obsloc(nr)
   integer i,k

      open(10,file=trim(fname)//'.dat')
         do k=1,nrt
         if (obstimes(k) == 0) exit
         do i=1,nr
            write(10,'(2i5)')obsloc(i),obstimes(k)
         enddo
         enddo
      close(10)
end subroutine
end module
