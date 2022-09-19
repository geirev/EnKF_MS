module m_gnuplot
contains
subroutine gnuplot(fname,variable,nrt)
   use mod_dimensions
   use mod_state
   character(len=*), intent(in) :: fname
   integer, intent(in) :: nrt
   type(state), intent(in) :: variable(0:nrt)
   integer k,i

   print *,trim(fname)//'o.dat'
   open(10,file=trim(fname)//'o.dat')
      write(10,'(1025i5)')nx,(i,i=1,nx)
      do k=0,nrt
         write(10,'(i5,1024g12.4)')k,variable(k)%ocean(:)
      enddo
   close(10)

   print *,trim(fname)//'a.dat'
   open(10,file=trim(fname)//'a.dat')
      write(10,'(1025i5)')nx,(i,i=1,nx)
      do k=0,nrt
         write(10,'(i5,1024g12.4)')k,variable(k)%atmos(:)
      enddo
   close(10)

end subroutine
end module
