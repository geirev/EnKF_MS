module m_obsxloc
contains
subroutine obsxloc(nr,obsloc)
   use mod_dimensions
   use mod_state
   implicit none
   integer, intent(in) :: nr
   integer, intent(inout) :: obsloc(nr)
   integer i,m
   real deltaobs

   if (nr == 0) return

   deltaobs=real(nx)/real(nr)
   if (mod(nr,2) == 0) then
      do m=1,nr
         obsloc(m)=nint(deltaobs/2.0 + real(m-1)*deltaobs)
         print '(a,i4,i7)','m xloc: ',m,obsloc(m)
      enddo
   else
      m=1
      obsloc(m)=nint(real(nx)/2.0)
      print '(a,i4,i7)','m xloc: ',m,obsloc(m)
      do i=1,(nr-1)/2
         m=m+1
         obsloc(m)=nint(real(nx)/2.0 + real(i)*deltaobs)
         print '(a,i4,i7)','m xloc: ',m,obsloc(m)
         m=m+1
         obsloc(m)=nint(real(nx)/2.0 - real(i)*deltaobs)
         print '(a,i4,i7)','m xloc: ',m,obsloc(m)
      enddo
   endif
end subroutine
end module
