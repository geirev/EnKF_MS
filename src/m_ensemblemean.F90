module m_ensemblemean
contains
subroutine ensemblemean(A,ave,nrens)
   use mod_state
   implicit none
   integer, intent(in) :: nrens
   type(state), intent(in)  :: A(nrens)
   type(state), intent(out) :: ave
   integer j

   if (nrens == 1) then
      ave=A(1)
      return
   endif

   ave=A(1)
   do j=2,nrens
      ave=ave+A(j)
   enddo
   ave=(1.0/real(nrens))*ave

end subroutine 
end module 
