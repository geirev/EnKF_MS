module m_ensemblevariance
contains
subroutine ensemblevariance(A,ave,var,nrens)
   use mod_state
   implicit none
   integer, intent(in) :: nrens
   type(state), intent(in)    :: A(nrens)
   type(state), intent(in)    :: ave
   type(state), intent(inout) :: var
   integer j

   var=0.0
   do j=1,nrens
      var=var+(A(j)-ave)*(A(j)-ave)
   enddo
   var=(1.0/real(nrens-1))*var

end subroutine ensemblevariance
end module m_ensemblevariance
