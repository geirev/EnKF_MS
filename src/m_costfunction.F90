module m_costfunction
contains
real function costfunction(W,H,nrens,nrobs)
   implicit none
   integer, intent(in) :: nrens
   integer, intent(in) :: nrobs
   real, intent(in)    :: W(nrens,nrens)
   real, intent(in)    :: H(nrobs,nrens)
   integer i,j

   costfunction=0.0
   do j=1,nrens
   do i=1,nrens
      costfunction=costfunction+W(i,j)**2
   enddo
   enddo

   do j=1,nrens
   do i=1,nrobs
      costfunction=costfunction+H(i,j)**2
   enddo
   enddo


end function
end module
