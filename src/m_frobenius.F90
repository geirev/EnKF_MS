module m_frobenius
implicit none
contains
real function frobenius(A,m,n)
   integer, intent(in) :: n
   integer, intent(in) :: m
   real, intent(in) :: A(m,n)
   integer i,j

   frobenius=0.0
   do j=1,n
   do i=1,m
      frobenius=frobenius+A(i,j)**2
   enddo
   enddo
   frobenius=sqrt(frobenius)

end function
end module
