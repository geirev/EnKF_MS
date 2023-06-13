module m_scaling
contains
subroutine  scaling(X,E,m,n)
   integer, intent(in)    :: m
   integer, intent(in)    :: n
   real,    intent(in)    :: E(m,n)
   real,    intent(inout) :: X(m,n)
   real fac(m)
   integer i,j

   do i=1,m
      fac=sqrt(dot_product(E(i,:),E(i,:))/real(n-1))
      X(i,:)=X(i,:)/fac(i)
!      print '(a,i3,f10.3)','aa',i,fac(i)
   enddo

   return

   fac=0.0
   do j=1,n
      do i=1,m
         fac(i)=fac(i)+E(i,j)*E(i,j)
      enddo
   enddo
   do i=1,m
      fac(i)=sqrt(fac(i)/real(n-1))
!      print '(a,i3,f10.3)','bb',i,fac(i)
   enddo
   do j=1,n
      do i=1,m
         X(i,j)=X(i,j)/fac(i)
      enddo
   enddo

end subroutine
end module
