module m_windowstat
contains
subroutine windowstat(win,nrw,nrens,mean,stdt)
   use mod_dimensions
   use mod_state
   implicit none
   integer, intent(in) :: nrw
   integer, intent(in) :: nrens

   type(state), intent(in)    :: win(0:nrw,nrens)
   type(state), intent(inout) :: mean(0:nrw)
   type(state), intent(inout) :: stdt(0:nrw)

   integer j,k

   do k=0,nrw
      mean(k)=0.0
      stdt(k)=0.0
   enddo

! mean
   do j=1,nrens
   do k=0,nrw
      mean(k)=mean(k)+win(k,j)
   enddo
   enddo
   do k=0,nrw
      mean(k)=mean(k)*(1.0/real(nrens))
   enddo

! variance
   do j=1,nrens
   do k=0,nrw
      stdt(k)=stdt(k)+(win(k,j)-mean(k))*(win(k,j)-mean(k))
   enddo
   enddo
   do k=0,nrw
      stdt(k)=stdt(k)*(1.0/real(nrens-1))
      stdt(k)=sqrt(stdt(k))
   enddo

end subroutine
end module
