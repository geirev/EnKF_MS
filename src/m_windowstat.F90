module m_windowstat
contains
subroutine windowstat(win,nrt,nrens,mean,stdt)
   use mod_dimensions
   use mod_state
   implicit none
   integer, intent(in) :: nrt
   integer, intent(in) :: nrens

   type(state), intent(in)    :: win(0:nrt,nrens)
   type(state), intent(inout) :: mean(0:nrt)
   type(state), intent(inout) :: stdt(0:nrt)

   integer j,k

   do k=0,nrt
      mean(k)=0.0
      stdt(k)=0.0
   enddo

   if (nrens==1) then
      do k=0,nrt
         mean(k)=win(k,1)
      enddo
      return
   endif

! mean
   do j=1,nrens
   do k=0,nrt
      mean(k)=mean(k)+win(k,j)
   enddo
   enddo
   do k=0,nrt
      mean(k)=mean(k)*(1.0/real(nrens))
   enddo

! variance
   do j=1,nrens
   do k=0,nrt
      stdt(k)=stdt(k)+(win(k,j)-mean(k))*(win(k,j)-mean(k))
   enddo
   enddo
   do k=0,nrt
      stdt(k)=stdt(k)*(1.0/real(nrens-1))
      stdt(k)=sqrt(stdt(k))
   enddo

end subroutine
end module
