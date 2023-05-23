module m_statistics
use mod_dimensions
use mod_state
real, allocatable :: fullave(:,:,:)
real, allocatable :: fullvar(:,:,:)
real, allocatable :: fullcor(:,:,:)
real, allocatable :: fullcoroa(:,:,:)
contains
subroutine statistics(full,nrt,nrens,nrf)
   implicit none
   integer, intent(in) :: nrt
   integer, intent(in) :: nrens
   integer, intent(in) :: nrf
   type(state), intent(in) :: full(0:nrt,nrens)
   integer i,j,k,ic,kc,stat

   allocate(fullave(nx,0:nrt,2), stat=stat) ; if (stat /= 0) stop "Memory allocation failed for fullave"
   allocate(fullvar(nx,0:nrt,2), stat=stat) ; if (stat /= 0) stop "Memory allocation failed for fullave"
   allocate(fullcor(nx,0:nrt,2), stat=stat) ; if (stat /= 0) stop "Memory allocation failed for fullave"
   allocate(fullcoroa(nx,0:nrt,2), stat=stat) ; if (stat /= 0) stop "Memory allocation failed for fullave"

   fullave = 0.0
   fullvar = 0.0
   fullcor = 0.0
   fullcoroa = 0.0


! mean
   do j=1,nrens
   do k=0,nrt
   do i=1,nx
      fullave(i,k,1)=fullave(i,k,1)+full(k,j)%ocean(i)
      fullave(i,k,2)=fullave(i,k,2)+full(k,j)%atmos(i)
   enddo
   enddo
   enddo
   fullave=fullave/real(nrens)

! variance
   do j=1,nrens
   do k=0,nrt
   do i=1,nx
      fullvar(i,k,1)=fullvar(i,k,1)+(full(k,j)%ocean(i)-fullave(i,k,1))**2
      fullvar(i,k,2)=fullvar(i,k,2)+(full(k,j)%atmos(i)-fullave(i,k,2))**2
   enddo
   enddo
   enddo
   fullvar=sqrt(fullvar/real(nrens-1))

! corariance oo and aa
   ic=nx/2
   kc=nrt/2
   do j=1,nrens
   do k=0,nrt
   do i=1,nx
      fullcor(i,k,1)=fullcor(i,k,1) + (full(k,j)%ocean(i)-fullave(i,k,1)) * (full(kc,j)%ocean(ic)-fullave(ic,kc,1))
      fullcor(i,k,2)=fullcor(i,k,2) + (full(k,j)%atmos(i)-fullave(i,k,2)) * (full(kc,j)%atmos(ic)-fullave(ic,kc,2))
   enddo
   enddo
   enddo
   fullcor=fullcor/real(nrens-1)
   do k=0,nrt
   do i=1,nx
      fullcor(i,k,1)=fullcor(i,k,1)/(fullvar(i,k,1)*fullvar(ic,kc,1))
      fullcor(i,k,2)=fullcor(i,k,2)/(fullvar(i,k,2)*fullvar(ic,kc,2))
   enddo
   enddo

! corariance oa and ao
   do j=1,nrens
   do k=0,nrt
   do i=1,nx
      fullcoroa(i,k,1)=fullcoroa(i,k,1) + (full(k,j)%atmos(i)-fullave(i,k,2)) * (full(kc,j)%ocean(ic)-fullave(ic,kc,1))
      fullcoroa(i,k,2)=fullcoroa(i,k,2) + (full(k,j)%ocean(i)-fullave(i,k,1)) * (full(kc,j)%atmos(ic)-fullave(ic,kc,2))
   enddo
   enddo
   enddo
   fullcoroa=fullcor/real(nrens-1)
   do k=0,nrt
   do i=1,nx
      fullcoroa(i,k,1)=fullcoroa(i,k,1)/(fullvar(i,k,2)*fullvar(ic,kc,1))
      fullcoroa(i,k,2)=fullcoroa(i,k,2)/(fullvar(i,k,1)*fullvar(ic,kc,2))
   enddo
   enddo



end subroutine
end module
