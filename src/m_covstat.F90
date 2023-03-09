module m_covstat
contains
subroutine covstat(full,nrt,nrens,mean,stdt,covo,cova,outdir)
   use mod_dimensions
   use mod_state
   implicit none
   integer, intent(in) :: nrt
   integer, intent(in) :: nrens
   character(len=25), intent(in) :: outdir

   type(state), intent(in)    :: full(0:nrt,nrens)
   type(state), intent(inout) :: mean(0:nrt)
   type(state), intent(inout) :: stdt(0:nrt)
   type(state), intent(inout) :: covo(0:nrt)
   type(state), intent(inout) :: cova(0:nrt)

   integer j,k,ic,kc

   ic=nx/2
   kc=nrt/2

   open(10,file=trim(outdir)//'/'//'obs.dat')
      write(10,'(2I5)')ic,kc
   close(10)

   do k=0,nrt
      covo(k)=0.0
      cova(k)=0.0
   enddo

! covariance function
   do j=1,nrens
   do k=0,nrt
      covo(k)=covo(k)+(full(k,j)-mean(k))*(full(kc,j)%ocean(ic)-mean(kc)%ocean(ic))
      cova(k)=cova(k)+(full(k,j)-mean(k))*(full(kc,j)%atmos(ic)-mean(kc)%atmos(ic))
   enddo
   enddo
   do k=0,nrt
      covo(k)=covo(k)*(1.0/real(nrens-1))
      cova(k)=cova(k)*(1.0/real(nrens-1))
   enddo

! correlation function
   do k=0,nrt
      covo(k)=covo(k)/(stdt(kc)%ocean(ic)*stdt(k))
      cova(k)=cova(k)/(stdt(kc)%atmos(ic)*stdt(k))
   enddo


end subroutine
end module
