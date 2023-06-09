module m_model
   use mod_dimensions
   use mod_state
   implicit none
   type(substate) u                      ! Model velocity
   type(substate) rh
   real :: lxo                           ! Length of ocean model domain
   real :: lxa                           ! Length of atmos model domain
   real, parameter :: dtout=1.0             ! Time between outputs
   real, parameter :: pi=3.14159265359

   integer advection

   type modpar
      real d1         ! diffusion
      real d2         ! biharmonic diffusion
      real oa         ! coupling
      real vback      ! background velocity
      real v          ! linear state dep velocity
   end type

   type(modpar) :: alpha,omega

contains

subroutine model(mem)
   use mod_dimensions
   use mod_state
   use mod_fftw3

   implicit none
   type(state),     intent(inout):: mem
   integer i,nrsteps,n
   real dt
   integer, save :: ifirst=1

   integer(8), save :: plan_c2r, plan_r2c
   real, save :: kappaa(0:nx/2)
   real, save :: kappao(0:nx/2)
   complex, save :: ddxa(0:nx/2)
   complex, save :: ddxo(0:nx/2)
   complex, save :: Ga(0:nx/2)
   complex, save :: Go(0:nx/2)
   real, save :: difa(0:nx/2)
   real, save :: difo(0:nx/2)
   real, save :: Aa(0:nx/2)
   real, save :: Ao(0:nx/2)
   real, save :: Ba(0:nx/2)
   real, save :: Bo(0:nx/2)
   real, save :: dt2,dt32

   real atmos(0:nx)
   real ocean(0:nx)
   complex  atmosfft(0:nx/2)
   complex  oceanfft(0:nx/2)

   real temp(0:nx)
   complex  tempfft(0:nx/2)

   complex :: nna(0:nx/2)
   complex :: nn1a(0:nx/2)
   complex :: nno(0:nx/2)
   complex :: nn1o(0:nx/2)

   complex :: lina(0:nx/2)
   complex :: lin1a(0:nx/2)
   complex :: lino(0:nx/2)
   complex :: lin1o(0:nx/2)


   dt=1.0/16.0
   nrsteps=nint(dtout/dt)
   dt=1.0/real(nrsteps)

   !print '(a,11f15.5)','A atmos:',mem%atmos(1:10),mem%atmos(nx)

   if (ifirst==1) then
      ifirst=0.0
!     Setting up the fftw3 plans
      call dfftw_plan_dft_r2c_1d(plan_r2c, nx, atmos, atmosfft, FFTW_MEASURE)
      call dfftw_plan_dft_c2r_1d(plan_c2r, nx, atmosfft, atmos, FFTW_MEASURE)

      dt2=dt/2.0
      dt32=3.0*dt/2.0

      do i=0,nx/2
         kappaa(i) =2.0*pi*real(i)/real(Lxa)                      ! Derivative of wave numbers kappa(k)=2*pi*k/(nx*dx)
         kappao(i) =2.0*pi*real(i)/real(Lxo)                      ! Derivative of wave numbers kappa(k)=2*pi*k/(nx*dx)
      enddo
      kappaa(nx/2)=0.0
      kappao(nx/2)=0.0

      do i=0,nx/2
         ddxa(i)=cmplx( 0.0, kappaa(i) )                          ! Spectral D=d/dx operator  img*kappa(i)
         ddxo(i)=cmplx( 0.0, kappao(i) )                          ! Spectral D=d/dx operator  img*kappa(i)
         Ga(i)   =-0.5*ddxa(i)                                    ! -0.5*D
         Go(i)   =-0.5*ddxo(i)                                    ! -0.5*D

         difa(i) =-alpha%d1*kappaa(i)**2 + alpha%d2*kappaa(i)**4    ! Diffusion operator in wave space
         difo(i) =-omega%d1*kappao(i)**2 + omega%d2*kappao(i)**4    ! Diffusion operator in wave space
         Aa(i)=1.0 + dt2*difa(i)
         Ao(i)=1.0 + dt2*difo(i)
         Ba(i)=1.0/(1.0 - dt2*difa(i))
         Bo(i)=1.0/(1.0 - dt2*difo(i))
      enddo
   endif

   atmos(1:nx)=mem%atmos(1:nx)
   atmos(0)=atmos(nx)
   call dfftw_execute_dft_r2c(plan_r2c, atmos, atmosfft)

   ocean(1:nx)=mem%ocean(1:nx)
   ocean(0)=ocean(nx)
   call dfftw_execute_dft_r2c(plan_r2c, ocean, oceanfft)

   lina(:)=(0.0,0.0)
   lino(:)=(0.0,0.0)
   nna(:)=(0.0,0.0)
   nno(:)=(0.0,0.0)


   do n=1,nrsteps
      nn1a(:)=nna(:)
      nn1o(:)=nno(:)
      lin1a(:)=lina(:)
      lin1o(:)=lino(:)

! Linear background advection
      lina(0:nx/2)=alpha%vback*Ga(0:nx/2)*atmosfft(0:nx/2)
      if (n==1) lin1a=lina

      lino(0:nx/2)=omega%vback*Go(0:nx/2)*oceanfft(0:nx/2)
      if (n==1) lin1o=lino

! Nonlinear advection
      tempfft=atmosfft
      call dfftw_execute_dft_c2r(plan_c2r, tempfft, temp)
      temp(:)=temp(:)/real(nx)
      temp(:)=temp(:)*temp(:)
      call dfftw_execute_dft_r2c(plan_r2c, temp, tempfft)
      do i=0,nx/2
         nna(i)=alpha%v*Ga(i)*tempfft(i)
      enddo
      if (n==1) nn1a=nna

      tempfft=oceanfft
      call dfftw_execute_dft_c2r(plan_c2r, tempfft, temp)
      temp(:)=temp(:)/real(nx)
      temp(:)=temp(:)*temp(:)
      call dfftw_execute_dft_r2c(plan_r2c, temp, tempfft)
      do i=0,nx/2
         nno(i)=omega%v*Go(i)*tempfft(i)
      enddo
      if (n==1) nn1o=nno

      do i=0,nx/2
         tempfft(i)=oceanfft(i)-atmosfft(i)
         atmosfft(i)=Ba(i)*(Aa(i)*atmosfft(i) + dt32*(nna(i)+lina(i)) - dt2*(nn1a(i)+lin1a(i)) + alpha%oa*tempfft(i) )
         oceanfft(i)=Bo(i)*(Ao(i)*oceanfft(i) + dt32*(nno(i)+lino(i)) - dt2*(nn1o(i)+lin1o(i)) - omega%oa*tempfft(i) )
      enddo

   enddo

   call dfftw_execute_dft_c2r(plan_c2r, atmosfft, atmos)
   call dfftw_execute_dft_c2r(plan_c2r, oceanfft, ocean)
   atmos(nx)=atmos(0)
   ocean(nx)=ocean(0)
   atmos=atmos/real(nx)
   ocean=ocean/real(nx)
   mem%atmos(1:nx)=atmos(1:nx)
   mem%ocean(1:nx)=ocean(1:nx)
!   print '(a,11f15.5)','E atmos:',mem%atmos(1:10),mem%atmos(nx)
!   call dfftw_destroy_plan(plan_r2c)
!   call dfftw_destroy_plan(plan_c2r)

end subroutine
end module
