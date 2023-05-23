module m_model
   use mod_dimensions
   use mod_state
   implicit none
   type(substate) u                      ! Model velocity
   type(substate) rh
   real, parameter :: dx=1.0             ! Horizontal grid spacing
   real, parameter :: dtout=1.0             ! Time between outputs

   type modpar
      real d1         ! diffusion
      real d2         ! biharmonic diffusion
      real oa         ! coupling
      real friction   ! friction
      real vback      ! background velocity
      real vlin       ! linear state dep velocity
   end type

   type(modpar) :: alpha,omega




contains
   subroutine model(mem,old,leuler)
   use mod_dimensions
   use mod_state
   implicit none
   type(state),     intent(inout):: mem
   type(state),     intent(inout):: old
   logical,         intent(inout):: leuler
   integer i,ia,ib,ia2,ib2,nrsteps,n
   type(state) new
   real dt
   real maxu


   dt=1.0

   if (alpha%d1 /= 0.0 ) dt=min(dt, dx**2/(2.0*abs(alpha%d1)))
   if (alpha%d2 /= 0.0 ) dt=min(dt, dx**4/(4.0*abs(alpha%d2)))
   if (omega%d1 /= 0.0 ) dt=min(dt, dx**2/(2.0*abs(omega%d1)))
   if (omega%d2 /= 0.0 ) dt=min(dt, dx**4/(4.0*abs(omega%d2)))


   maxu=0.0
   do i=1,nx
      u%atmos=alpha%vback + alpha%vlin*mem%atmos(i)
      maxu=max(maxu,abs(u%atmos))
      u%ocean=omega%vback + omega%vlin*mem%ocean(i)
      maxu=max(maxu,abs(u%ocean))
   enddo
   dt=min(dt,0.25*dx/maxu)
   nrsteps=nint(1.0/dt)
   dt=1.0/real(nrsteps)
   leuler=.true.

   print *,'n=',nrsteps,dt,real(nrsteps)*dt
   do n=1,nrsteps
      if (leuler) then
         do i=1,nx
            ia=mod(i-2+nx,nx)+1
            ib=mod(i,nx)+1
            ia2=mod(i-3+nx,nx)+1
            ib2=mod(i+1,nx)+1

            u%atmos=alpha%vback + alpha%vlin*mem%atmos(i)

            if (u%atmos > 0.0) then
               new%atmos(i) = mem%atmos(i) - dt*u%atmos*(mem%atmos(i)-mem%atmos(ia))/dx
            elseif (u%atmos < 0.0) then
               new%atmos(i) = mem%atmos(i) - dt*u%atmos*(mem%atmos(ib)-mem%atmos(i))/dx
            else
               new%atmos(i) = mem%atmos(i)
            endif
            new%atmos(i) = new%atmos(i)                                                                                      &
                + dt*alpha%d1*(mem%atmos(ia)-2.0*mem%atmos(i)+mem%atmos(ib))/(dx**2)                                         &
                + dt*alpha%d2*(mem%atmos(ia2)-4.0* mem%atmos(ia)+6.0*mem%atmos(i)-4.0*mem%atmos(ib)+mem%atmos(ib2))/(dx**4)  &
                + dt*alpha%oa*(mem%ocean(i)-mem%atmos(i))                                                                    &
                - dt*alpha%friction*mem%atmos(i)

            u%ocean=omega%vback + omega%vlin*mem%ocean(i)
            if (u%ocean > 0.0) then
               new%ocean(i) = mem%ocean(i) - dt*u%ocean*(mem%ocean(i)-mem%ocean(ia))/dx
            elseif (u%ocean < 0.0) then
               new%ocean(i) = mem%ocean(i) - dt*u%ocean*(mem%ocean(ib)-mem%ocean(i))/dx
            else
               new%ocean(i) = mem%ocean(i)
            endif
            new%ocean(i) = new%ocean(i)                                                                                      &
                + dt*omega%d1*(mem%ocean(ia)-2.0*mem%ocean(i)+mem%ocean(ib))/(dx**2)                                         &
                + dt*omega%d2*(mem%ocean(ia2)-4.0* mem%ocean(ia)+6.0*mem%ocean(i)-4.0*mem%ocean(ib)+mem%ocean(ib2))/(dx**4)  &
                + dt*omega%oa*(mem%atmos(i)-mem%ocean(i))                                                                    &
                - dt*omega%friction*mem%ocean(i)
         enddo
      else
         do i=1,nx
            ia=mod(i-2+nx,nx)+1
            ib=mod(i,nx)+1
            ia2=mod(i-3+nx,nx)+1
            ib2=mod(i+1,nx)+1

            u%atmos=alpha%vback + alpha%vlin*mem%atmos(i)
            new%atmos(i) = old%atmos(i) - dt*u%atmos*(mem%atmos(ib)-mem%atmos(ia))/dx                                           &
                + 2.0*dt*alpha%d1*(mem%atmos(ia)-2.0*mem%atmos(i)+mem%atmos(ib))/(dx**2)                                        &
                + 2.0*dt*alpha%d2*(mem%atmos(ia2)-4.0*mem%atmos(ia)+6.0*mem%atmos(i)-4.0*mem%atmos(ib)+mem%atmos(ib2))/(dx**4) &
                + 2.0*dt*alpha%oa*(mem%ocean(i)-mem%atmos(i))                                                                   &
                - 2.0*dt*alpha%friction*mem%atmos(i)

            u%ocean=omega%vback + omega%vlin*mem%ocean(i)
            new%ocean(i) = old%ocean(i) - dt*u%ocean*(mem%ocean(ib)-mem%ocean(ia))/dx                                           &
                + 2.0*dt*omega%d1*(mem%ocean(ia)-2.0*mem%ocean(i)+mem%ocean(ib))/(dx**2)                                        &
                + 2.0*dt*omega%d2*(mem%ocean(ia2)-4.0*mem%ocean(ia)+6.0*mem%ocean(i)-4.0*mem%ocean(ib)+mem%ocean(ib2))/(dx**4) &
                + 2.0*dt*omega%oa*(mem%atmos(i)-mem%ocean(i))                                                                   &
                - 2.0*dt*omega%friction*mem%ocean(i)
         enddo
      endif
      old=mem
      mem=new
   enddo

end subroutine
end module
