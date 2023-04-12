module m_model
   use mod_dimensions
   use mod_state
   type(substate) u                      ! Model velocity
   type(substate) rh                     ! Decorrelation lengths of ocean and atmosphere
   real, parameter :: dx=1.0             ! Horizontal grid spacing
   real, parameter :: dt=1.0             ! Time step of atmospheric model

   real :: o2a
   real :: a2o
   real :: b1
   real :: b2

contains
   subroutine model(mem,old,leuler)
   use mod_dimensions
   use mod_state
   implicit none
   type(state),     intent(inout):: mem
   type(state),     intent(inout):: old
   logical,         intent(in)   :: leuler
   integer i,ia,ib
   type(state) new
   if (leuler) then
      do i=1,nx
         ia=mod(i-2+nx,nx)+1
         ib=mod(i,nx)+1
         !new%atmos(i) = mem%atmos(i) - dt*u%atmos*(mem%atmos(ib)-mem%atmos(ia))/(2.0*dx) &
         new%atmos(i) = mem%atmos(ia) &
                      + dt*o2a*mem%ocean(i) - dt*b1*mem%atmos(i)
         new%ocean(i) = mem%ocean(i) - dt*u%ocean*(mem%ocean(ib)-mem%ocean(ia))/(2.0*dx) &
                      + dt*a2o*mem%atmos(i) - dt*b2*mem%ocean(i)
      enddo
   else
      do i=1,nx
         ia=mod(i-2+nx,nx)+1
         ib=mod(i,nx)+1
         !new%atmos(i) =  old%atmos(i) - dt*u%atmos*(mem%atmos(ib)-mem%atmos(ia))/dx &
         new%atmos(i) =  mem%atmos(ia) &
                      +2.0*dt*o2a*mem%ocean(i) - 2.0*dt*b1*mem%atmos(i)
         new%ocean(i) = old%ocean(i) - dt*u%ocean*(mem%ocean(ib)-mem%ocean(ia))/dx &
                      +2.0*dt*a2o*mem%atmos(i) - 2.0*dt*b2*mem%ocean(i)
      enddo
   endif
   old=mem
   mem=new

end subroutine
end module
