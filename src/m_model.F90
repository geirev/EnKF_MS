module m_model
contains
   subroutine model(mem,old,u,dx,dt,leuler)
   use mod_dimensions
   use mod_state
   implicit none
   type(state),     intent(inout):: mem
   type(state),     intent(inout):: old
   type(substate),  intent(in)   :: u
   real,            intent(in)   :: dx
   real,            intent(in)   :: dt
   logical,         intent(in)   :: leuler
   integer i,ia,ib
   type(state) new
   real, parameter :: o2a=0.01
   real, parameter :: a2o=0.00
   real, parameter :: b1=0.001
   real, parameter :: b2=0.001
   if (leuler) then
      do i=1,nx
         ia=mod(i-2+nx,nx)+1
         ib=mod(i,nx)+1
         !new%atmos(i) = mem%atmos(i) - dt*u%atmos*(mem%atmos(ib)-mem%atmos(ia))/2. &
         new%atmos(i) = mem%atmos(ia) &
                      + dt*o2a*mem%ocean(i) - dt*b1*mem%atmos(i)
         new%ocean(i) = mem%ocean(i) - dt*u%ocean*(mem%ocean(ib)-mem%ocean(ia))/2. &
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
