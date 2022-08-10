pi = 3.1415927
a1 = 0.5
a2 = 0.1
b1 = 0.
b2 = 0.0

    if (first == 0) then
    do i=1,nx
       ia=mod(i-2+nx,nx)+1
       ib=mod(i,nx)+1
       new%atmos(i) = mem%atmos(i) - dt*u%atmos*(mem%atmos(ib)-mem%atmos(ia))/2. &
                    + dt*a1*mem%ocean(i) - dt*b1*mem%atmos(i)
       new%ocean(i) = mem%ocean(i) - dt*u%ocean*(mem%ocean(ib)-mem%ocean(ia))/2. &
                    + dt*a2*mem%atmos(i) - dt*b2*mem%ocean(i)
    enddo
    else
   do i=1,nx
      ia=mod(i-2+nx,nx)+1
      ib=mod(i,nx)+1
      new%atmos(i) = old%atmos(i) - dt*u%atmos*(mem%atmos(ib)-mem%atmos(ia))/dx +2.0*dt*a1*mem%ocean(i) - dt*b1*mem%atmos(i)
      new%ocean(i) = old%ocean(i) - dt*u%ocean*(mem%ocean(ib)-mem%ocean(ia))/2. +2.0*dt*a2*mem%atmos(i) - dt*b2*mem%ocean(i)
   enddo

   endif

