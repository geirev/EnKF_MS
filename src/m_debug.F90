module m_debug
contains
subroutine debug(mem,ref,nrens,outdir,lwin,iter)
   use mod_dimensions
   use mod_state
   use mod_observation
   implicit none
   integer, intent(in) :: nrens
   type(state), intent(in) :: mem(nrens)
   type(state), intent(in) :: ref
   integer, intent(in) :: iter
   integer, intent(in) :: lwin
   integer i,j
   type(state) ave
   character(len=25), intent(in) :: outdir
   character(len=2) tagiter, tagwin

   ave=0.0
   do j=1,nrens
      ave=ave + mem(j)
   enddo
   ave=ave*(1.0/real(nrens))

   write(tagwin,'(i2.2)')lwin
   write(tagiter,'(i2.2)')iter
   open(10,file=trim(outdir)//'/sol'//tagwin//'_'//tagiter//'.dat')
      write(10,'(22a11)')'         i','  Ref_ocean','  Ref_atmos',&
                                      '  Est_ocean','  Est_atmos',&
                                      '  mem1ocean','  mem1atmos',&
                                      '  mem2ocean','  mem2atmos',&
                                      '  mem3ocean','  mem3atmos',&
                                      '  mem4ocean','  mem4atmos',&
                                      '  mem5ocean','  mem5atmos',&
                                      '  mem6ocean','  mem6atmos',&
                                      '  mem7ocean','  mem7atmos',&
                                      '  mem8ocean','  mem8atmos',&
                                      '  mem9ocean','  mem9atmos'
      do i=1,nx
         write(10,'(i11,22f11.4)')i,ref%ocean(i),ref%atmos(i),&
                                   ave%ocean(i),ave%atmos(i),&
                                   mem(1)%ocean(i),mem(1)%atmos(i),&
                                   mem(2)%ocean(i),mem(2)%atmos(i),&
                                   mem(3)%ocean(i),mem(3)%atmos(i),&
                                   mem(4)%ocean(i),mem(4)%atmos(i),&
                                   mem(5)%ocean(i),mem(5)%atmos(i),&
                                   mem(6)%ocean(i),mem(6)%atmos(i),&
                                   mem(7)%ocean(i),mem(7)%atmos(i),&
                                   mem(8)%ocean(i),mem(8)%atmos(i),&
                                   mem(9)%ocean(i),mem(9)%atmos(i)
      enddo
   close(10)

end subroutine debug
end module m_debug


