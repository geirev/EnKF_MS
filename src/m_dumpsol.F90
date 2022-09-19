module m_dumpsol
contains
subroutine dumpsol(time,ana,ave,var,cov,nx,dx,obs,nrobs,mem,nrens,xx)
   use mod_state
   use mod_observation
   implicit none
   integer, intent(in) :: nx
   integer, intent(in) :: nrens
   integer, intent(in) :: nrobs
   real, intent(in) :: time
   real, intent(in) :: dx
   type(state), intent(in) :: mem(nrens)
   type(state), intent(in) :: ana
   type(state), intent(in) :: ave
   type(state), intent(in) :: var
   type(state), intent(in) :: cov(2)
   type(state) :: std
   type(observation), intent(in) :: obs(nrobs)
   integer i,m
   character(len=1) xx
   character(len=4) outtag
   logical ex

   inquire(file='Solution',exist=ex)
   if ( .not.ex ) call system('mkdir Solution')

   inquire(file='Members',exist=ex)
   if ( .not.ex ) call system('mkdir Ensemble')

   write(outtag,'(i4.4)')nint(time)

   open(10,file='Solution/sol_'//trim(outtag)//xx//'.dat')
      std=sqrt(var)
      write(10,'(11a11)')'         x','  Ref_Ocean','  Ref_Atmos',&
                                      '  Est_Ocean','  Est_Atmos',&
                                      '  Std_Ocean','  Std_Atmos',&
                                      ' Cov1_Ocean',' Cov1_Atmos',&
                                      ' Cov2_Ocean',' Cov2_Atmos'
      do i=1,nx
         write(10,'(11f11.4)')real(i-1)*dx,ana%ocean(i),ana%atmos(i),&
                                           ave%ocean(i),ave%atmos(i),&
                                           std%ocean(i),std%atmos(i),&
                                           cov(1)%ocean(i),cov(1)%atmos(i),&
                                           cov(2)%ocean(i),cov(2)%atmos(i)
      enddo
   close(10)

   if (nrobs > 0) then
      open(10,file='Solution/oceanobs_'//trim(outtag)//'.dat')
      open(11,file='Solution/atmosobs_'//trim(outtag)//'.dat')
         do m=1,nrobs
            if (obs(m)%observed(1:5) == 'ocean') then
               write(10,'(2i4,2f10.4)')obs(m)%tloc,obs(m)%xloc-1,obs(m)%d,sqrt(obs(m)%var) !,sqrt(obs(m)%var)
            endif
            if (obs(m)%observed(1:5) == 'atmos') then
               write(11,'(2i4,2f10.4)')obs(m)%tloc,obs(m)%xloc-1,obs(m)%d,sqrt(obs(m)%var) !,sqrt(obs(m)%var)
            endif
         enddo
      close(10)
      close(11)
   endif

!   open(10,file='Members/ens_'//trim(outtag)//xx//'.dat')
!      do i=1,nx
!         write(10,'(202f10.4)')real(i-1)*dx,mem(1:min(10,nrens))%ocean(i),mem(1:min(10,nrens))%atmos(i)
!      enddo
!   close(10)

end subroutine dumpsol
end module m_dumpsol


