module m_dumpsol
contains
subroutine dumpsol(time,ana,ave,var,cov,nx,dx,obs,nroceanobs,nratmosobs,mem,nrens,xx)
   use mod_state
   use mod_observation
   implicit none
   integer, intent(in) :: nx
   integer, intent(in) :: nrens
   integer, intent(in) :: nroceanobs
   integer, intent(in) :: nratmosobs
   real, intent(in) :: time
   real, intent(in) :: dx
   type(state), intent(in) :: mem(nrens)
   type(state), intent(in) :: ana
   type(state), intent(in) :: ave
   type(state), intent(in) :: var
   type(state), intent(in) :: cov(2)
   type(state) :: std
   type(observation), intent(in) :: obs(nroceanobs+nratmosobs)
   integer i,m
   character(len=1) xx
   character(len=4) outtag
   logical ex


!   inquire(iolength=reclA)time,ana,ave,var
!   open(10,file='solution.uf',form='unformatted',access='direct',recl=reclA)
!      write(10,rec=iprt)time,ana,ave,sqrt(var+0.000001)
!   close(10)


   inquire(file='Solution',exist=ex)
   if ( .not.ex ) call system('mkdir Solution')

   inquire(file='Members',exist=ex)
   if ( .not.ex ) call system('mkdir Ensemble')

   write(outtag,'(i4.4)')nint(time)

   open(10,file='Solution/sol_'//trim(outtag)//xx//'.dat')
      std=sqrt(var)
      write(10,'(11a10)')'          x','  Ref_Ocean','  Ref_Atmos',&
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

   if (time > 0.0) then
      open(10,file='Solution/oceanobs_'//trim(outtag)//'.dat')
         write(10,'(3a10)')'         x',' Ocean_obs','   Std_obs'
         do m=1,nroceanobs
            write(10,'(4f10.4)')(obs(m)%pos-1)*dx,obs(m)%d,sqrt(obs(m)%var)
         enddo
      close(10)

      open(10,file='Solution/atmosobs_'//trim(outtag)//'.dat')
         i=nroceanobs
         write(10,'(3a10)')'         x',' Atmos_obs','   Std_obs'
         do m=1,nratmosobs
            write(10,'(5f10.4)')(obs(i+m)%pos-1)*dx,obs(i+m)%d,-sqrt(obs(i+m)%var),sqrt(obs(i+m)%var)
         enddo
      close(10)
   endif

   open(10,file='Members/ens_'//trim(outtag)//xx//'.dat')
      do i=1,nx
         write(10,'(202f10.4)')real(i-1)*dx,mem(1:min(10,nrens))%ocean(i),mem(1:min(10,nrens))%atmos(i)
      enddo
   close(10)

end subroutine dumpsol
end module m_dumpsol


