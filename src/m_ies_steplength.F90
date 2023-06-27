module m_ies_steplength
contains
subroutine ies_steplength(steplength,costf,nrwindows,nmda,W,Wold,D,Y,Yold,nrens,nrobs,iter,l)
   use m_costfunction
   use m_ansi_colors
   implicit none
   real, intent(inout) :: steplength
   integer, intent(in) :: nrwindows
   integer, intent(in) :: nmda
   real, intent(inout) :: costf(nmda,nrwindows)
   integer, intent(in) :: nrobs
   integer, intent(in) :: nrens
   real, intent(inout) :: W(nrens,nrens)
   real, intent(inout) :: Wold(nrens,nrens)
   real, intent(inout) :: D(nrobs,nrens)
   real, intent(inout) :: Y(nrobs,nrens)
   real, intent(inout) :: Yold(nrobs,nrens)
   integer, intent(in) :: iter
   integer, intent(in) :: l
   character(len=10) cstep

   costf(iter,l)=costfunction(W,D-Y,nrens,nrobs)
   if ((iter > 1) .and. (costf(iter,l) > costf(iter-1,l))) then
      steplength=0.5*steplength
      W=WOld
      Y=YOld
      write(cstep,'(f10.4)')steplength
      print '(tr5,2a)',color(" New steplength is :",c_red),color(cstep,c_red)
   endif
   Wold=W
   Yold=Y
   costf(iter,l)=costfunction(W,D-Y,nrens,nrobs)

end subroutine
end module
