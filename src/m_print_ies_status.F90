module m_print_ies_status
contains
integer function print_ies_status(fac,steplength,costf,Wold,W,nrens,iter)
   use m_frobenius
   use m_ansi_colors
   implicit none
   integer, intent(in) :: nrens
   integer, intent(in) :: iter
   real, intent(in) :: fac
   real, intent(in) :: steplength
   real, intent(in) :: costf
   real, intent(in) :: W(nrens,nrens)
   real, intent(in) :: Wold(nrens,nrens)

   character(len=8)  cfac8
   character(len=10) ftag10
   character(len=10) cstep
   character(len=13) ctag13

   write(cstep,'(f10.4)')steplength
   write(cfac8,'(f8.2)')fac
   write(ctag13,'(g13.7)')costf
   write(ftag10,'(f10.4)')frobenius(Wold-W,nrens,nrens)
   print '(tr5,a,i3,9a)','main: iter=',iter,' -> ',color(" steplength:",c_yellow),color(cstep,c_yellow),&
                                                                  color(" LMfac:",c_yellow),color(cfac8,c_yellow),&
                                                                  color(" wdiff=",c_green), color(ftag10,c_green),&
                                                                  color(" costf=",c_cyan),  color(ctag13,c_cyan)
   print_ies_status=1
end function
end module
