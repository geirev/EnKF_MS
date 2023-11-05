module m_assimilation
contains
subroutine assimilation(obs,win,win0,DA,E,D,Y,S,X,W,Wold,Yold,costf,&
                        tini,tfin,nrobs,iter,obsoloc,obsaloc,obsotimes,obsatimes,lwin,steplength)
   use mod_dimensions
   use mod_state
   use mod_observation
   use m_readinfile
   use m_pertE
   use m_prepY
   use m_scaling
   use m_ies_steplength
   use m_ies
   use m_print_ies_status
   implicit none
   integer, intent(in) :: nrobs
   integer, intent(in) :: iter
   type(state),    intent(inout) :: win(0:nrt,nrens)
   type(state),    intent(in) :: win0(0:nrt,nrens)
   type(observation),  intent(in):: obs(nrobs)
   integer, intent(in) :: lwin
   integer, intent(in) :: tini
   integer, intent(in) :: tfin
   integer, intent(in) :: obsoloc(nro)
   integer, intent(in) :: obsaloc(nra)
   integer, intent(in) :: obsotimes(nrt)
   integer, intent(in) :: obsatimes(nrt)
   real,    intent(in) :: DA(nrobs,nrens)
   real,    intent(inout) :: steplength
   integer ldw
   integer iprt
   integer j

   real, intent(inout) :: costf(nmda,nrwindows)
   real, intent(inout) :: E(nrobs,nrens)
   real, intent(inout) :: D(nrobs,nrens)
   real, intent(inout) :: Y(nrobs,nrens)
   real, intent(inout) :: S(nrobs,nrens)
   real, intent(inout) :: X(nrens,nrens)             ! The X matrix :-)
   real, intent(inout) :: W(nrens,nrens)             ! ies iteration matrix
   real, intent(inout) :: Wold(nrens,nrens)          ! ies iteration matri for reducing steplength
   real, intent(inout) :: Yold(nrobs,nrens)          ! ies predicted measurements for reducing steplength

   if ((cmethod(1:3)=='MDA') .or. ((cmethod(1:3)=='IES').and.(iter==1))) then
!     Using ESMDA: We simullate a new E for every DA step and then scale it with sqrt(nmda),
!                  thus, D=DA+sqrt(nmda)*E is updated every MDA step.
!     Using IES  : We simulate E only in the first iteration and compute D=DA+E which we use unchanged
!                  in all following IES iterations.  We are passing the predicted measurement Y and the
!                  perturbed measurements D to IES
      print '(tr5,a,i3,a)','main: iter=',iter,' -> Calling pertE to get E'
      call pertE(E,nrobs,tini,tfin,obsoloc,obsaloc,obsotimes,obsatimes)
      if (cmethod(1:3) == 'MDA') then
         E=sqrt(real(nmda))*E
         W=0.0
         steplength=1.0
      endif
      D=DA+E
      call scaling(D,E,nrobs,nrens)
   endif

   print '(tr5,a,i3,a)','main: iter=',iter,' -> Calling prepY to get Y'
   call prepY(Y,win,nrobs,tini,tfin,obsoloc,obsaloc,obsotimes,obsatimes)
   call scaling(Y,E,nrobs,nrens)

   if (cmethod(1:3) == 'IES') then
      call ies_steplength(steplength,costf,nrwindows,nmda,W,Wold,D,Y,Yold,nrens,nrobs,iter,lwin)
      if (steplength < 0.01) return
   endif

   print '(tr5,a,2(i3,a))','main: iter=',iter,' -> Calling ies with mode: ',mode_analysis
   call ies(Y,D,W,nrens,nrobs,steplength,mode_analysis)

   if (cmethod(1:3) == 'IES') iprt=print_ies_status(steplength,costf(iter,lwin),Wold,W,nrens,iter)

! X = I + W/sqrt(N-1)
   X=W/sqrt(real(nrens-1))
   do j=1,nrens
      X(j,j)=X(j,j)+1.0
   enddo

! Window update: for last iteration and in ES we are updating the whole window not just the intial condition.
! For linear dynamics, this doesn't change anything, but for strongly unstable dynamics it may give a better
! posterior estimate over the window, and better starting point for the next window.
   write(*,'(tr5,a,i3,a)',advance='no')'main: iter=',iter,' -> Ensemble update'
   if ((iter==nmda).and.(.not.lsim)) then
      ldw=ndim*(nrw+1)
      call dgemm('N','N',ldw,nrens,nrens,1.0,win0(tini:tfin,:),ldw,X,nrens,0.0,win(tini:tfin,:),ldw)
   else
      ldw=ndim
      call dgemm('N','N',ldw,nrens,nrens,1.0,win0(tini,:),ldw,X,nrens,0.0,win(tini,:),ldw)
   endif
   print '(a,i3,a)','..... done'

end subroutine
end module
