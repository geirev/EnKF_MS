module m_ies
contains
subroutine ies(Y,D,W,nrens,nrobs,steplength,mode_analysis,fac)
   use mod_anafunc       ! from EnKF_analysis
   use m_scaling         ! to rescale matrices to measurement std dev
   implicit none
   integer, intent(in) :: nrens            ! Number of realizations
   integer, intent(in) :: nrobs            ! Number of measurements
   real, intent(in)    :: Y(nrobs,nrens)   ! Ensemble of predicted measurements
   real, intent(in)    :: D(nrobs,nrens)   ! Ensemble of perturbed measurements D=d+E
   real, intent(inout) :: W(nrens,nrens)   ! Coefficient matrix
   real, intent(in)    :: steplength       ! Steplength in update
   integer, intent(in) :: mode_analysis    !
   real, intent(in)    :: fac              ! LM factor

   integer i
   real :: DB(nrobs,nrens)     ! D work array
   real :: YB(nrobs,nrens)     ! Y work array
   real :: E(nrobs,nrens)      ! E
   real :: S(nrobs,nrens)      ! Y*PI
   real :: H(nrobs,nrens)      ! D-Y
   real :: Omega(nrens,nrens)  ! I+W*PI
   real :: X5(nrens,nrens)

   real, allocatable :: Z(:,:)
   real, allocatable :: X3(:,:)
   real, allocatable :: eig(:)

   real :: truncation=0.99

   integer ipiv(nrens),info
   integer nrmin

   YB=Y
   if (fac /= 1.0) then
      E=proj(D,nrobs,nrens)*sqrt(real(nrens-1))
      DB=D+(fac-1.0)*E
      E=fac*E
      call scaling(YB,E,nrobs,nrens)
      call scaling(DB,E,nrobs,nrens)
   endif
   DB=D

! E = D*PI
   E=proj(DB,nrobs,nrens)

! Y = Y*PI
   S=proj(YB,nrobs,nrens)

! Remember AA projection for n < N-1

! Omega = I + W*PI
   Omega=proj(W,nrens,nrens)
   do i=1,nrens
      Omega(i,i)=Omega(i,i)+1.0
   enddo

! Solve Omega^T*S^T = Y^T for S
  call dgesv(nrens,nrobs,transpose(Omega),nrens,ipiv,transpose(S),nrens,info)

! H = S*W + D - Y
   H=DB-YB

! H=H+matmul(S,W)
   call dgemm('N','N',nrobs,nrens,nrens,1.0,S,nrobs,W,nrens,1.0,H,nrobs)

   if (mode_analysis == 10) then
      S=sqrt(real(nrens-1))*S
!      print '(a)','IES S: '
!      print '(10f13.4)',S(1:10,1:10)
!      print '(a)','IES H: '
!      print '(10f13.4)',H(1:10,1:10)
      call exact_diag_inversion(S,H,X5,nrens,nrobs)
!      print '(a)','IES X5 from exact_diag_inversion: '
!      print '(10f13.4)',X5(1:10,1:10)

   elseif (mode_analysis == 13) then
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! X5=S'*(S*S'+E*E')^+*H
      nrmin=min(nrobs,nrens)
      allocate(eig(nrmin))
      allocate(Z(nrobs,nrmin))
      allocate(X3(nrobs,nrens))

      call lowrankE(S,E,nrobs,nrens,nrmin,Z,eig,truncation,1)

      if (nrobs > 1) then
         call genX3(nrens,nrobs,nrmin,eig,Z,H,X3)
      else
         X3=H*eig(1)
      endif

      call dgemm('t','n',nrens,nrens,nrobs,1.0,S,nrobs,X3,nrobs,0.0,X5,nrens)
      deallocate(Z,X3,eig)
   else
      print '(a,i3)','invalid mode_analysis=',mode_analysis
      stop
   endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   W = W - steplength*(W - X5)
   W = (1.0-steplength)*W + steplength*X5
!   print '(a)','IES W: '
!   print '(10f13.4)',W(1:10,1:10)

end subroutine

function proj(X,n,nrens)
   implicit none
   integer, intent(in) :: n
   integer, intent(in) :: nrens
   real proj(n,nrens)
   real, intent(in) :: X(n,nrens)

   integer i
   real mean

   do i=1,n
      mean=sum(X(i,1:nrens))/real(nrens)
      proj(i,:)=X(i,:)-mean
   enddo
   proj=proj/sqrt(real(nrens-1))

end function

end module
