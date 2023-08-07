program main
   use mod_dimensions                      ! Defines state dimension
   use mod_state                           ! Defines model state
   use mod_observation                     ! Defines observation state
   use m_readinfile                        ! Reading infile.in
   use m_set_random_seed2                  ! Sets new random seed unless 'seed.dat' exists
   use m_model                             ! Model time stepping
   use m_frobenius                         ! Frobenius norm between two matrices
   use m_pseudo1D                          ! Samples pseudo-random fields in 1-D
   use m_fixsample1D                       ! Ensures zero mean and unit variance of sampled ensembles
   use m_obsxloc                           ! Computes x position of measurements
   use m_obstloc                           ! Computes time of measurements
   use m_obscount                          ! Counts the number of measurements in a DA window
   use m_obspoints                         ! Saves the measurement locations in space and time for plotting
   use m_ensemblemean                      ! Compute ensemble mean at current time (from mem)
   use m_ensemblevariance                  ! Compute ensemble variance at current time (from mem)
   use m_ensemblecovariance                ! Compute ensemble covariance at current time (from mem)
   use m_enkfprep                          ! Setting up predicted measurements and matrices for analysis
   use m_prepY                             ! Setting up predicted measurements for IES analysis
   use m_scaling                           ! Setting up predicted measurements for IES analysis
   use m_prepD                             ! Setting up predicted measurements for IES analysis
   use m_pertE                             ! Setting up predicted measurements for IES analysis
   use m_random                            ! Generate random normal numbers
   use m_tecfld                            ! Tecplot output (not used)
   use mod_shapiro                         ! Shapiro filter in case it is needed
   use m_windowstat                        ! Compute ensemble mean and variance in a DA window (from win)
   use m_covstat                           ! Computes space time covariances in case you have a lot of memory
   use m_gnuplot                           ! Generate space time plots for plotting in gnuplot (c.gnu)
   use m_ies                               ! The ies analysis routine
   use m_print_ies_status                  ! Prints status of ies convergence
   use m_ies_steplength                    ! Computes the new step lengths for ies
   use, intrinsic :: omp_lib
   use m_ansi_colors
   implicit none

   real, allocatable        :: tmp1(:,:),tmp2(:,:)
   type(state), allocatable :: full(:,:)   ! The ensemble of realizations over the whole simulation
   type(state), allocatable :: win(:,:)    ! The ensemble of realizations over an assimilation window.
   type(state), allocatable :: win0(:,:)   ! Initial ensemble of realizations over an assimilation window.
   type(state), allocatable :: winref(:)   ! The reference solution over an assimilation window.
   type(state), allocatable :: mem(:)      ! The ensemble of realizations
   type(state), allocatable :: sysnoise(:) ! The ensemble of sampled system noise updated every timestep
   type(state) ref                         ! reference solution
   type(state) ave                         ! ensemble average

   real, allocatable :: samples(:,:)       ! work array used when sampling in pseudo1D

! Spacew time statistics diagnostic variables
   type(state), allocatable :: refout(:)   ! ensemble average as a function of space and time
   type(state), allocatable :: mean(:)     ! ensemble average as a function of space and time
   type(state), allocatable :: stdt(:)     ! ensemble std dev as a function of space and time
   type(state), allocatable :: covo(:)     ! ensemble std dev as a function of space and time
   type(state), allocatable :: cova(:)     ! ensemble std dev as a function of space and time
   type(substate), allocatable :: rmse(:)  ! time series of rmse values (mean - referece)
   type(substate), allocatable :: rmss(:)  ! time series of rms std values

! Observation location variables
   integer, allocatable :: obsoloc(:)      ! location of ocean observations in space
   integer, allocatable :: obsaloc(:)      ! location of atmos observations in space
   integer, allocatable :: obsotimes(:)    ! location of ocean observations in time
   integer, allocatable :: obsatimes(:)    ! location of atmos observations in time

   type(observation), allocatable :: obsA(:)
   type(observation), allocatable :: obsO(:)

! EnKF analysis variables
   real, allocatable :: YO(:,:)             ! Predicted meaurements
   real, allocatable :: EO(:,:)             ! Measurement perturbations
   real, allocatable :: DDO(:,:)             ! Innovations
   real, allocatable :: DAO(:,:)            ! Innovations
   real, allocatable :: WO(:,:)             ! ies iteration matrix
   real, allocatable :: WoldO(:,:)          ! ies iteration matri for reducing steplength
   real, allocatable :: YoldO(:,:)          ! ies predicted measurements for reducing steplength
   real, allocatable :: XO(:,:)            ! The X matrix :-)
   real, allocatable :: YA(:,:)             ! Predicted meaurements
   real, allocatable :: EA(:,:)             ! Measurement perturbations
   real, allocatable :: DDA(:,:)             ! Innovations
   real, allocatable :: DAA(:,:)            ! Innovations
   real, allocatable :: WA(:,:)             ! ies iteration matrix
   real, allocatable :: WoldA(:,:)          ! ies iteration matri for reducing steplength
   real, allocatable :: YoldA(:,:)          ! ies predicted measurements for reducing steplength
   real, allocatable :: XA(:,:)            ! The X matrix :-)


! other variables
   integer tini                            ! start time of an assimilation window
   integer tfin                            ! end   time of an assimilation window
   integer nrobs                           ! Total number of measurement per assimilation window
   integer nrobsO
   integer nrobsA

   integer j,k,l,iter,ldw,iprt,i,n
   real time
   real fac
   real :: dxsamp=1.0
   real :: steplength
   real, allocatable :: costf(:,:)

   call readinfile()

   call set_random_seed2

   call system('rm -f eigenvalues.dat')
   call system('rm -f obsloc?.dat')
!   call system('touch obsloca.dat obsloco.dat')

! Now that we know all dimensions, allocate the main arrays
   if (lglobstat) allocate (full(0:nrt,nrens))
   allocate (win(0:nrw,nrens))
   allocate (win0(0:nrw,nrens))
   allocate (winref(0:nrw))
   allocate (refout(0:nrt))
   allocate (mean(0:nrt))
   allocate (stdt(0:nrt))
   allocate (covo(0:nrt))
   allocate (cova(0:nrt))
   allocate (mem(nrens))
   allocate (sysnoise(nrens))
   allocate (samples(nx,nrens))
   allocate (rmse(0:nrt))
   allocate (rmss(0:nrt))

   allocate (WO(nrens,nrens))
   allocate (WA(nrens,nrens))
   allocate (XO(nrens,nrens))
   allocate (XA(nrens,nrens))
   allocate (WoldO(nrens,nrens))
   allocate (WoldA(nrens,nrens))
   allocate (YoldO(nrens,nrens))
   allocate (YoldA(nrens,nrens))

   allocate (obsoloc(nro))
   allocate (obsaloc(nra))
   allocate (obsotimes(nrt))
   allocate (obsatimes(nrt))

   allocate(costf(nmda,nrwindows))
   costf=0.0
   time=0.0

!   call omp_set_num_threads(10)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Uniformly distributed measurements for ocean and atmosphere in space and time
   print *,'Ocean observation locations'
   call obsxloc(nro,obsoloc)
   call obstloc(nrt,obst0o,obsdto,obsotimes)
   call obspoints(trim(outdir)//'/obsloco',obsotimes,nrt,obsoloc,nro)

   print *,'Atmos observation locations'
   call obsxloc(nra,obsaloc)
   call obstloc(nrt,obst0a,obsdta,obsatimes)
   call obspoints(trim(outdir)//'/obsloca',obsatimes,nrt,obsaloc,nra)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! The true solution smooth pseudo random field drawn from  N(0,1,rh).
   call pseudo1D(ref%ocean,nx,1,rh%ocean,dxsamp,nx)
   call pseudo1D(ref%atmos,nx,1,rh%atmos,dxsamp,nx)
   print *,'main: ref ok'

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! First guess solution stored in ave
   call pseudo1D(ave%ocean,nx,1,rh%ocean,dxsamp,nx)
   call pseudo1D(ave%atmos,nx,1,rh%atmos,dxsamp,nx)
   ave=(ave + ref)*(1.0/sqrt(2.0))
   print *,'main: fg ok'

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Initialization of ensemble
   call pseudo1D(samples,nx,nrens,rh%ocean,dxsamp,nx)
   if (samp_fix) call fixsample1D(samples,nx,nrens)
   do j=1,nrens
      mem(j)%ocean(1:nx)=samples(1:nx,j)
   enddo

   call pseudo1D(samples,nx,nrens,rh%atmos,dxsamp,nx)
   if (samp_fix) call fixsample1D(samples,nx,nrens)
   do j=1,nrens
      mem(j)%atmos(1:nx)=samples(1:nx,j)
   enddo

   do j=1,nrens
      mem(j)%ocean=ave%ocean + sqrt(inivar%ocean)*mem(j)%ocean
      mem(j)%atmos=ave%atmos + sqrt(inivar%atmos)*mem(j)%atmos
   enddo
   print *,'main: ensemble ok'


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Loop over assimilation windows
   print '(a)','---------------------------------------------------------------------------'
   do l=1,nrwindows                                   ! Integrate model over the assimilation window
      tini=(l-1)*nrw                                  ! time at beginning of assimilation window
      tfin=l*nrw                                      ! time at end of assimilation window
      print *
      print '(a)','---------------------------------------------------------------------------'
      print '(a,i3.3,a,i5.5,a,i5.5)','DA window number ',l,' running from ',tini,' to ',tfin
      print *

! Reference solution over assimilation window stored in refout
      winref(0)=ref
      do k=1,nrw
         call model(ref)
         winref(k)=ref
      enddo
      refout(tini:tfin)=winref

! prior ensemble integration over assimilation window
      win(0,:)=mem(:)                                 ! store ensemble at tini in win
      print '(tr5,a,i3,a,f10.2,a,i3)','main: prior ensemble simulation:  window=',l
!$OMP PARALLEL DO
      do j=1,nrens
         do k=1,nrw
            call model(mem(j))
            win(k,j)=mem(j)
         enddo
      enddo
!$OMP END PARALLEL DO
      win0=win ! Store prior ensemble over window for IES algorithm as win=win0*X

! Counting the number of measurements in the current DA window
      nrobs=obscount(nrt,tini,tfin,obsotimes,obsatimes,nro,nra)
      nrobsA=obscount(nrt,tini,tfin,obsotimes,obsatimes,0,nra)
      nrobsO=obscount(nrt,tini,tfin,obsotimes,obsatimes,nro,0)
      print '(tr5,a,3i5)','main: Total number of measurements in the DA window= ',nrobs,nrobsA,nrobsO
      print *

      if (nrobs > 0) then
! Allocate DA variables
         if (nrobsO > 0) then
            allocate(obsO(nrobsO))
            allocate(EO(nrobsO,nrens))
            allocate(DDO(nrobsO,nrens))
            allocate(DAO(nrobsO,nrens))
            allocate(YO(nrobsO,nrens))
         endif
         if (nrobsA > 0) then
            allocate(obsA(nrobsA))
            allocate(EA(nrobsA,nrens))
            allocate(DDA(nrobsA,nrens))
            allocate(DAA(nrobsA,nrens))
            allocate(YA(nrobsA,nrens))
         endif

! Generate the observations in obs and DA
         winref=refout(tini:tfin)
         print '(tr5,a)','main: -> Calling prepD to get DA'
         if (nrobsO > 0) call prepD(obsO,DAO,winref,nrobsO,l,tini,tfin,obsoloc,obsaloc,obsotimes,obsatimes,iter,&
                                    nro,0,nrt,nrw,nrens,obsvar)
         if (nrobsA > 0) call prepD(obsA,DAA,winref,nrobsA,l,tini,tfin,obsoloc,obsaloc,obsotimes,obsatimes,iter,&
                                    0,nra,nrt,nrw,nrens,obsvar)

! Loop over nmda steps
         steplength=steplength0
         WO=0.0
         WA=0.0
         fac=1.0
         do iter=1,nmda
            print '(tr5,a,i3,a)','main: iter=',iter,' -> '
            if ((cmethod(1:3)=='MDA') .or. ((cmethod(1:3)=='IES').and.(iter==1))) then
!                 Using ESMDA: We simullate a new E for every DA step and then scale it with sqrt(nmda),
!                              thus, D=DDA+sqrt(nmda)*E is updated every MDA step.
!                 Using IES  : We simulate E only in the first iteration and compute D=DDA+E which we use unchanged
!                              in all following IES iterations.  We are passing the predicted measurement Y and the
!                              perturbed measurements D to IES
               print '(tr5,a,i3,a)','main: iter=',iter,' -> Calling pertE'
               if (nrobsO > 0) call pertE(EO,nrobsO,l   ,tini,tfin,obsoloc,obsaloc,obsotimes,obsatimes,nro,0,nrt,nrens,&
                                          obsvar,covmodel,dx,rd,nrw)
               if (nrobsA > 0) call pertE(EA,nrobsA,l   ,tini,tfin,obsoloc,obsaloc,obsotimes,obsatimes,0,nra,nrt,nrens,&
                                          obsvar,covmodel,dx,rd,nrw)
               if (cmethod(1:3) == 'MDA') then
                  if (nrobsO > 0) EO=sqrt(real(nmda))*EO
                  if (nrobsA > 0) EA=sqrt(real(nmda))*EA
                  if (nrobsO > 0) WO=0.0
                  if (nrobsA > 0) WA=0.0
                  steplength=1.0
               endif
            print '(tr5,a,i3,a)','main: iter=',iter,' -> Calling D scaling'
               if (nrobsO > 0) DDO=DAO+EO
               if (nrobsA > 0) DDA=DAA+EA
               if (nrobsO > 0) call scaling(DDO,EO,nrobsO,nrens)
               if (nrobsA > 0) call scaling(DDA,EA,nrobsA,nrens)
            endif

            if ((cmethod(1:3)=='IES').and.LM) then
            ! Levenberg–Marquardt algorithm with IES
               fac=10.0+100.0*exp(-10.0*real(iter-1)/(real(nmda)))
            endif

            print '(tr5,a,i3,a)','main: iter=',iter,' -> Calling prepY'
            if (nrobsO > 0) call prepY(YO,win,nrobsO,l,tini,tfin,obsoloc,obsaloc,obsotimes,obsatimes,nro,0,nrt,nrw,nrens)
            if (nrobsA > 0) call prepY(YA,win,nrobsA,l,tini,tfin,obsoloc,obsaloc,obsotimes,obsatimes,0,nra,nrt,nrw,nrens)
            if (nrobsO > 0) call scaling(YO,EO,nrobsO,nrens)
            if (nrobsA > 0) call scaling(YA,EA,nrobsA,nrens)

            if (cmethod(1:3) == 'IES') then
               print '(tr5,a,i3,a)','main: iter=',iter,' -> Calling steplength'
               if (nrobsO > 0) call ies_steplength(steplength,costf,nrwindows,nmda,WO,WoldO,DDO,YO,YoldO,nrens,nrobsO,iter,l)
               if (nrobsA > 0) call ies_steplength(steplength,costf,nrwindows,nmda,WA,WoldA,DDA,YA,YoldA,nrens,nrobsA,iter,l)
               if (steplength < 0.01) exit
            endif

            print '(tr5,a,i3,a,i3)','main: iter=',iter,' -> Calling iesO with mode: ',mode_analysis
            if (nrobsO > 0) call ies(YO,DDO,WO,nrens,nrobsO,steplength,mode_analysis,fac)
            print '(tr5,a,i3,a,i3)','main: iter=',iter,' -> Calling iesA with mode: ',mode_analysis
            if (nrobsA > 0) call ies(YA,DDA,WA,nrens,nrobsA,steplength,mode_analysis,fac)

            print '(tr5,a,i3,a)','main: iter=',iter,' -> Calling print_ies_status'
            if (cmethod(1:3) == 'IES') iprt=print_ies_status(fac,steplength,costf(iter,l),WoldO,WO,nrens,iter)

! X = I + W/sqrt(N-1)
            if (nrobsO > 0) then
               XO=WO/sqrt(real(nrens-1))
               do j=1,nrens
                  XO(j,j)=XO(j,j)+1.0
               enddo
            endif
            if (nrobsA > 0) then
               XA=WA/sqrt(real(nrens-1))
               do j=1,nrens
                  XA(j,j)=XA(j,j)+1.0
               enddo
            endif

! Window update: for last iteration and in ES we are updating the whole window not just the intial condition.
! For linear dynamics, this doesn't change anything, but for strongly unstable dynamics it may give a better
! posterior estimate over the window, and better starting point for the next window.
            write(*,'(tr5,a,i3,a)',advance='no')'main: iter=',iter,' -> Ensemble update'
            if ((iter==nmda).and.(.not.lsim)) then
               print *,'YY'
               ldw=nx*(nrw+1)
               allocate( tmp1(ldw,nrens), tmp2(ldw,nrens) )
! Ocean
               if (nrobsO > 0) then
                  do j=1,nrens
                     do k=0,nrw
                        do i=1,nx
                           n=k*nx+i
                           tmp1(n,j)=win0(k,j)%Ocean(i)
                        enddo
                     enddo
                  enddo

                  call dgemm('N','N',ldw,nrens,nrens,1.0,tmp1,ldw,XO,nrens,0.0,tmp2,ldw)

                  do j=1,nrens
                     do k=0,nrw
                        do i=1,nx
                           n=k*nx+i
                           win(k,j)%Ocean(i)=tmp2(n,j)
                        enddo
                     enddo
                  enddo
               endif
! Atmos
               if (nrobsA > 0) then
                  do j=1,nrens
                     do k=0,nrw
                        do i=1,nx
                           n=k*nx+i
                           tmp1(n,j)=win0(k,j)%Atmos(i)
                        enddo
                     enddo
                  enddo

                  call dgemm('N','N',ldw,nrens,nrens,1.0,tmp1,ldw,XA,nrens,0.0,tmp2,ldw)

                  do j=1,nrens
                     do k=0,nrw
                        do i=1,nx
                           n=k*nx+i
                           win(k,j)%Atmos(i)=tmp2(n,j)
                        enddo
                     enddo
                  enddo
               endif

               deallocate(tmp1,tmp2)


            else
               allocate( tmp1(nx,nrens), tmp2(nx,nrens) )
! Ocean
               if (nrobsO > 0) then
                  do j=1,nrens
                     do i=1,nx
                        tmp1(i,j)=win0(0,j)%Ocean(i)
                     enddo
                  enddo

                  call dgemm('N','N',nx,nrens,nrens,1.0,tmp1,nx,XO,nrens,0.0,tmp2,nx)

                  do j=1,nrens
                     do i=1,nx
                        win(0,j)%Ocean(i)=tmp2(i,j)
                     enddo
                  enddo
               endif
! Atmos
               if (nrobsA > 0) then
                  do j=1,nrens
                     do i=1,nx
                        tmp1(i,j)=win0(0,j)%Atmos(i)
                     enddo
                  enddo

                  call dgemm('N','N',nx,nrens,nrens,1.0,tmp1,nx,XA,nrens,0.0,tmp2,nx)

                  do j=1,nrens
                     do i=1,nx
                        win(0,j)%Atmos(i)=tmp2(i,j)
                     enddo
                  enddo
               endif

               deallocate(tmp1,tmp2)

            endif
            print '(a,i3,a)','..... done'

! timestep loop over assimilation window
            if ((iter < nmda).or.lsim) then
               write(*,'(tr5,a,i3,a,i3,i5,a,i5)',advance='no')'main: iter=',iter,' -> post-update-RESIM: window=',l,tini,'->',tfin
               mem(:)=win(0,:)
!$OMP PARALLEL DO
               do j=1,nrens
                  do k=1,nrw
                     call model(mem(j))
                     win(k,j)=mem(j)
                  enddo
               enddo
!$OMP END PARALLEL DO
               print '(a)','..... done'
               print *
            else
               write(*,'(tr5,a,i3,a,i3,i5,a,i5)',advance='no')'main: iter=',iter,' -> post-update-NOSIM: window=',l,tini,'->',tfin
            endif

! update prior ensemble in MDA
            if (cmethod(1:3)=="MDA") win0=win

         enddo ! end loop over mda steps
            print '(a)','..... p'
         if (nrobsA > 0) deallocate(obsA,YA,EA,DDA,DAA)
            print '(a)','..... pp'
         if (nrobsO > 0) deallocate(obsO,YO,EO,DDO,DAO)
            print '(a)','..... ppp'
      endif

! Continue integration in next assimilation window from final value of current window
      mem(:)=win(nrw,:)

! Only needed to compute space-time covariances
      if (lglobstat) full(tini:tfin,:)=win(0:nrw,:)

! Compute ensemble mean and variance over current DA window
      call windowstat(win,nrw,nrens,mean(tini:tfin),stdt(tini:tfin))
   enddo ! End loop over assimilation windows

! Compute root-mean-squared errors and std. dev.
   do k=0,nrt
      rmse(k)=sqrt(average((mean(k)-refout(k))*(mean(k)-refout(k))))
      rmss(k)=sqrt(average(stdt(k)*stdt(k)))
   enddo
   open(10,file=trim(outdir)//'/rmse.dat')
      do k=0,nrt
         write(10,'(i5,1024g12.4)')k,rmse(k)%Atmos,rmse(k)%Ocean,rmss(k)%Atmos,rmss(k)%Ocean
      enddo
   close(10)

! Print the cost functions values for each window.
   if (cmethod == 'IES') then
      open(10,file=trim(outdir)//'/costf.dat')
         do iter=1,nmda
            write(10,'(i5,100g12.4)')iter,costf(iter,1:nrwindows)
         enddo
      close(10)
   endif

! Dumping reference solution, mean and standard deviations for plotting
   call gnuplot('gnu_ref',refout,nrt,outdir)
   call gnuplot('gnu_ave',mean,nrt,outdir)
   call gnuplot('gnu_std',stdt,nrt,outdir)
   if (lglobstat) then
      print *,'calling full covariance statistics'
      call covstat(full,nrt,nrens,mean,stdt,covo,cova,outdir)
      call gnuplot('gnu_covo',covo,nrt,outdir)
      call gnuplot('gnu_cova',cova,nrt,outdir)
   endif

end program main

