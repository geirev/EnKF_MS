module m_readinfile
! Variables read from infile
   use mod_state
   use mod_shapiro
   use m_model

   integer nrt                           ! Number of timesteps
   integer nrw                           ! Number of timesteps within a data assimilation window
   integer nrwindows                     ! Number of data assimilation windows


   type(substate) inivar
   type(substate) sysvar
   type(substate) obsvar

   integer obst0o                        ! time of first ocean observations
   integer obst0a                        ! time of first atmos observations
   integer obsdto                        ! time between ocean observations
   integer obsdta                        ! time between atmos observations
   integer nro                           ! Number of ocean measurement per assimilation time
   integer nra                           ! Number of atmos measurement per assimilation time

   integer nrens                         ! ensemble size
   integer nmda                          ! Number of mda steps (1=ES)
   real    steplength0                    ! IES steplength
   character(len=3) cmethod              ! MDA or IES
   integer mode_analysis                 ! 1 standard, 2 fixed R
   logical samp_fix
   logical :: lupdate_randrot=.true.
   logical :: lrandrot=.true.            ! random rotation in SQRT schemes
   logical :: lsymsqrt=.true.             ! Symetrical square root in SQRT schemes

! inflation
   integer inflate                       ! 0--no inflation, 1--constant inflation, 2--adaptive inflation
   real infmult                          ! constant inflation or adjustment of adaptive inflation

! local analysis
   integer local                         ! 0-no localization, 1-distance based, 2-adaptive
   real obs_radius                       ! Number of grid cells for including measurements (distance based)
   real obs_truncation                   ! Correlations for truncating measurements in adaptive scheme

! parameters for analysis scheme
   real truncation                       ! Truncation of singular values
   character(len=8) covmodel             ! Diagonal or Gaussian measurement error covariance model
   real rd                               ! Horizontal correlation of observation errors in Gaussian case
   logical Rexact                        ! Use exact(true) or lowrank(false) R matrix
   real dx

   integer nsh
   real, allocatable :: sh(:)

! Diagnostics
   logical lglobstat                     ! Dump global statistics output
   character(len=25) outdir
   logical oldana

   logical :: lm                         ! Levenberg-Marquardt in IES algorithm
   logical :: lsim                       ! Ensemble simulation over window in final MDA or IES step

   contains
   subroutine readinfile
   use mod_dimensions
   use mod_state
   use m_model
   implicit none

   character(len=2) ca
   character(len=3) :: version='1.0'
   character(len=3) ver
   logical ex
   character(len=9) :: cmd='mkdir -p '
   character(len=50) :: cpinfile
   character(len=50) :: cpgnu
   integer i

! reading input data
   inquire(file='infile.in',exist=ex)
   if (.not.ex) then
      print '(a)','Did not find inputfile infile.in...'
      stop
   endif
   print '(a)','--------------------------------------------------------------------------------'
   open(10,file='infile.in')
      read(10,'(tr1,a)')ver
      if (version /= ver) then
         print *,'Update version of infile.in to:',version,ver
         stop
      endif
      read(10,*)outdir             ; print '(a,tr1,a)', 'outdir=     ',trim(outdir)
      read(10,*)nrens              ; print '(a,tr3,i5)','nrens=         ',nrens
      read(10,*)nrt                ; print '(a,tr3,i5)','nrt=           ',nrt
      read(10,*)nrw                ; print '(a,tr3,i5)','nrw=           ',nrw
      nrwindows=nint(real(nrt-50)/real(nrw))+1
      if (50+(nrwindows-1)*nrw < nrt) nrwindows=nrwindows+1
                                     print '(a,tr3,i5)','nrwindows=     ',nrwindows
      nrt=50+(nrwindows-1)*nrw            ; print '(a,tr3,i5)','nrt updated=   ',nrt


      read(10,'(1x,l1)')lglobstat  ; print '(a,tr7,l1)','lglobstat=     ',lglobstat
      print '(a)','--------------------------------------------------------------------------------'

      read(10,'(a)')ca
      if (ca /= '#1') then
         print *,'#1: error in infile.in'
         stop
      endif
      read(10,*)rh%ocean           ; print '(a,tr2,f6.2)', 'rh%ocean=      ',rh%ocean
      read(10,*)rh%atmos           ; print '(a,tr2,f6.2)', 'rh%atmos=      ',rh%atmos
      read(10,*)nsh                ; print '(a,tr3,I5)',   'nshapiro=      ',nsh

      if (nsh >0) then
         allocate(sh(nsh))
         call shfact(nsh,sh)
      endif
      print '(a)','--------------------------------------------------------------------------------'

      read(10,'(a)')ca
      if (ca /= '#2') then
         print *,'#2: error in infile.in'
         stop
      endif
      read(10,'(1x,l1)')samp_fix   ; print '(a,tr7,l1)',   'samp_fix=      ',samp_fix
      read(10,*)inivar%ocean       ; print '(a,tr2,f6.2)','inistd%ocean=  ',inivar%ocean ; inivar%ocean=inivar%ocean**2
      read(10,*)inivar%atmos       ; print '(a,tr2,f6.2)','inistd%atmos=  ',inivar%atmos ; inivar%atmos=inivar%atmos**2
      read(10,*)sysvar%ocean       ; print '(a,tr2,f6.2)','sysstd%ocean=  ',sysvar%ocean ; sysvar%ocean=sysvar%ocean**2
      read(10,*)sysvar%atmos       ; print '(a,tr2,f6.2)','sysstd%atmos=  ',sysvar%atmos ; sysvar%atmos=sysvar%atmos**2
      read(10,*)obsvar%ocean       ; print '(a,tr2,f6.2)','obsstd%ocean=  ',obsvar%ocean ; obsvar%ocean=obsvar%ocean**2
      read(10,*)obsvar%atmos       ; print '(a,tr2,f6.2)','obsstd%atmos=  ',obsvar%atmos ; obsvar%atmos=obsvar%atmos**2
      print '(a)','--------------------------------------------------------------------------------'

      read(10,'(a)')ca
      if (ca /= '#3') then
         print *,'#3: error in infile.in'
         stop
      endif

      read(10,*)nro                ; print '(a,tr1,i5)','nro=             ',nro
      read(10,*)nra                ; print '(a,tr1,i5)','nra=             ',nra
      read(10,*)obst0o             ; print '(a,tr1,i5)','obst0o=          ',obst0o
      read(10,*)obst0a             ; print '(a,tr1,i5)','obst0a=          ',obst0a
      read(10,*)obsdto             ; print '(a,tr1,i5)','obsdto=          ',obsdto
      read(10,*)obsdta             ; print '(a,tr1,i5)','obsdta=          ',obsdta
      print '(a)','--------------------------------------------------------------------------------'

      read(10,'(a)')ca
      if (ca /= '#4') then
         print *,'#4: error in infile.in'
         stop
      endif
      read(10,*)oldana             ; print '(a,tr9,l1)',   'oldana=      ',oldana
      read(10,*)cmethod            ; print '(a,tr7,a)',    'cmethod=     ',cmethod
      read(10,*)nmda               ; print '(a,tr7,i3 )',  'nmda=        ',nmda
      read(10,*)steplength0        ; print '(a,tr4,f6.2)', 'steplength=  ',steplength0
      read(10,*)mode_analysis      ; print '(a,tr8,i2)',   'mode_ana=    ',mode_analysis
      read(10,'(1x,l1)')lm         ; print '(a,tr9,l1)',   'Leveberg M=  ',lm
      read(10,'(1x,l1)')lsim       ; print '(a,tr9,l1)',   'lsim      =  ',lsim
      read(10,*)truncation         ; print '(a,tr4,f6.3)', 'truncation=  ',truncation
      read(10,'(1x,a)')covmodel    ; print '(a,tr2,a)',    'covmodel=    ',trim(covmodel)
      read(10,*)rd                 ; print '(a,tr4,f6.2)', 'rd=          ',rd
      read(10,'(1x,l1)')Rexact     ; print '(a,tr9,l1)',   'Rexact=      ',Rexact
      read(10,'(1x,l1)')lrandrot   ; print '(a,tr9,l1)',   'lrandrot=    ',lrandrot
      read(10,*)inflate,infmult    ; print '(a,tr5,i5,tr2,f6.3)','inflation=   ',inflate,infmult
      read(10,*)local,obs_radius,obs_truncation; print '(a,tr5,i5,tr2,2f6.3)','localization=',local,obs_radius,obs_truncation
      print '(a)','--------------------------------------------------------------------------------'
      call execute_command_line (cmd//outdir, exitstat=i)

      read(10,'(a)')ca
      if (ca /= '#5') then
         print *,'#5: error in infile.in'
         stop
      endif

      read(10,*)lxa                    ; print '(a,tr2,f8.2)','Length of ocean=',lxa
      read(10,*)alpha%d1               ; print '(a,tr2,f8.4)','alpha%d1       =',alpha%d1
      read(10,*)alpha%d2               ; print '(a,tr2,f8.4)','alpha%d2       =',alpha%d2
      read(10,*)alpha%oa               ; print '(a,tr2,f8.4)','alpha%oa       =',alpha%oa
      read(10,*)alpha%vback            ; print '(a,tr2,f8.4)','alpha%vback    =',alpha%vback
      read(10,*)alpha%v                ; print '(a,tr2,f8.4)','alpha%v        =',alpha%v

      read(10,*)lxo                    ; print '(a,tr2,f8.2)','Length of ocean=',lxo
      read(10,*)omega%d1               ; print '(a,tr2,f8.4)','omega%d1       =',omega%d1
      read(10,*)omega%d2               ; print '(a,tr2,f8.4)','omega%d2       =',omega%d2
      read(10,*)omega%oa               ; print '(a,tr2,f8.4)','omega%oa       =',omega%oa
      read(10,*)omega%vback            ; print '(a,tr2,f8.4)','omega%vback    =',omega%vback
      read(10,*)omega%v                ; print '(a,tr2,f8.4)','omega%v        =',omega%v
      print '(a)','--------------------------------------------------------------------------------'

      omega%vback=omega%vback*real(lxo)/real(lxa)
      print *,'omega%vback corrected by factor (lxo/lxa): ',real(lxo)/real(lxa)
   close(10)
   print '(a)','--------------------------------------------------------------------------------'

   cpinfile='cp infile.in infile.'//trim(outdir)
   call execute_command_line (trim(cpinfile), exitstat=i)

   cpinfile='cp infile.in '//trim(outdir)
   call execute_command_line (trim(cpinfile), exitstat=i)

!   cpgnu='./p.sh > '//trim(outdir)//'/p.gnu'
!   call execute_command_line (trim(cpgnu), exitstat=i)

   cpgnu='cp c.gnu '//trim(outdir)
   call execute_command_line (trim(cpgnu), exitstat=i)

   cpgnu='cp cpdf.gnu '//trim(outdir)
   call execute_command_line (trim(cpgnu), exitstat=i)

! We assume constant atmospheric velocity equal to one
   u%atmos=1.0

   dx=lxa/nx
end subroutine
end module
