module m_readinfile
! Variables read from infile
   use mod_state
   use mod_shapiro

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

! Shapiro filter variables
   integer nsh
   real, allocatable :: sh(:)

! Diagnostics
   logical lglobstat                     ! Dump global statistics output

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
      read(10,*)nrens              ; print *,'nrens=       ',nrens
      read(10,*)nrt                ; print *,'nrt=         ',nrt
      read(10,*)nrw                ; print *,'nrw=         ',nrw
      nrwindows=nint(real(nrt)/real(nrw))
      if (nrwindows*nrw < nrt) nrwindows=nrwindows+1
      print *,'nrwindows=   ',nrwindows
      nrt=nrwindows*nrw            ; print *,'nrt updated= ',nrt


      read(10,'(1x,l1)')lglobstat  ; print *,'lglobstat=   ',lglobstat

      read(10,'(a)')ca
      if (ca /= '#1') then
         print *,'#1: error in infile.in'
         stop
      endif
      read(10,*)u%ocean            ; print *,'u%ocean=     ',u%ocean
      read(10,*)rh%ocean           ; print *,'rh%ocean=    ',rh%ocean
      read(10,*)rh%atmos           ; print *,'rh%atmos=    ',rh%atmos
      read(10,*)nsh                ; print *,'nshapiro=    ',nsh

      if (nsh >0) then
         allocate(sh(nsh))
         call shfact(nsh,sh)
      endif

      read(10,'(a)')ca
      if (ca /= '#2') then
         print *,'#2: error in infile.in'
         stop
      endif
      read(10,'(1x,l1)')samp_fix   ; print *,'samp_fix=    ',samp_fix
      read(10,*)inivar%ocean       ; print *,'inistd%ocean=',inivar%ocean ; inivar%ocean=inivar%ocean**2
      read(10,*)inivar%atmos       ; print *,'inistd%atmos=',inivar%atmos ; inivar%atmos=inivar%atmos**2
      read(10,*)sysvar%ocean       ; print *,'sysstd%ocean=',sysvar%ocean ; sysvar%ocean=sysvar%ocean**2
      read(10,*)sysvar%atmos       ; print *,'sysstd%atmos=',sysvar%atmos ; sysvar%atmos=sysvar%atmos**2
      read(10,*)obsvar%ocean       ; print *,'obsstd%ocean=',obsvar%ocean ; obsvar%ocean=obsvar%ocean**2
      read(10,*)obsvar%atmos       ; print *,'obsstd%atmos=',obsvar%atmos ; obsvar%atmos=obsvar%atmos**2

      read(10,'(a)')ca
      if (ca /= '#3') then
         print *,'#3: error in infile.in'
         stop
      endif

      read(10,*)nro                ; print *,'nro=  ',nro
      read(10,*)nra                ; print *,'nra=  ',nra
      read(10,*)obst0o             ; print *,'obst0o=      ',obst0o
      read(10,*)obst0a             ; print *,'obst0a=      ',obst0a
      read(10,*)obsdto             ; print *,'obsdto=      ',obsdto
      read(10,*)obsdta             ; print *,'obsdta=      ',obsdta

      read(10,'(a)')ca
      if (ca /= '#4') then
         print *,'#4: error in infile.in'
         stop
      endif
      read(10,*)mode_analysis      ; print *,'mode_ana=    ',mode_analysis
      read(10,*)truncation        ; print *,'truncation=  ',truncation
      read(10,'(1x,a)')covmodel   ; print *,'covmodel=    ',trim(covmodel)
      read(10,*)rd                ; print *,'rd      =    ',rd
      read(10,'(1x,l1)')Rexact    ; print *,'Rexact=      ',Rexact
      read(10,'(1x,l1)')lrandrot  ; print *,'lrandrot=    ',lrandrot
      read(10,*)inflate,infmult   ; print *,'inflation=   ',inflate,infmult
      read(10,*)local,obs_radius,obs_truncation; print *,'localization=',local,obs_radius,obs_truncation
   close(10)


! We assume constant atmospheric velocity equal to one
   u%atmos=1.0
end subroutine
end module
