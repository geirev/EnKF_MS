module mod_state
! Model state definition
   use mod_dimensions
   private
   public :: state, operator(+), operator(-), operator(*), operator(/), assignment(=), sqrt, average, abs
   public :: substate

   type substate
      real Ocean
      real Atmos
   end type substate

   type state
! Solution data
      real Ocean(nx)
      real Atmos(nx)
   end type state

! Overloaded and generic operators
   interface operator(+)
      module procedure add_state,&
                       add_substate
   end interface

   interface operator(-)
      module procedure subtract_state
   end interface

   interface operator(*)
      module procedure state_real_mult,&
                       real_state_mult,&
                       state_state_mult,&
                       substate_substate_mult
   end interface

   interface operator(/)
      module procedure state_state_div,&
                       substate_real_div
   end interface

   interface assignment(=)
      module procedure assign_state
      module procedure assign_substate
   end interface

   interface sqrt
      module procedure sqrt_state
      module procedure sqrt_substate
   end interface

   interface abs
      module procedure abs_state
   end interface

   interface average
      module procedure average_state
   end interface

contains

!   function add_constant_to_state(A,r)
!      type(state) add_constant_to_state
!      type(state), intent(in) :: A
!      real, intent(in) :: r
!      add_constant_to_state%Atmos(:)       = A%Atmos(:) + r
!      add_constant_to_state%Ocean(:)       = A%Ocean(:) + r
!   end function add_constant_to_state

   function average_state(A)
      type(substate) average_state
      type(state), intent(in) :: A
      average_state%Atmos       = sum(A%Atmos(:))/real(nx)
      average_state%Ocean       = sum(A%Ocean(:))/real(nx)
   end function average_state

   function abs_state(A)
      type(state) abs_state
      type(state), intent(in) :: A
      abs_state%Atmos       = abs(A%Atmos)
      abs_state%Ocean       = abs(A%Ocean)
   end function abs_state

   function sqrt_state(A)
      type(state) sqrt_state
      type(state), intent(in) :: A
      real :: eps=0.1E-14
      sqrt_state%Atmos       = sqrt(A%Atmos+eps)
      sqrt_state%Ocean       = sqrt(A%Ocean+eps)
   end function sqrt_state

   function sqrt_substate(A)
      type(substate) sqrt_substate
      type(substate), intent(in) :: A
      real :: eps=0.1E-14
      sqrt_substate%Atmos       = sqrt(A%Atmos+eps)
      sqrt_substate%Ocean       = sqrt(A%Ocean+eps)
   end function sqrt_substate

   function add_state(A,B)
      type(state) add_state
      type(state), intent(in) :: A
      type(state), intent(in) :: B
      add_state%Atmos       = A%Atmos  + B%Atmos
      add_state%Ocean       = A%Ocean  + B%Ocean
   end function add_state

   function add_substate(A,B)
      type(substate) add_substate
      type(substate), intent(in) :: A
      type(substate), intent(in) :: B
      add_substate%Atmos       = A%Atmos  + B%Atmos
      add_substate%Ocean       = A%Ocean  + B%Ocean
   end function add_substate

   function subtract_state(A,B)
      type(state) subtract_state
      type(state), intent(in) :: A
      type(state), intent(in) :: B
      subtract_state%Atmos       = A%Atmos  - B%Atmos
      subtract_state%Ocean       = A%Ocean  - B%Ocean
   end function subtract_state

   function state_real_mult(A,B)
      type(state) state_real_mult
      type(state), intent(in) :: A
      real, intent(in) :: B
      state_real_mult%Atmos       = B*A%Atmos
      state_real_mult%Ocean       = B*A%Ocean
   end function state_real_mult

   function real_state_mult(B,A)
      type(state) real_state_mult
      type(state), intent(in) :: A
      real, intent(in) :: B
      real_state_mult%Atmos       = B*A%Atmos
      real_state_mult%Ocean       = B*A%Ocean
   end function real_state_mult

   function state_state_mult(A,B)
      type(state) state_state_mult
      type(state), intent(in) :: A
      type(state), intent(in) :: B
      state_state_mult%Atmos       = A%Atmos  * B%Atmos
      state_state_mult%Ocean       = A%Ocean  * B%Ocean
   end function state_state_mult

   function substate_substate_mult(A,B)
      type(substate) substate_substate_mult
      type(substate), intent(in) :: A
      type(substate), intent(in) :: B
      substate_substate_mult%Atmos       = A%Atmos  * B%Atmos
      substate_substate_mult%Ocean       = A%Ocean  * B%Ocean
   end function substate_substate_mult

   function state_state_div(A,B)
      type(state) state_state_div
      type(state), intent(in) :: A
      type(state), intent(in) :: B
      integer i
      do i=1,nx
         state_state_div%Atmos(i)       = A%Atmos(i) / B%Atmos(i)
         state_state_div%Ocean(i)       = A%Ocean(i) / B%Ocean(i)
      enddo
   end function state_state_div


   function substate_real_div(A,r)
      type(substate) substate_real_div
      type(substate), intent(in) :: A
      real, intent(in) :: r
      substate_real_div%Atmos       = A%Atmos / r
      substate_real_div%Ocean       = A%Ocean / r
   end function substate_real_div

   subroutine assign_state(A,r)
      type(state), intent(out) :: A
      real, intent(in) :: r
      A%Atmos       = r
      A%Ocean       = r
   end subroutine assign_state

   subroutine assign_substate(A,r)
      type(substate), intent(out) :: A
      real, intent(in) :: r
      A%Atmos       = r
      A%Ocean       = r
   end subroutine assign_substate

end module mod_state

