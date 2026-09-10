!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine get_background_dimension_from_balance (dim_p)

implicit none;

integer, intent(out) :: dim_p;

print *, 'warning: dummy subroutine get_background_dimension_from_balance_code () is called!'

dim_p = 0;

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine get_background_profiles_from_balance (dim_p, r_p, q_p, n_p, Ti_p, Te_p, Vth_p, Vz_p, Er_p)

implicit none;

integer, intent(out) :: dim_p;
real(8), dimension(*), intent(out) :: r_p;
real(8), dimension(*), intent(out) :: q_p;
real(8), dimension(*), intent(out) :: n_p;
real(8), dimension(*), intent(out) :: Ti_p;
real(8), dimension(*), intent(out) :: Te_p;
real(8), dimension(*), intent(out) :: Vth_p;
real(8), dimension(*), intent(out) :: Vz_p;
real(8), dimension(*), intent(out) :: Er_p;

print *, 'warning: dummy subroutine get_background_profiles_from_balance_code () is called!'

dim_p = 0;

! The dummy supplies zero points, so there are no profile elements to fill.

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
