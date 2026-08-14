program test_ntv_feedback
   use QLBalance_kinds, only: dp
   use baseparam_mod, only: am, btor, c, ev, p_mass, rtor, Z_i
   use grid_mod, only: npoib, npoic, nbaleqs, rb, rc, Ercov, ipbeg, ipend, &
                       deriv_coef, reint_coef, sqrt_g_times_B_theta_over_c
   use neort_interface, only: prepare_profile_data_for_neort
   use plasma_parameters, only: params, qsaf
   use rhs_balance_m, only: compute_unsmoothed_radial_electric_field
   use wave_code_data, only: q, Vth

   implicit none

   integer :: i, num_failed
    real(dp), dimension(3) :: E_r_a, E_r_b
    real(dp), dimension(3) :: r, s_tor
    real(dp), dimension(3, 2) :: profile_a, profile_b, profile_repeat
   real(dp) :: expected_mach, vth_ion

   num_failed = 0
   call setup_fixture()

   params(2, :) = 1.0e4_dp
   call compute_unsmoothed_radial_electric_field(params, E_r_a)
   Ercov = -huge(1.0_dp)
   call prepare_profile_data_for_neort(profile_a, r, s_tor)

   vth_ion = sqrt(2.0_dp*params(4, 1)/(am*p_mass))
   do i = 1, size(r)
      expected_mach = c*E_r_a(i)/(r(i)*btor/qsaf(i))*rtor/vth_ion
      call assert_equal(profile_a(i, 2), expected_mach, "Mach from accepted E_r")
   end do

   Ercov = huge(1.0_dp)
   call prepare_profile_data_for_neort(profile_repeat, r, s_tor)
   do i = 1, size(r)
      call assert_equal(profile_repeat(i, 2), profile_a(i, 2), "global Ercov independence")
   end do

   params(2, :) = 2.0e4_dp
   call compute_unsmoothed_radial_electric_field(params, E_r_b)
   call prepare_profile_data_for_neort(profile_b, r, s_tor)
   do i = 1, size(r)
      call assert_equal(E_r_b(i), 2.0_dp*E_r_a(i), "immediate E_r update")
      call assert_equal(profile_b(i, 2), 2.0_dp*profile_a(i, 2), "immediate Mach update")
   end do

   params(2, :) = -1.0e4_dp
   call compute_unsmoothed_radial_electric_field(params, E_r_b)
   call prepare_profile_data_for_neort(profile_b, r, s_tor)
   do i = 1, size(r)
      call assert_equal(profile_b(i, 2), -profile_a(i, 2), "Mach sign")
   end do

   call teardown_fixture()

   if (num_failed > 0) then
      print '(A,I0)', "FAILED assertions: ", num_failed
      stop 1
   end if
   print *, "ALL TESTS PASSED"

contains

   subroutine setup_fixture()
      npoib = 3
      npoic = 3
      nbaleqs = 4

      allocate (rb(npoib), rc(npoic), Ercov(npoib))
      allocate (ipbeg(npoib), ipend(npoib))
      allocate (deriv_coef(1, npoib), reint_coef(1, npoib))
      allocate (sqrt_g_times_B_theta_over_c(npoib))
      allocate (params(nbaleqs, npoic), qsaf(npoic))
      allocate (q(npoib), Vth(npoib))

      rb = [10.0_dp, 20.0_dp, 30.0_dp]
      rc = rb
      r = rb
      s_tor = [0.1_dp, 0.5_dp, 0.9_dp]
      ipbeg = [1, 2, 3]
      ipend = ipbeg
      deriv_coef = 0.0_dp
      reint_coef = 1.0_dp
      sqrt_g_times_B_theta_over_c = 1.0_dp
      q = 1.0_dp
      Vth = 0.0_dp
      qsaf = 1.0_dp

      params(1, :) = 1.0e13_dp
      params(2, :) = 0.0_dp
      params(3, :) = 100.0_dp*ev
      params(4, :) = 100.0_dp*ev

      am = 2.0_dp
      Z_i = 1.0_dp
      btor = 2.0e4_dp
      rtor = 100.0_dp
   end subroutine setup_fixture

   subroutine teardown_fixture()
      deallocate (rb, rc, Ercov)
      deallocate (ipbeg, ipend, deriv_coef, reint_coef)
      deallocate (sqrt_g_times_B_theta_over_c)
      deallocate (params, qsaf, q, Vth)
   end subroutine teardown_fixture

   subroutine assert_equal(actual, expected, label)
      real(dp), intent(in) :: actual, expected
      character(len=*), intent(in) :: label
      real(dp), parameter :: tolerance = 1.0e-12_dp
      real(dp) :: scale

      scale = max(abs(expected), 1.0_dp)
      if (abs(actual - expected) > tolerance*scale) then
         print '(A,A)', "FAIL: ", label
         print '(A,ES22.14)', "  expected: ", expected
         print '(A,ES22.14)', "  actual:   ", actual
         num_failed = num_failed + 1
      end if
   end subroutine assert_equal

end program test_ntv_feedback
