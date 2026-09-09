subroutine kim_prepare_resonances

    use kim_resonances_m, only: iunit_res, r_res
    use config_m, only: type_of_run
    use setup_m, only: m_mode, n_mode, type_br_field
    use species_m, only: plasma
    use KIM_kinds_m, only: dp
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite

    implicit none

    integer :: j
    real(dp) :: qres, fraction
    logical :: found

    iunit_res = 157
    r_res = 0.0_dp
    if (n_mode == 0) then
        write (*, *) 'No rational resonance for m=', m_mode, ', n=', n_mode
        return
    end if
    if (.not. allocated(plasma%q) .or. .not. allocated(plasma%r_grid)) &
        error stop 'resonance location requires q and radial profiles'
    if (plasma%grid_size < 2 .or. size(plasma%q) /= plasma%grid_size .or. &
        size(plasma%r_grid) /= plasma%grid_size) &
        error stop 'resonance profile grid mismatch'
    if (.not. all(ieee_is_finite(plasma%q)) .or. &
        .not. all(ieee_is_finite(plasma%r_grid))) &
        error stop 'resonance profiles must be finite'
    if (any(plasma%r_grid(2:) <= plasma%r_grid(:plasma%grid_size - 1))) &
        error stop 'resonance radial grid must be strictly increasing'

    ! The cylindrical Fourier convention gives k_parallel=0 at q=-m/n.
    ! Search in radial order so reversed shear and plateaus select the
    ! innermost crossing deterministically, without a dimensionless fallback.
    qres = -real(m_mode, dp) / real(n_mode, dp)
    found = .false.
    do j = 1, plasma%grid_size
        if (plasma%q(j) == qres) then
            r_res = plasma%r_grid(j)
            found = .true.
            exit
        end if
        if (j == plasma%grid_size) cycle
        if (qres > min(plasma%q(j), plasma%q(j + 1)) .and. &
            qres < max(plasma%q(j), plasma%q(j + 1))) then
            fraction = (qres - plasma%q(j)) / (plasma%q(j + 1) - plasma%q(j))
            r_res = (1.0_dp - fraction) * plasma%r_grid(j) + fraction * plasma%r_grid(j + 1)
            found = .true.
            exit
        end if
    end do
    if (.not. found) then
        write (*, *) 'No resonance for m=', m_mode, ', n=', n_mode, ', signed q=-m/n=', qres
        return
    end if

    ! Keep the artificial point-charge placement for legacy nonperiodic
    ! solvers. A periodic window must always be centered on the physical root.
    if (type_br_field == 2 .and. trim(type_of_run) /= 'electrostatic_periodic') &
        r_res = plasma%r_grid(plasma%grid_size) / 2.0_dp

    write (*, *) 'resonant radius: ', r_res

end subroutine kim_prepare_resonances
