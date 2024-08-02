module diag_mod

    use control_mod

    implicit none
    logical :: write_diag,write_diag_b
    integer :: iunit_diag,iunit_diag_b
    integer :: i_mn_loop
    logical :: toggle
    double precision :: relax_ff
    double precision, dimension(:), allocatable :: field_fac_old

    integer, dimension(1) :: ind_dqle, ind_dqli

    double precision :: timscal_dql, rate_dql, timscal_dqli

    contains

    subroutine determineDqlDiagnostic

        use time_evolution, only: dqle11_prev, dqli11_prev, timstep
        use grid_mod, only: dqle11, dqli11

        implicit none

        timscal_dql = maxval(abs(dqle11_prev - dqle11))/maxval(dqle11_prev + dqle11)
        ind_dqle = maxloc(abs(dqle11_prev - dqle11))
        timscal_dqli = maxval(abs(dqli11_prev - dqli11))/maxval(dqli11_prev + dqli11)
        ind_dqli = maxloc(abs(dqli11_prev - dqli11))
        rate_dql = timscal_dql/timstep

    end subroutine

end module diag_mod
