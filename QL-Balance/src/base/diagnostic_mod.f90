module QLbalance_diag

    use control_mod
    use QLBalance_kinds, only: dp

    implicit none

    logical :: write_diag = .false.
    logical :: write_diag_b = .false.
    logical :: toggle

    integer :: iunit_diag,iunit_diag_b
    integer :: i_mn_loop
    integer, dimension(1) :: ind_dqle, ind_dqli

    real(dp) :: timscal_dql, rate_dql, timscal_dqli

end module QLbalance_diag
