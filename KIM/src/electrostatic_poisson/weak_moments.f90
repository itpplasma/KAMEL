module weak_moments_m
    use KIM_kinds_m, only: dp
    implicit none
    private
    public :: weak_moment_projector_t

    ! Physical moments use the raw dr mass, with no field boundary conditions.
    type :: weak_moment_projector_t
        private
        real(dp), allocatable :: diagonal(:), offdiagonal(:)
    contains
        procedure :: init => projector_init
        procedure :: project => projector_project
    end type weak_moment_projector_t

contains

    subroutine projector_init(self, r)
        class(weak_moment_projector_t), intent(inout) :: self
        real(dp), intent(in) :: r(:)
        real(dp) :: h
        integer :: i, n, info
        external :: dpttrf

        n = size(r)
        if (n < 2) error stop 'Weak moment grid requires at least two points'
        if (any(r(2:) <= r(:n-1))) error stop 'Weak moment grid must increase'
        if (allocated(self%diagonal)) deallocate(self%diagonal, self%offdiagonal)
        allocate(self%diagonal(n), self%offdiagonal(n-1))
        self%diagonal = 0.0_dp
        do i = 1, n-1
            h = r(i+1)-r(i)
            self%diagonal(i:i+1) = self%diagonal(i:i+1)+h/3.0_dp
            self%offdiagonal(i) = h/6.0_dp
        end do
        call dpttrf(n, self%diagonal, self%offdiagonal, info)
        if (info /= 0) error stop 'Raw moment mass factorization failed'
    end subroutine projector_init

    subroutine projector_project(self, weak_load, moment)
        class(weak_moment_projector_t), intent(in) :: self
        complex(dp), intent(in) :: weak_load(:)
        complex(dp), allocatable, intent(out) :: moment(:)
        real(dp), allocatable :: rhs(:, :)
        integer :: n, info
        external :: dpttrs

        if (.not. allocated(self%diagonal)) error stop 'Moment projector not initialized'
        n = size(self%diagonal)
        if (size(weak_load) /= n) error stop 'Weak load and moment grid differ'
        allocate(rhs(n, 2))
        rhs(:, 1) = real(weak_load, dp)
        rhs(:, 2) = aimag(weak_load)
        call dpttrs(n, 2, self%diagonal, self%offdiagonal, rhs, n, info)
        if (info /= 0) error stop 'Raw moment mass solve failed'
        moment = cmplx(rhs(:, 1), rhs(:, 2), dp)
    end subroutine projector_project

end module weak_moments_m
