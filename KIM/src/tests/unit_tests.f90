module unit_tests

    implicit none
    integer :: ierr = 0

    contains

    subroutine test_all(ierr)

        implicit none
        integer, intent(inout) :: ierr

        call test_debye_kernel(ierr)

    end subroutine

    subroutine test_debye_kernel(ierr)

        use debye_kernel
        implicit none
        double complex :: result
        integer, intent(inout) :: ierr

        write(*,*) '(T): Testing debye kernel'
        
        result = func_debye_kernel(0.0d0, 1.0d0)
        if (result /= 0.0d0) then
            write(*,*) '(T): Test failed'
            ierr = ierr + 1
        else
            write(*,*) '(T): Test passed'
        end if

        result = func_debye_kernel(1.0d0, 1.0d0)

        if (result /= 0.0d0) then
            write(*,*) '(T): Test passed'
        else
            write(*,*) '(T): Test failed'
            ierr = ierr + 1
        end if

    end subroutine

    subroutine test_sparse_solver(ierr)

        use sparse_mod, only: sp2fullComplex, sparse_solveComplex_b1, column_pointer2full, sparse_solve_suitesparseComplex_b1, &
                          sparse_solve_method
        use poisson_solver, only: dense_to_sparse
        use config, only: output_path
        implicit none
        integer, intent(in) :: ierr
        integer, dimension(:), allocatable :: irow, pcol, icol
        integer :: nz_out, nrow, ncol
        integer :: iopt

        integer :: dim = 128
        integer :: i,j
        logical :: ex

        double complex, dimension(:), allocatable :: A_nz ! non-zero elements of A matrix
        double complex, dimension(:,:), allocatable :: A_mat ! A matrix
        double complex, dimension(:,:), allocatable :: A_sparse_check ! A matrix reconfigured from sparse matrix

        double complex, dimension(:), allocatable :: b_vec ! b vector and x vector

        allocate(A_mat(dim, dim))
        allocate(b_vec(dim))

        A_mat(i,j) = 0.0d0

        do i= 1, dim
            A_mat(i,i) = -1e8
        end do

        b_vec(dim/2) = 1.0d0

        sparse_solve_method = 1
        call dense_to_sparse(A_mat, irow, pcol, A_nz, nrow, ncol, nz_out)
        call sparse_solveComplex_b1(nrow, ncol, nz_out, irow, pcol, A_nz, b_vec, iopt)

        inquire(file=trim(output_path)//'test_outputs', exist=ex)
        if (.not. ex) then
            call system('mkdir -p '//trim(output_path)//'test_outputs')
        end if
        open(unit = 79, file=trim(output_path)//'test_outputs/x_re.dat')
        open(unit = 80, file=trim(output_path)//'test_outputs/x_im.dat')
        do i = 1, dim
            write(79,*) i, real(b_vec(i))
            write(80,*) i, dimag(b_vec(i))
        end do
        close(79)
        close(80)

    end subroutine

end module