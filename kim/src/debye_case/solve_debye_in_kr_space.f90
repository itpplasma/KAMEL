subroutine solve_debye_in_kr_space

    use debye_kernel
    use grid, only: kr_grid, rg_grid
    use constants, only: e_charge, com_unit, pi
    use sparse_mod, only: sparse_solveComplex_b1, column_pointer2full,&
     sparse_solveComplex_A_b1, sparse_solve_method, sparse_solve
    use config, only: output_path

    implicit none

    integer :: i,j
    double complex :: A_matrix(kr_grid%npts_b, kr_grid%npts_b)
    double complex :: b_vec(kr_grid%npts_b)
    double complex :: phi_in_rg(rg_grid%npts_b)

    integer :: sparse_solver_option

    A_matrix = 0.0d0
    b_vec = 0.0d0

    call calculate_debye_length

    open(unit = 80, file = trim(output_path)//'fields/A_mat.dat')
    open(unit = 81, file = trim(output_path)//'fields/b_vec.dat')
    do i = 1, kr_grid%npts_b
        A_matrix(i, i) = kr_grid%xb(i)**2.0d0 + 1.0d0 / lambda_D**2
        write(80,*) kr_grid%xb(i), real(A_matrix(i, i)), dimag(A_matrix(i, i))
        b_vec(i) = 2.0d0 * e_charge * exp(-com_unit * rg_grid%max_val/2.0d0 * kr_grid%xb(i))
        write(81,*) kr_grid%xb(i), real(b_vec(i)), dimag(b_vec(i))
    end do
    close(80)
    close(81)

    !call dense_to_sparse(A_matrix, irow, pcol, A_nz, nrow, ncol, nz_out)
    sparse_solver_option = 3
    sparse_solve_method = 3
    call sparse_solveComplex_A_b1(A_matrix, b_vec, sparse_solver_option)
    !call sparse_solve(A_matrix, b_vec)
    !call sparse_solveComplex_b1(nrow, ncol, nz_out, irow, pcol, A_nz, b_vec, sparse_solver_option)

    open(unit = 80, file = trim(output_path)//'fields/phi_k.dat')
    do i =1, kr_grid%npts_b
        write(80,*) kr_grid%xb(i), real(b_vec(i)), dimag(b_vec(i))
    end do
    close(80)

    phi_in_rg = cmplx(0.0d0, 0.0d0)
    do i = 1, kr_grid%npts
        do j =1, rg_grid%npts_b
            phi_in_rg(j) = phi_in_rg(j) + b_vec(i) * exp(com_unit * rg_grid%xb(j) * kr_grid%xb(i))
        end do
    end do


    open(unit = 80, file = trim(output_path)//'fields/phi_in_rg.dat')
    do i =1, rg_grid%npts_b
        write(80,*) rg_grid%xb(i), real(phi_in_rg(i)), dimag(phi_in_rg(i))
    end do
    close(80)


end subroutine

