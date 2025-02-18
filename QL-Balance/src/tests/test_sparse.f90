program test_sparse

    use sparse_mod

    implicit none

    sparse_talk = .true.
    sparse_solve_method = 3

    print *, "Sparse talk: ", sparse_talk

end program