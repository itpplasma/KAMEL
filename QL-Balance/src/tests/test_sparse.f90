program test_sparse

    use sparse_mod

    implicit none

    sparse_talk = .true.
    call sparse_example(1)
    print *, "This should be the answer: "
    print *, "{-1.17021, 0.234043, -0.0212766, -2.28723, 0.255319}"

end program