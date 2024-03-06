/*! \file
    \brief C declarations of some of LAPACK Fortran subroutines.
*/

#if defined(__cplusplus)
extern "C" {
#endif

/*-----------------------------------------------------------------*/

void zgeqrf_ (int *, int *, double *, int *, double *, double *, int *, int *);

void zungqr_ (int *, int *, int *, double *, int *, double *, double *, int *, int *);

void zunmqr_ (char *, char *, int *, int *, int *, double *, int *, double *, int *, double *, int *, int *);

void zgemm_ (char *, char *, int *, int *, int *, double *, double *, int *, double *, int *, double *, double *, int *);

void ztrtri_ (char *, char *, int *, double *, int *, int *);

void zgemv_ (char *, int *, int *, double *, double *, int *, double *, int *, double *, double *, int *);

void ztrmm_ (char *, char *, char *, char *, int *, int *, double *, double *, int *, double *, int *);

void zgeevr_ (char *, char *, int *, double *, int *, double *, double *, int *, double *, int *, double *, int *, double *, int *);

void zgetrf_ (int *, int *, double *, int *, int *, int *);

void zgetri_ (int *, double *, int *, int *, double *, int *, int *);

void zgesv_ (int *, int *, double *, int *, int *, double *, int *, int *);

void dgesv_ (int *, int *, double *, int *, int *, double *, int *, int *);

/*-----------------------------------------------------------------*/

#if defined(__cplusplus)
}
#endif
