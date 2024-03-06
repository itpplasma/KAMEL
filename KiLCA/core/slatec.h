/*! \file slatec.h
    \brief Some declarations of SLATEC Fortran subroutines.
*/

#if defined(__cplusplus)
extern "C" {
#endif

void dpsort_ (double *, int *, int *, int *, int *);

int ddeabm_ (void (*df)(double *, double *, double *, double *, int *),
                const int *neq, double *t, double *y, double *tout, int *info,
                double *rtol, double *atol, int *idid, double *rwork,
                const int *lrw, int *iwork, const int *liw, double *rpar,
                int *ipar);

int ddassl_ (void (*res)(double *r, double *u, double *du, double *res, int *ires, double *rpar, int *ipar),
int *, double *, double *, double *, double *, int *, double *, double *,
int *, double *, const int *, int *, const int *, double *, int *, int *);

int dnsqe_ (void (*nlsys)(int *Neq, double *u, double *res, int *flag), int *,  int *,  const int *, double *, double *, double *,  int *,  int *, double *,  const int *);


#if defined(__cplusplus)
}
#endif
