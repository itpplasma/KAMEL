/*! \file
    \brief The implementation of functions declared in solver.h.
*/

#include "solver.h"
#include "rhs_func.h"

#include <cvode/cvode.h>               /* main integrator header file */
#include <nvector/nvector_serial.h>    /* serial N_Vector types, fct. and macros */
#include <sundials/sundials_math.h>    /* contains the macros ABS, SQR, and EXP */
#include <sundials/sundials_types.h>   /* definition of realtype */
#include <sunlinsol/sunlinsol_dense.h> // new dense solver, sundials_dense is deprecated
#include <sunmatrix/sunmatrix_dense.h>

SysRHSFcn rhs_mat;

/*-----------------------------------------------------------------*/
static int check_flag(void* flagvalue, const char* funcname, int opt);
/*-----------------------------------------------------------------*/

/* Functions Called by the Solver */
int func(realtype r, N_Vector y, N_Vector ydot, void* Dmat) {
    rhs_mat((double)r, N_VGetArrayPointer(y), N_VGetArrayPointer(ydot), Dmat);
    return 0;
}

/*-----------------------------------------------------------------*/

int integrate_basis_vecs(SysRHSFcn f, int Nfs, int Nw, int dim, double* rvec, double* Smat,
    solver_settings* ss, void* params) {
    /*
f: a rhs function;
Nfs: a number of fundamental solutions to be integrated simulaniously;
Nw: dimension of the problem (number of waves)
dim: dimension of r-grid;
rvec: a grid points where solution has to be found
Smat: a solution; on entrance Smat[0..Neq-1] must be filled by starting values (re,im,re,im,...), on
exit - contains basis vectors at rvec points packed one after another.
*/

    SUNMatrix A;
    SUNLinearSolver LS;

    rhs_mat = f;

    int Neq = 2 * Nfs * Nw;

    int Nort = ss->Nort;

    // a memory block used for internal calulations, of length Nort*(1+Neq+2*Nfs):
    double* mem = new double[Nort * (1 + Neq + 2 * Nfs)];

    // void *cvode_mem = CVodeCreate (CV_BDF, CV_NEWTON);
    //  Call CvodeCreate to create CVode memory block and specify the
    //  ADAMS differentiation formula (for nonstiff problems)
    void* cvode_mem = CVodeCreate(CV_ADAMS);
    if (check_flag((void*)cvode_mem, "CVodeCreate", 0))
        return 1;

    int i, step = 0;

    /* for cvode: keeps the values obtained by integration and ortonorm. */
    double *rdata, *ydata, *taudata;

    rdata = mem;
    ydata = rdata + Nort;
    taudata = ydata + Neq * Nort;

    /* initial values: */
    *rdata = rvec[0];
    for (i = 0; i < Neq; i++)
        ydata[i] = Smat[i];

    N_Vector y = N_VMake_Serial(Neq, ydata);

    if (!y) {
        fprintf(stderr, "\nerror: int_basis_vecs: y vector allocation failed!..");
        return 1;
    }

    N_Vector yval = N_VMake_Serial(Neq, ydata);
    if (!yval) {
        fprintf(stderr, "\nerror: int_basis_vecs: yval vector allocation failed!..");
        return 1;
    }

    realtype reltol = ss->eps_rel, abstol = ss->eps_abs;

    int flag;

    flag = CVodeInit(cvode_mem, func, (realtype)rvec[0], y);
    if (flag != CV_SUCCESS) {
        fprintf(stderr, "\nerror: int_basis_vecs: CVodeInit failed!..");
        return 1;
    }

    flag = CVodeSetMaxOrd(cvode_mem, 12);
    if (flag != CV_SUCCESS) {
        fprintf(stderr, "\nerror: int_basis_vecs: CVodeSetMaxOrd failed!..");
        return 1;
    }

    flag = CVodeSStolerances(cvode_mem, reltol, abstol);
    if (flag != CV_SUCCESS) {
        fprintf(stderr, "\nerror: int_basis_vecs: CVodeSStolerances failed!..");
        return 1;
    }

    // Create dense SUNMatrix for use in linear solver
    A = SUNDenseMatrix(Neq, Neq);
    if (check_flag((void*)A, "SUNDenseMatrix", 0))
        return 1;

    // Create dense linear solver for use by CVode
    LS = SUNLinSol_Dense(y, A);
    if (check_flag((void*)LS, "SUNLinSol_Dense", 0))
        return 1;

    // Attach the linear solver and matrix to CVode
    flag = CVodeSetLinearSolver(cvode_mem, LS, A);
    if (check_flag((void*)&flag, "CVodeSetLinearSolver", 1))
        return 1;

    // Call CVDense to specify the CVDENSE dense linear solver:
    //  flag = CVodeSetLinearSolver(cvode_mem, LS, NULL);
    // flag = CVDense(cvode_mem, Neq);
    // if (flag != CV_SUCCESS)
    //{
    //     fprintf(stderr, "\nerror: int_basis_vecs: cvdense failed!..");
    //     return 1;
    // }

    realtype rf;
    rf = (realtype)rvec[dim - 1];

    flag = CVodeSetStopTime(cvode_mem, rf);
    if (flag != CV_SUCCESS) {
        fprintf(stderr, "\nerror: int_basis_vecs: cvodestoptime failed!..");
        return 1;
    }

    // flag = CVodeSetFdata (cvode_mem, params);
    flag = CVodeSetUserData(cvode_mem, params);
    if (flag != CV_SUCCESS) {
        fprintf(stderr, "\nerror: int_basis_vecs: CVodeSetUserData failed!..\n");
        return (1);
    }

// Jacobian settings:
#if false // USE_JACOBIAN_IN_ODE_SOLVER == 1

flag = CVDlsSetDenseJacFn(cvode_mem, Jacobian);
if (flag != CV_SUCCESS)
{
    fprintf(stderr, "\nerror: int_basis_vecs: CVDenseSetJacFn failed!..\n");
    return(1);
}

#endif

    // We everywhere represent fortran complex arrays as double ones of double size

    int ort_flag, INFO, LWORK = -1;

    double* WORK = (double*)xmalloc(2 * Nfs * sizeof(double));

    /* call to determine an optimal LWORK: ydata should'nt change! */
    zgeqrf_(&Nw, &Nfs, ydata, &Nw, taudata, WORK, &LWORK, &INFO);

    if (INFO) {
        fprintf(stderr, "\nerror: int_basis_vecs: zgeqrf_ failed!: %d", INFO);
        return 1;
    }

    LWORK = (int)WORK[0];

    // fprintf(stdout, "\nint_basis_vecs_: optimal LWORK=%d\n", LWORK);

    WORK = (double*)xrealloc(WORK, 2 * LWORK * sizeof(double));

    double max, min, mod;

    int ind, rpind, flag2, tot_steps;

    double* adr;

    rpind = 0;     /*index of a reached point in rvec */
    tot_steps = 0; /* total number of CVode steps */

    int dirint = signum(rvec[dim - 1] - rvec[0]);

    while (1) {
        if (step == Nort - 1) {
            fprintf(stderr,
                "\nerror: int_basis_vecs: maximum number of ortonormalization steps is reached: %d",
                step);
            fflush(stderr);
            break;
        }

        /*makes one time step, data in y points to a Smat+step*Neq:*/
        flag = CVode(cvode_mem, rf, y, rdata, CV_ONE_STEP);
        tot_steps++;

        /*check for problems:*/
        if (!(flag == CV_SUCCESS || flag == CV_TSTOP_RETURN)) {
            fprintf(stderr, "\nerror: int_basis_vecs: cvode failed!: t=%g\tflag=%d", *rdata, flag);
            fflush(stderr);
            break;
        }

        /* checks for grid points: */
        for (i = rpind + 1; i < dim; i++) {
            if (dirint * (rvec[i] - (*rdata)) <= 0.0) /*determine y values:*/
            {
                NV_DATA_S(yval) = Smat + i * Neq;
                flag2 = CVodeGetDky(cvode_mem, rvec[i], 0, yval);

                if (flag2 != CV_SUCCESS) {
                    fprintf(stderr, "\nerror: int_basis_vecs_: cvodegetdky failed!: t=%g\tflag=%d",
                        *rdata, flag2);
                    fflush(stderr);
                    return 1;
                }
                rpind = i;
            } else
                break;
        }

        /* check for finish: */
        if (flag == CV_TSTOP_RETURN) {
            if (rpind != dim - 1) {
                fprintf(stderr,
                    "\nerror: int_basis_vecs: a wrong index is detected: rpind=%d\tdim=%d", rpind,
                    dim);
            }
            break;
        }

        /*makes QR to check ortoganality; if still Ok, proceed:*/
        /*better to estimate a scale of damping and do not check each step!!!*/

        adr = ydata + Neq; /* a memory pointer for next iteration */
        for (i = 0; i < Neq; i++)
            adr[i] = ydata[i];

        zgeqrf_(&Nw, &Nfs, adr, &Nw, taudata, WORK, &LWORK, &INFO);

        if (INFO) {
            fprintf(stderr, "\nerror: int_basis_vecs: zgeqrf_ failed!: %d", INFO);
            fflush(stderr);
            break;
        }

        /*check for orthogonality:*/
        max = sqrt(adr[0] * adr[0] + adr[1] * adr[1]);
        min = max;

        for (i = 1; i < Nfs; i++) {
            ind = 2 * i * (Nw + 1); /*re(i,j): 2*Nwaves*j+2*i */
            mod = sqrt(adr[ind] * adr[ind] + adr[ind + 1] * adr[ind + 1]);
            if (mod > max)
                max = mod;
            if (mod < min)
                min = mod;
        }

        // if (min < (ss->norm_fac)*max) ort_flag = 0; else ort_flag = 1;
        if (max > (ss->norm_fac) * min)
            ort_flag = 0;
        else
            ort_flag = 1;

        if (ort_flag == 1)
            continue; /*ydata and rdata are inchanged, just next step:*/

        if (ss->debug > 1) {
            fprintf(stdout, "\nint_basis_vecs: QRstep = %5d  r = %9.6f\tmin = %8.3e\tmax = %8.3e",
                step, *rdata, min, max);
        }

        /*if orthoganalization is neccessary: store the solution:*/
        for (i = 0; i < Neq; i++)
            ydata[i] = adr[i]; /*stores the ort. data*/

        /*taudata is alredy in Tvec array*/

        /*preparing to go further: next r-grid pointers:*/
        rdata += 1;
        *rdata = *(rdata - 1);

        /*ydata should be a new orthoganal set of vecs:*/
        ydata += Neq;
        zungqr_(&Nw, &Nfs, &Nfs, ydata, &Nw, taudata, WORK, &LWORK, &INFO);
        if (INFO) {
            fprintf(stderr, "\nerror: int_basis_vecs: zungqr_ failed!: %d", INFO);
            break;
        }

        NV_DATA_S(y) = ydata;

        taudata += 2 * Nfs;
        step++; /* increase of an index of the orthonogalization */

        /*restart solver with new start values*/

        flag = CVodeReInit(cvode_mem, (realtype)(*rdata), y);
        if (flag != CV_SUCCESS) {
            fprintf(stderr, "\nerror: int_basis_vecs: cvodereinit failed!: flag=%d\n", flag);
            break;
        }
    }

    if (ss->debug > 0) {
        fprintf(stdout, "\ninformation: int_basis_vecs: Nsteps = %d\tNorts = %d r = %f", tot_steps,
            step, *rdata);
    }

    /*renormalization of basis vectors: we pass the data for last orth. step*/
    if ((rdata - 1 != mem + step - 1) || (ydata - Neq != mem + (Nort) + Neq * (step - 1)) ||
        (taudata - 2 * Nfs != mem + (Nort) + Neq * (Nort) + 2 * Nfs * (step - 1))) {
        fprintf(stderr, "\nerror: int_basis_vecs: wrong pointers to renormalization info:");
        fprintf(stderr, "\nrdata1=%p rdata2=%p", rdata - 1, mem + step - 1);
        fprintf(stderr, "\nydata1=%p ydata2=%p", ydata - Neq, mem + (Nort) + Neq * (step - 1));
        fprintf(stderr, "\ntdata1=%p tdata2=%p", taudata - 2 * Nfs,
            mem + (Nort) + Neq * (Nort) + 2 * Nfs * (step - 1));
        return 1;
    }

    // fprintf(stdout, "\ncheck: integrate_basis_vecs_: rdata=%d ydata=%d udata=%d", (int)mem,
    // (int)(mem+(Nort)), (int)Smat);

    Nort = step;

    /* pass the pointers to the solution at the last point: */
    renorm_basis_vecs_(Nfs, Nw, dim, rvec, Smat + (dim - 1) * Neq, Nort, rdata - 1, ydata - Neq,
        taudata - 2 * Nfs);

    N_VDestroy_Serial(y);
    N_VDestroy_Serial(yval);
    CVodeFree(&cvode_mem);

    free(WORK);

    delete[] mem;

    return (0);
}

/*-----------------------------------------------------------------*/

int renorm_basis_vecs_(int Nfs, int Nw, int dim, double* rvec, double* udata, int Nort,
    double* rdata, double* ydata, double* taudata) {
    /*
all complex arrays are stored as double 1D ones of double length:
*/

    int Neq = 2 * Nfs * Nw;

    double* tmp = (double*)xmalloc(Neq * sizeof(double));
    double* buf = (double*)xmalloc(Neq * sizeof(double));

    int i, k, ind;

    double alpha[2] = {1.0, 0.0}, beta[2] = {0.0, 0.0};

    char side = 'L', trans = 'N', uplo = 'U', diag = 'N';

    int info;

    /* buf is the nearly unit matrix after last orth.: */
    for (k = 0; k < Nw; k++) {
        for (i = 0; i < Nfs; i++) {
            ind = 2 * Nw * i + 2 * k; /*index of a real part of (k,i) element*/
            buf[ind] = 0.0;
            buf[ind + 1] = 0.0;
            if (i == k)
                buf[ind] = 1.0;
        }
    }

    /* main loop over r-points: in direction back to integration */
    /* for last points buf = unit matrix */

    int dirint = signum(rvec[dim - 1] - rvec[0]);
    int step, ropind = Nort;

    for (k = dim - 1; k > -1; k--) /* index of r grid point */
    {
        for (step = ropind - 1; step > -1; step--) {
            if (dirint * (rvec[k] - (*rdata)) > 0.0)
                break;

            ropind--;

            /*a new transformation needed: buf=tmp*buf*/
            /*
                fprintf(stdout, "\nrenorm_basis_vecs_: new matrix: ropind=%4d r[%4d]=%6.3f
               *rort=%6.3f ydata=%d\n", ropind, k, rvec[k], *rdata, (int)ydata);
                */
            /* finds inverse of R(k) and owerwrites itself in ydata */
            ztrtri_(&uplo, &diag, &Nfs, ydata, &Nw, &info);
            if (info) {
                fprintf(stderr,
                    "\nerror: renorm_basis_vecs_: ztrtri_ failed!: info=%d k=%d r=%f rdata=%p "
                    "ydata=%p",
                    info, k, rvec[k], rdata, ydata);
                /*return 1;*/
            }

            /*buf is a inv(R(k+1))*inv(R(k+2))...*inv(R(Nort-1))*/
            /* multiplies inv(R(k)) on buf to get a new buf */
            ztrmm_(&side, &uplo, &trans, &diag, &Nfs, &Nfs, alpha, ydata, &Nw, buf, &Nw);
            if (info) {
                fprintf(stderr, "\nerror: renorm_basis_vecs_: ztrmm_ failed!: %d", info);
                /*return 1;*/
            }
            rdata--;
            ydata -= Neq;
            taudata -= 2 * Nfs;
        }

        /*multiplies U(k) on tmp: this is final numbers to store in udata*/
        zgemm_(&trans, &trans, &Nw, &Nfs, &Nfs, alpha, udata, &Nw, buf, &Nw, beta, tmp, &Nw);

        /* fills Smat for point k by values stored in tmp: */
        for (i = 0; i < Neq; i++)
            udata[i] = tmp[i];
        /*
        fprintf(stdout, "\rrenorm_basis_vecs_: ropind=%4d r[%4d]=%6.3f *rort=%6.3f ydata=%d",
        ropind, k, rvec[k], *rdata, (int)ydata);
        */
        udata -= Neq; /* next point */
    }

    // fprintf(stdout, "\ncheck: renorm_basis_vecs_: rdata=%d ydata=%d udata=%d", (int)(rdata+1),
    // (int)(ydata+Neq), (int)(udata+Neq));

    free(buf);
    free(tmp);

    return 0;
}

/*-----------------------------------------------------------------*/

int superpose_basis_vecs_(int Nfs, int Nw, int dim, double* Smat, double* Cvec, double* Svec) {
    /*
all complex arrays are stored as a double 1D ones of double length:
dim Smat = 2*Nw*Nfs*nsteps
dim Cvec = 2*Nfs
dim Svec = 2*Nw*nsteps
*/

    int Neq = 2 * Nfs * Nw;

    double *udata, *sdata;
    udata = Smat; /* initial values for reference */
    sdata = Svec; /* initial values for reference */

    double alpha[2] = {1.0, 0.0}, beta[2] = {0.0, 0.0};
    char trans = 'N';

    int k, incx = 1, incy = 1;

    /* main loop over r-points: */
    for (k = 0; k < dim; k++) {
        /* multiply U on C to superpose */
        zgemv_(&trans, &Nw, &Nfs, alpha, udata, &Nw, Cvec, &incx, beta, sdata, &incy);
        udata += Neq;
        sdata += 2 * Nw;
    }

    return 0;
}

/*-----------------------------------------------------------------*/
static int check_flag(void* flagvalue, const char* funcname, int opt) {

    int* errflag;

    // Check if sundials function returned NULL pointer - no memory allocated
    if (opt == 0 && flagvalue == NULL) {
        fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n", funcname);
        return 1;
    }

    // Check if flag < 0
    else if (opt == 1) {
        errflag = (int*)flagvalue;
        if (*errflag < 0) {
            fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed with flag = %d\n\n", funcname, *errflag);
            return 1;
        }
    }

    // Check if fnc returned NULL pointer - no memory allocated
    else if (opt == 2 && flagvalue == NULL) {
        fprintf(stderr, "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n", funcname);
        return 1;
    }

    return 0;
}
