#include <cstdio>

#include <cvode/cvode.h>
#include <nvector/nvector_serial.h>
#include <sundials/sundials_types.h>
#include <sunlinsol/sunlinsol_dense.h>
#include <sunmatrix/sunmatrix_dense.h>

static int f(realtype t, N_Vector y, N_Vector ydot, void* user_data);
static int check_flag(void* flagvalue, const char* funcname, int opt);
static void PrintFinalStats(void* cvode_mem);

extern "C" {
int cvodeint_(int* Neqp, double* x1, double* x2, double* y, double* eps);
void rhs_balance_(double*, const double*, const double*);
}

/******************************************************************************/

int cvodeint_(int* Neqp, double*, double* x2, double* y, double*) {
    int const neq = *Neqp;

    N_Vector yv = N_VMake_Serial(neq, y);
    if (!yv) {
        fprintf(stderr, "\nerror: cvodeint: y vector allocation failed!..");
        return 1;
    }

    realtype const tfinal = *x2;
    realtype t;
    int flag;

    // INTEGRATION METHOD -----------------------------------------------------
    // old method
    // void *cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);
    void* cvode_mem = CVodeCreate(CV_ADAMS);
    if (check_flag(cvode_mem, "CVodeCreate", 0))
        return 1;

    SUNMatrix A;
    SUNLinearSolver LS;

    // Create dense SUNMatrix for use in linear solver
    A = SUNDenseMatrix(neq, neq);
    if (check_flag((void*)A, "SUNDenseMatri", 0))
        return 1;

    // Create dense linear solver for use by CVODE
    LS = SUNLinSol_Dense(yv, A);
    if (check_flag((void*)LS, "SUNLinSol_Dense", 0))
        return 1;

    flag = CVodeSetLinearSolver(cvode_mem, LS, A);
    if (check_flag((void*)&flag, "CVodeSetLinearSolver", 1))
        return 1;

    flag = CVode(cvode_mem, tfinal, yv, &t, CV_NORMAL);

    if (flag != CV_SUCCESS) {
        fprintf(stderr, "\nerror: cvode failed!");
    }

    PrintFinalStats(cvode_mem);

    CVodeFree(&cvode_mem);

    return 0;
}

/******************************************************************************/

static int f(realtype t, N_Vector y, N_Vector ydot, void* user_data) {
    rhs_balance_(&t, N_VGetArrayPointer(y), N_VGetArrayPointer(ydot));
    return 0;
}

/******************************************************************************/

static int check_flag(void* flagvalue, const char* funcname, int opt) {
    /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
    if (0 == opt && nullptr == flagvalue) {
        fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n", funcname);
        return 1;
    }

    /* Check if flag < 0 */
    else if (1 == opt) {
        int const errflag = *reinterpret_cast<int*>(flagvalue);
        if (errflag < 0) {
            fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed with flag = %d\n\n", funcname, errflag);
            return 1;
        }
    }

    /* Check if function returned NULL pointer - no memory allocated */
    else if (2 == opt && nullptr == flagvalue) {
        fprintf(stderr, "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n", funcname);
        return 1;
    }

    return 0;
}

/******************************************************************************/

static void PrintFinalStats(void* cvode_mem) {
    long int nst, nfe, nsetups, nje, nfeLS, nni, ncfn, netf, nge;
    int flag;

    flag = CVodeGetNumSteps(cvode_mem, &nst);
    check_flag(&flag, "CVodeGetNumSteps", 1);
    flag = CVodeGetNumRhsEvals(cvode_mem, &nfe);
    check_flag(&flag, "CVodeGetNumRhsEvals", 1);
    flag = CVodeGetNumLinSolvSetups(cvode_mem, &nsetups);
    check_flag(&flag, "CVodeGetNumLinSolvSetups", 1);
    flag = CVodeGetNumErrTestFails(cvode_mem, &netf);
    check_flag(&flag, "CVodeGetNumErrTestFails", 1);
    flag = CVodeGetNumNonlinSolvIters(cvode_mem, &nni);
    check_flag(&flag, "CVodeGetNumNonlinSolvIters", 1);
    flag = CVodeGetNumNonlinSolvConvFails(cvode_mem, &ncfn);
    check_flag(&flag, "CVodeGetNumNonlinSolvConvFails", 1);

    flag = CVodeGetNumJacEvals(cvode_mem, &nje);
    check_flag(&flag, "CVDlsGetNumJacEvals", 1);
    flag = CVodeGetNumRhsEvals(cvode_mem, &nfeLS);
    check_flag(&flag, "CVDlsGetNumRhsEvals", 1);

    flag = CVodeGetNumGEvals(cvode_mem, &nge);
    check_flag(&flag, "CVodeGetNumGEvals", 1);

    printf("\nFinal Statistics:\n");
    printf("nst = %-6ld nfe  = %-6ld nsetups = %-6ld nfeLS = %-6ld nje = %ld\n", nst, nfe, nsetups,
        nfeLS, nje);
    printf("nni = %-6ld ncfn = %-6ld netf = %-6ld nge = %ld\n \n", nni, ncfn, netf, nge);
}

/******************************************************************************/
