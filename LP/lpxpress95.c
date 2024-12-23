/****************************************************************************/
/*                                                                          */
/*            Interface Routines to Xpress LP solver                        */
/*                                                                          */
/*                                                                          */
/*  NOTE: Use this code in place of lp_none.c to access the Xpress          */
/*   libraries. You will also need to link the Xpress library via the       */
/*   makefile.                                                              */
/*                                                                          */
/*  Written by:  Francesco Cavaliere                                        */
/*  Date: Dec 20, 2024 (cava)                                               */
/*                                                                          */
/*  Original Interface written by:  Applegate, Bixby, Chvatal, and Cook     */
/*  Date: April 17, 1997                                                    */
/*        June 19, 1997 (bico, REB)                                         */
/*        January 13, 2002 (bico)                                           */
/*                                                                          */
/*    EXPORTED FUNCTIONS:                                                   */
/*                                                                          */
/*  int CClp_init (CClp **lp)                                               */
/*    INITIALIZES the LP.                                                   */
/*                                                                          */
/*  int CClp_force_perturb (CClp *lp)                                       */
/*    Forces a perturbation in the LP simplex solves                        */
/*                                                                          */
/*  int CClp_tune_small (CClp *lp)                                          */
/*    SETS solver options for tiny problems.                                */
/*                                                                          */
/*  int CClp_disable_presolve (CClp *lp)                                    */
/*    DISABLES the solvers presolve.                                        */
/*                                                                          */
/*  void CClp_free (CClp **lp)                                              */
/*    FREES the LP, both the allocation in CClp_init () and the             */
/*    allocation in CClp_loadlp.                                            */
/*                                                                          */
/*  void CClp_freelp (CClp **lp)                                            */
/*    FREES the LP loaded in CClp_loadlp (or CClp_create), but does         */
/*    not free the data allocated by CClp_init.                             */
/*                                                                          */
/*  int CClp_loadlp (CClp *lp, const char *name, int ncols, int nrows,      */
/*      int objsense, double *obj, double *rhs, char *sense,                */
/*      int *matbeg, int *matcnt,                                           */
/*      int *matind, double *matval,                                        */
/*      double *lb, double *ub)                                             */
/*    LOADS the data into the LP.                                           */
/*      -name attaches a name to the LP (it can be used by the LP solver    */
/*       in io routines)                                                    */
/*      -ncols and nrows give the number of columns and rows in the LP      */
/*      -objsense should be 1 for minimize and -1 for maximize              */
/*      -obj and rhs are arrays giving the objective function and rhs       */
/*      -sense is an array specifying 'L', 'E', or 'G' for each of the      */
/*       rows                                                               */
/*      -matbeg, matcnt, matind, and matval give the coefficients of        */
/*       the constraint matrix in column by column order. matbeg gives      */
/*       gives the index of the start of each column; matcnt gives the      */
/*       number of coefficients in each column; matind gives the indices    */
/*       of the rows where the coefficients are located in the constraint   */
/*       matrix (so for column j, the indices are given in matcnt[j]        */
/*       locations starting at matind[matbeg[j]]; and matval gives the      */
/*       actual coefficients (organized like matind).                       */
/*      -lb and ub are arrays giving the upper and lower bounds of          */
/*       the variables.                                                     */
/*                                                                          */
/*  int CClp_create (CClp *lp, const char *name)                            */
/*    CREATES an empty lp.  This supports an alternative to CClp_loadlp     */
/*    for loading a problem.                                                */
/*      -name attaches a name to the LP (it can be used by the LP solver    */
/*       in io routines)                                                    */
/*                                                                          */
/*  int CClp_new_row (CClp *lp, char sense, double rhs)                     */
/*    ADDS a new empty row to the lp                                        */
/*      -sense is 'L', 'E', or 'G' for a <=, =, or >= constraint            */
/*      -rhs is the right-hand side of the row                              */
/*                                                                          */
/*  int CClp_change_sense (CClp *lp, int row, char sense)                   */
/*    CHANGES the sense of a row                                            */
/*      -row is the row number to change                                    */
/*      -sense is 'L', 'E', or 'G' to change to <=, =, or >=                */
/*                                                                          */
/*  int CClp_opt (CClp *lp, int method)                                     */
/*    CALLS designated LP solution method.                                  */
/*                                                                          */
/*  int CClp_limited_dualopt (CClp *lp, int lim, int *status,               */
/*      double *upperbound)                                                 */
/*    CALLS the dual simplex method with a limit on the number of pivots.   */
/*      -upperbound it is used to cutoff the dual simplex method (when      */
/*       the objective value reaches upperbound); it can be NULL            */
/*      -status returns the status of the optimization (it can be NULL)     */
/*                                                                          */
/*  int CClp_addrows (CClp *lp, int newrows, int newnz, double *rhs,        */
/*      char *sense, int *rmatbeg, int *rmatind, double *rmatval)           */
/*    ADDS the rows to the LP.                                              */
/*      -newrows is the number of rows to be added                          */
/*      -newnz is the number of nonzero coefficients in the new rows        */
/*      -rhs is an array of the rhs values for the new rows                 */
/*      -sense is 'L', 'E', or 'G' for each of the new rows                 */
/*      -rmatbeg, rmatind, and rmatval give the coefficients of the         */
/*       new rows in sparse format. The arrays can be freed after the       */
/*       call.                                                              */
/*                                                                          */
/*  int CClp_addcols (CClp *lp, int newcols, int newnz, double *obj,        */
/*      int *cmatbeg, int *cmatind, double *cmatval,                        */
/*      double *lb, double *ub)                                             */
/*    ADDS the columns to the LP.                                           */
/*                                                                          */
/*  int CClp_delete_row (CClp *lp, int i)                                   */
/*    DELETES row i of the LP.                                              */
/*                                                                          */
/*  int CClp_delete_set_of_rows (CClp *lp, int *delstat)                    */
/*    DELETES the rows corresponding to 1 entries in delstat.               */
/*      -delstat is a 0/1 array having an entry for each row                */
/*                                                                          */
/*  int CClp_delete_column (CClp *lp, int i)                                */
/*    DELETES column i from the LP.                                         */
/*                                                                          */
/*  int CClp_delete_set_of_columns (CClp *lp, int *delstat)                 */
/*    DELETES the columns corresponding to the 1 entries in delstat.        */
/*      -delstat is a 0/1 array having an entry for each column             */
/*                                                                          */
/*  int CClp_setbnd (CClp *lp, int col, char lower_or_upper, double bnd)    */
/*    SETS the bound on the variable index by col.                          */
/*      -lower_or_upper should be either 'L' or 'U'                         */
/*                                                                          */
/*  int CClp_get_warmstart (CClp *lp, CClp_warmstart **w)                   */
/*    SAVES information for efficiently resolving the current lp in w,      */
/*    for example, basis or norm information                                */
/*                                                                          */
/*  int CClp_load_warmstart (CClp *lp, CClp_warmstart *w)                   */
/*    RESTORES the warmstart information in w.                              */
/*                                                                          */
/*  int CClp_build_warmstart (CClp_warmstart **w, CClp_info *i)             */
/*    BUILDS some warmstart information from the row/column information     */
/*    in i.                                                                 */
/*                                                                          */
/*  void CClp_free_warmstart (CClp_warmstart **w)                           */
/*    FREES the memory used by w.                                           */
/*                                                                          */
/*  int CClp_sread_warmstart (CC_SFILE *f, CClp_warmstart **w)              */
/*    READS warmstart information from the f.                               */
/*                                                                          */
/*  int CClp_swrite_warmstart (CC_SFILE *f, CClp_warmstart *w)              */
/*    WRITES warmstart information from the f.                              */
/*                                                                          */
/*  int CClp_get_info (CClp *lp, CClp_info **i)                             */
/*    BUILDS information useful for efficiently answering questions         */
/*    about the status of rows and columns                                  */
/*                                                                          */
/*  int CClp_create_info (CClp_info **i, int rcount, int ccount)            */
/*    CREATES a structure for storing information about the status of       */
/*    rows and columns.                                                     */
/*                                                                          */
/*  int CClp_is_col_active (CClp_info *i, int c)                            */
/*    returns 1 if column e is active, 0 otherwise.                         */
/*    "active" means participating in the current solution (for example,    */
/*    it could mean basic or nonbasic at upper bound)                       */
/*                                                                          */
/*  int CClp_is_row_active (CClp_info *i, int r)                            */
/*    returns 1 if row e is active, 0 otherwise.                            */
/*    "active" means participating in the current solution (for example,    */
/*    it could mean the row's slack is non-basic)                           */
/*                                                                          */
/*  void CClp_set_col_active (CClp_info *i, int c)                          */
/*    marks column e as active (for eventual CClp_build_warmstart)          */
/*                                                                          */
/*  void CClp_set_col_inactive (CClp_info *i, int c)                        */
/*    marks column e as inactive (for eventual CClp_build_warmstart)        */
/*                                                                          */
/*  void CClp_set_col_upper (CClp_info *i, int c)                           */
/*    marks column e as active at upper bound (for eventual                 */
/*    CClp_build_warmstart)                                                 */
/*                                                                          */
/*  void CClp_set_row_active (CClp_info *i, int r)                          */
/*    marks row r as active (for eventual CClp_build_warmstart)             */
/*                                                                          */
/*  void CClp_set_row_inactive (CClp_info *i, int r)                        */
/*    marks row r as inactive (for eventual CClp_build_warmstart)           */
/*                                                                          */
/*  void CClp_free_info (CClp_info **i)                                     */
/*    FREES the memory used by i.                                           */
/*                                                                          */
/*  int CClp_x (CClp *lp, double *x)                                        */
/*    RETURNS the current LP solution.                                      */
/*      -x should be an array of length at least ncols                      */
/*                                                                          */
/*  int CClp_rc (CClp *lp, double *rc)                                      */
/*    RETURNS the current reduced costs.                                    */
/*      -rc should be an array of length at least ncols                     */
/*                                                                          */
/*  int CClp_pi (CClp *lp, double *pi)                                      */
/*    RETURNS the dual values on the constraints.                           */
/*      -pi should be an array of length at least nrows                     */
/*    NOTES: If the lp and dual lp are feasible, these pi values are        */
/*      the traditional dual solution.  If the dual is unbounded, these     */
/*      pi satisfy                                                          */
/*                                                                          */
/*       pi_i <= 0  for <= constraints                                      */
/*       pi_i >= 0  for >= constraints                                      */
/*                                                                          */
/*       pi'b - sum (pi'A_j * u_j: pi'A_j > 0)                              */
/*            - sum (pi'A_j * l_j: pi'A_j < 0) > 0                          */
/*                                                                          */
/*    where b is the rhs vector, u_j is the upper bound on variable x_j,    */
/*    l_j the lower bound, and A_j the constraint matrix column for x_j.    */
/*                                                                          */
/*  int CClp_objval (CClp *lp, double *obj)                                 */
/*    RETURNS the objective value of the lp.                                */
/*                                                                          */
/*  int CClp_nrows (CClp *lp)                                               */
/*    RETURNS the number of rows in the LP.                                 */
/*                                                                          */
/*  int CClp_ncols (CClp *lp)                                               */
/*    RETURNS the number of columns in the LP.                              */
/*                                                                          */
/*  int CClp_nnonzeros (CClp *lp)                                           */
/*    RETURNS the number of nonzeros in the LP.                             */
/*                                                                          */
/*  int CClp_status (CClp *lp, int *status)                                 */
/*    CHECKS whether the current lp is infeasible or whether an optimal     */
/*     solution has been found. It returns an error if the LP has not       */
/*     not been optimized.                                                  */
/*      -lp is the lp                                                       */
/*      -status returns 0 if the lp has an optimal solution and 1 if it     */
/*       is infeasible.                                                     */
/*                                                                          */
/*  int CClp_getweight (CClp *lp, int nrows, int *rmatbeg, int *rmatind,    */
/*      double *rmatval, double *weight)                                    */
/*    NOTE: NOT IMPLEMENTED                                                 */
/*    COMPUTES the duals of the steepest edge norms for the n rows          */
/*     specified in rmatbeg, rmatind, and rmatval.                          */
/*      -weight returns the array of weights; the array should be at        */
/*       least nrows long                                                   */
/*                                                                          */
/*  int CClp_dump_lp (CClp *lp, const char *fname)                          */
/*    WRITES the LP to file fname.                                          */
/*                                                                          */
/*  int CClp_getgoodlist (CClp *lp, int *goodlist, int *goodlen_p,          */
/*      double *downpen, double *uppen)                                     */
/*    RETURNS an array of the column indices corresponding to variables     */
/*     that move in both the up and down directions. This is a useful       */
/*     list of candidates for strong branching.                             */
/*      -goodlist, downpen and uppen should be arrays of length at          */
/*       least ncols.                                                       */
/*                                                                          */
/*  int CClp_strongbranch (CClp *lp, int *candidatelist, int ncand,         */
/*      double *downpen, double *uppen, int iterations,                     */
/*      double upperbound)                                                  */
/*    RETURNS estimates of the lp values obtained by setting each of the    */
/*      ncand variables listed in candidatelist to 0 and 1. The estimates   */
/*      are obtained by performing iterations pivots of dual simplex        */
/*      method. upperbound is used to cutoff the dual simplex method.       */
/*      downpen and uppen should never be > upperbound.                     */
/*       -downpen and uppen should be arrays of length at least ncand       */
/*                                                                          */
/*                                                                          */
/****************************************************************************/


#include <math.h>

#include "lp.h"
#include "util.h"
#include "xprs.h"

#undef CC_XPRESS_DISPLAY

#undef XP_DEBUG
#ifdef XP_DEBUG
#define XP_IF_DEBUG
#else
#define XP_IF_DEBUG if (0)
#endif

/* Without macros, error handling would likely make up more than 80% of the
 * source code. While the macros add some complexity, they significantly reduce
 * the amount of duplicated blocks.*/

#define SOLVER_WARMSTART_NAME "CPL5"

#define _XP_STR(X)   _XP__STR(X)
#define _XP__STR(X)  #X
#define _XP_LINE_STR _XP_STR(__LINE__)

/* Current location in source file */
#define XP_LOC_STR __FILE__ ":" _XP_LINE_STR " > "

/* Classic min ternary operator */
#define XP_MIN(A, B) ((A) < (B) ? (A) : (B))

/* Check if EXPR return a non zero value, if it does go to CLEANUP */
#define XP_CHK(EXPR)                                                       \
    do {                                                                   \
        rval = (EXPR);                                                     \
        if ((rval) != 0) {                                                 \
            fprintf(stderr, XP_LOC_STR "'" #EXPR "' returned %d\n", rval); \
            goto CLEANUP;                                                  \
        }                                                                  \
    } while (0)

/* During cleanup we cannot go to CLEANUP again, so we cumulate errors */
#define XP_CUMUL_CHK(EXPR)                                                 \
    do {                                                                   \
        rval += (EXPR);                                                    \
        if ((rval) != 0)                                                   \
            fprintf(stderr,                                                \
                    XP_LOC_STR "'" #EXPR "' returned %d, continuing...\n", \
                    rval);                                                 \
    } while (0)

/* Allocate, zero and check some space in memory */
#define XP_CHK_CALLOC(PTR, COUNT, TYPE)                                     \
    do {                                                                    \
        (PTR) = CC_SAFE_MALLOC(COUNT, TYPE);                                \
        if ((PTR) == (TYPE*)NULL) {                                         \
            fprintf(stderr,                                                 \
                    XP_LOC_STR "Out of memory: " #PTR " = " #TYPE "[%d]\n", \
                    COUNT);                                                 \
            rval = 1;                                                       \
            goto CLEANUP;                                                   \
        }                                                                   \
        memset(PTR, 0, (COUNT) * sizeof(TYPE));                             \
    } while (0)

/* Work with Xpress controls, OP={get, set}, TYPE={int, dbl} */
#define XP_DO_CTRL(OP, TYPE, CTRL, VALUE) \
    XP_CHK(XPRS##OP##TYPE##control(lp->xp_prob, CTRL, VALUE))

/* Get Xpress attribute */
#define XP_GET_ATTR(TYPE, CTRL, VALUE) \
    XP_CHK(XPRSget##TYPE##attrib(lp->xp_prob, CTRL, VALUE))

/* During cleanup we cannot go to CLEANUP again, so we cumulate errors */
#define XP_CUMUL_CTRL(OP, TYPE, CTRL, VALUE) \
    XP_CUMUL_CHK(XPRS##OP##TYPE##control(lp->xp_prob, CTRL, VALUE))

/****************************************************************************/
/*    EXPORTED TYPE DEFINITIONS                                             */
/****************************************************************************/

struct CClp {
    XPRSprob xp_prob;
};

struct CClp_info {
    int  rcount;
    int  ccount;
    int* rstat;
    int* cstat;
};

struct CClp_warmstart {
    struct CClp_info i;
};

/****************************************************************************/
/*    PRIVATE HELPER FUNCTIONS                                              */
/****************************************************************************/

static int  priv_check_lp_status(CClp* lp);
static int  priv_print_stop_reason(CClp* lp);
static int  priv_check_sol_status(CClp* lp, int* sol_status, int silent);
static int  priv_primalopt(CClp* lp);
static int  priv_dualopt(CClp* lp);
static int  priv_baropt(CClp* lp);
static int  priv_lpopt(CClp* lp, char pd, char pd_backup);
static void priv_free_info(CClp_info* i);
static int  priv_init_info(CClp_info* i, int ccount, int rcount);
static void priv_messagecb(XPRSprob cbprob, void* cbdata, char const* msg,
                           int len, int msgtype);
static int  priv_extract_nz(int const* vec, int vcnt, int** nz_vec, int* nzcnt);

static int priv_check_lp_status(CClp* lp) {
    int rval = 0;

    int solve_status = 0;
    XP_GET_ATTR(int, XPRS_SOLVESTATUS, &solve_status);

    if (solve_status != XPRS_SOLVESTATUS_COMPLETED) {
        fprintf(stderr, "Xpress optimization not completed\n");

        if (solve_status == XPRS_SOLVESTATUS_FAILED)
            fprintf(stderr, "Xpress optimization failed\n");
        if (solve_status == XPRS_SOLVESTATUS_STOPPED)
            XP_CHK(priv_print_stop_reason(lp));
        if (solve_status == XPRS_SOLVESTATUS_UNSTARTED)
            fprintf(stderr, "Xpress optimization not started\n");
        rval = 1;
    }

CLEANUP:
    return rval;
}

static int priv_print_stop_reason(CClp* lp) {
    int rval        = 0;
    int stop_status = 0;

    fprintf(stderr, "Xpress optimization was stopped\n");
    XP_GET_ATTR(int, XPRS_STOPSTATUS, &stop_status);
    if (stop_status == XPRS_STOP_NONE)
        fprintf(stderr, "Xpress no interruption - the solve completed "
                        "normally\n");
    if (stop_status == XPRS_STOP_TIMELIMIT)
        fprintf(stderr, "Xpress time limit hit\n");
    if (stop_status == XPRS_STOP_CTRLC)
        fprintf(stderr, "Xpress sigint received\n");
    if (stop_status == XPRS_STOP_NODELIMIT)
        fprintf(stderr, "Xpress node limit hit in branch-and-bound\n");
    if (stop_status == XPRS_STOP_ITERLIMIT)
        fprintf(stderr, "Xpress LP iteration limit hit\n");
    if (stop_status == XPRS_STOP_MIPGAP)
        fprintf(stderr, "Xpress MIP gap is sufficiently small\n");
    if (stop_status == XPRS_STOP_SOLLIMIT)
        fprintf(stderr, "Xpress solution limit hit\n");
    if (stop_status == XPRS_STOP_MEMORYERROR)
        fprintf(stderr, "Xpress insufficient memory\n");
    if (stop_status == XPRS_STOP_USER)
        fprintf(stderr, "Xpress user interrupt\n");
    if (stop_status == XPRS_STOP_SOLVECOMPLETE)
        fprintf(stderr, "Xpress MIP TREE not fully explored\n");
    if (stop_status == XPRS_STOP_LICENSELOST)
        fprintf(stderr, "Xpress license has been lost\n");
    if (stop_status == XPRS_STOP_NUMERICALERROR)
        fprintf(stderr, "Xpress solver encountered numerical errors\n");

CLEANUP:
    return rval;
}

static int priv_check_sol_status(CClp* lp, int* sol_status, int silent) {
    int rval   = 0;
    int status = 0;

    XP_GET_ATTR(int, XPRS_LPSTATUS, &status);

    if (!silent && status != XPRS_LP_OPTIMAL) {
        if (status == XPRS_LP_INFEAS)
            fprintf(stderr, "Xpress LP: Infeasible\n");
        else if (status == XPRS_LP_UNBOUNDED)
            fprintf(stderr, "Xpress LP: Unbounded\n");
        else if (status == XPRS_LP_UNSTARTED)
            fprintf(stderr, "Xpress LP: Unstarted\n");
        else if (status == XPRS_LP_CUTOFF)
            fprintf(stderr, "Xpress LP: Objective worse than cutoff\n");
        else if (status == XPRS_LP_UNFINISHED)
            fprintf(stderr, "Xpress LP: Unfinished\n");
        else if (status == XPRS_LP_CUTOFF_IN_DUAL)
            fprintf(stderr, "Xpress LP: Cutoff in dual\n");
        else if (status == XPRS_LP_UNSOLVED)
            fprintf(stderr, "Xpress LP: Problem encountered numerical "
                            "issues\n");
        else if (status == XPRS_LP_NONCONVEX)
            fprintf(stderr, "Xpress LP: Problem not convex\n");
    }

    if (sol_status != (int*)NULL)
        *sol_status = status;

CLEANUP:
    return rval;
}

static int priv_lpopt(CClp* lp, char pd, char pd_backup) {
    int rval      = 0;
    int solstatus = 0;

#ifdef XP_DEBUG
    double obj;
    int    ncols = CClp_ncols(lp);
    int    nrows = CClp_nrows(lp);
#endif

    XP_CHK(XPRSlpoptimize(lp->xp_prob, &pd));
    XP_CHK(priv_check_lp_status(lp));
    XP_CHK(priv_check_sol_status(lp, &solstatus, 1));
    if ((solstatus == XPRS_LP_INFEAS || solstatus == XPRS_LP_UNBOUNDED) &&
        pd_backup != '\0')
        XP_CHK(priv_lpopt(lp, pd_backup, '\0'));

#ifdef XP_DEBUG
    XP_CHK(XPRSgetdblattrib(lp->xp_prob, XPRS_LPOBJVAL, &obj));
    printf("CClp_opt: obj = %f (%d x %d)\n", obj, nrows, ncols);
#endif

CLEANUP:
    return rval;
}

static int priv_primalopt(CClp* lp) {
    return priv_lpopt(lp, 'p', '\0');
}

static int priv_dualopt(CClp* lp) {
    return priv_lpopt(lp, 'd', '\0');
}

static int priv_baropt(CClp* lp) {
    return priv_lpopt(lp, 'b', 'd');
}

static void priv_free_info(CClp_info* i) {
    if (i != (CClp_info*)NULL) {
        CC_IFFREE(i->cstat, int);
        CC_IFFREE(i->rstat, int);
    }
}

static int priv_init_info(CClp_info* i, int ccount, int rcount) {
    int rval = 0;

    if (ccount == 0 || rcount == 0) {
        fprintf(stderr, "No row or columns to intialized LP info\n");
        rval = 1;
        goto CLEANUP;
    }

    priv_free_info(i);
    i->ccount = ccount;
    i->rcount = rcount;

    XP_CHK_CALLOC(i->cstat, i->ccount, int);
    XP_CHK_CALLOC(i->rstat, i->rcount, int);

    return 0;

CLEANUP:
    priv_free_info(i);
    return rval;
}

static void priv_messagecb(XPRSprob cbprob, void* cbdata, char const* msg,
                           int len, int msgtype) {
    (void)cbprob; /* unused */
    (void)cbdata; /* unused */

    if (msgtype == 4) {
        fprintf(stderr, "XPR> %*s\n", len, msg);
        fflush(stdout);
    }

    XP_IF_DEBUG {
        if (msgtype == 3) {
            fprintf(stderr, "XPR> %*s\n", len, msg);
            fflush(stdout);
        }
    }

#ifdef CC_XPRESS_DISPLAY
    if (msgtype == 1 || msgtype == 2) {
        fprintf(stdout, "XPR> %*s\n", len, msg);
        fflush(stdout);
    }
#endif
}

static int priv_extract_nz(int const* vec, int vcnt, int** nz_vec, int* nzcnt) {
    int rval = 0, i = 0, j = 0;

    CC_IFFREE(*nz_vec, int);
    *nzcnt = 0;

    for (i = 0; i < vcnt; i++)
        if (vec[i])
            (*nzcnt)++;

    if (*nzcnt == 0)
        return 0;

    XP_CHK_CALLOC(*nz_vec, *nzcnt, int);
    for (i = 0, j = 0; i < vcnt; i++)
        if (vec[i])
            (*nz_vec)[j++] = i;

CLEANUP:
    return rval;
}

/****************************************************************************/
/*    EXPORTED FUNCTIONS                                                    */
/****************************************************************************/

int xprs_prob_count = 0; /* Global Xpress problem counter */

int CClp_init(CClp** lp_ptr) {
    int  rval = 0;
    char message[512];

    CClp_free(lp_ptr);
    XP_CHK_CALLOC(*lp_ptr, 1, CClp);
    if (xprs_prob_count++ == 0) {
        rval = XPRSinit("");
        if (rval != 0) {
            fprintf(stderr, "XPRSinit failed, return code %d\n", rval);
            XPRSgetlicerrmsg(message, sizeof(message));
            fprintf(stderr, "Licensing error: %s\n", message);
            goto CLEANUP;
        }
    }

    XP_CHK(XPRScreateprob(&(*lp_ptr)->xp_prob));

CLEANUP:
    return rval;
}

void CClp_freelp(CClp** lp_ptr) {
    CClp* lp = *lp_ptr;

    if (lp == (CClp*)NULL)
        return;

    if (lp->xp_prob != (XPRSprob)NULL)
        XPRSdestroyprob(lp->xp_prob);
    lp->xp_prob = (XPRSprob)NULL;
    --xprs_prob_count;
}

void CClp_free(CClp** lp_ptr) {

    CClp_freelp(lp_ptr);
    CC_IFFREE(*lp_ptr, CClp);

    if (xprs_prob_count == 0)
        XPRSfree();
}

int CClp_create(CClp* lp, char const* name) {
    int rval = 0;

    if (lp->xp_prob == (XPRSprob)NULL)
        XP_CHK(XPRScreateprob(&lp->xp_prob));
    XP_CHK(XPRSsetprobname(lp->xp_prob, name));
    XP_CHK(lp->xp_prob == NULL ? 2 : 0);

    XP_CHK(XPRSaddcbmessage(lp->xp_prob, priv_messagecb, NULL, 0));
#ifdef CC_XPRESS_DISPLAY
    XP_DO_CTRL(set, int, XPRS_OUTPUTLOG, XPRS_OUTPUTLOG_FULL_OUTPUT);
    XP_DO_CTRL(set, int, XPRS_LPLOGSTYLE, 1);
#endif

    XP_DO_CTRL(set, int, XPRS_KEEPBASIS, 1);
    XP_DO_CTRL(set, int, XPRS_DUALGRADIENT, 1);
    XP_DO_CTRL(set, int, XPRS_PRICINGALG, 1);
    XP_DO_CTRL(set, dbl, XPRS_OPTIMALITYTOL, 1.0E-9);
    XP_DO_CTRL(set, dbl, XPRS_FEASTOL, 1.0E-9);
    XP_DO_CTRL(set, int, XPRS_RANDOMSEED, 0);

CLEANUP:
    return rval;
}

int CClp_loadlp(CClp* lp, char const* name, int ncols, int nrows, int objsense,
                double* obj, double* rhs, char* sense, int* matbeg, int* matcnt,
                int* matind, double* matval, double* lb, double* ub) {

    int     rval = 0;
    double* rng  = (double*)NULL;

    XP_CHK(CClp_create(lp, name));
    XP_CHK(XPRSloadlp(lp->xp_prob, name, ncols, nrows, sense, rhs, rng, obj,
                      matbeg, matcnt, matind, matval, lb, ub));
    XP_CHK(XPRSchgobjsense(lp->xp_prob, objsense));

CLEANUP:
    return rval;
}

int CClp_force_perturb(CClp* lp) {
    int rval = 0;

    XP_DO_CTRL(set, int, XPRS_AUTOPERTURB, 1);

CLEANUP:
    return 0;
}

int CClp_tune_small(CClp* lp) {
    int rval = 0;

    XP_DO_CTRL(set, int, XPRS_DUALGRADIENT, 0);
    XP_DO_CTRL(set, int, XPRS_PRICINGALG, 0);

CLEANUP:
    return rval;
}

int CClp_disable_presolve(CClp* lp) {
    int rval = 0;

    XP_DO_CTRL(set, int, XPRS_PRESOLVE, 0);

CLEANUP:
    return rval;
}

int CClp_new_row(CClp* lp, char sense, double rhs) {
    int*    start   = (int*)NULL;
    int*    colind  = (int*)NULL;
    double* rng     = (double*)NULL;
    double* rowcoef = (double*)NULL;

    return XPRSaddrows(lp->xp_prob, 1, 0, &sense, &rhs, rng, start, colind,
                       rowcoef);
}

int CClp_change_sense(CClp* lp, int row, char sense) {
    return XPRSchgrowtype(lp->xp_prob, 1, &row, &sense);
}

int CClp_opt(CClp* lp, int method) {
    switch (method) {
    case CClp_METHOD_PRIMAL:
        return priv_primalopt(lp);
    case CClp_METHOD_DUAL:
        return priv_dualopt(lp);
    case CClp_METHOD_BARRIER:
        return priv_baropt(lp);
    default:
        fprintf(stderr, "Nonexistent method in CClp_opt\n");
    }
    return 1;
}

#define CC_MAX_REFACTORFREQ 150

int CClp_limited_dualopt(CClp* lp, int iterationlim, int* status,
                         double* upperbound) {
    int    rval            = 0;
    int    cleanup_level   = 0;
    int    solstatus       = 0;
    int    old_autoperturb = 0;
    int    old_lpiterlimit = 0;
    double old_abscutoff   = 0;
    int    old_presolve    = 0;
    int    old_invertfreq  = 0;

    XP_DO_CTRL(get, int, XPRS_AUTOPERTURB, &old_autoperturb);
    XP_DO_CTRL(get, int, XPRS_LPITERLIMIT, &old_lpiterlimit);
    XP_DO_CTRL(get, dbl, XPRS_MIPABSCUTOFF, &old_abscutoff);
    XP_DO_CTRL(get, int, XPRS_PRESOLVE, &old_presolve);
    XP_DO_CTRL(get, int, XPRS_INVERTFREQ, &old_invertfreq);

    XP_DO_CTRL(set, int, XPRS_AUTOPERTURB, 0);
    ++cleanup_level;
    XP_DO_CTRL(set, int, XPRS_LPITERLIMIT, iterationlim);
    ++cleanup_level;
    XP_DO_CTRL(set, dbl, XPRS_MIPABSCUTOFF, *upperbound);
    ++cleanup_level;
    XP_DO_CTRL(set, int, XPRS_PRESOLVE, 0);
    ++cleanup_level;
    XP_DO_CTRL(set, int, XPRS_INVERTFREQ, iterationlim + 1);
    ++cleanup_level;

    XP_CHK(XPRSlpoptimize(lp->xp_prob, "d"));
    XP_CHK(priv_check_sol_status(lp, &solstatus, 1));

    if (status == (int*)NULL)
        goto CLEANUP;

    if (solstatus == XPRS_LP_OPTIMAL)
        *status = CClp_SUCCESS;
    else if (solstatus == XPRS_LP_INFEAS)
        *status = CClp_INFEASIBLE;
    else
        *status = CClp_UNKNOWN;

CLEANUP:
    if (cleanup_level >= 5)
        XP_CUMUL_CTRL(set, int, XPRS_INVERTFREQ, old_invertfreq);
    if (cleanup_level >= 4)
        XP_CUMUL_CTRL(set, int, XPRS_PRESOLVE, old_presolve);
    if (cleanup_level >= 3)
        XP_CUMUL_CTRL(set, dbl, XPRS_MIPABSCUTOFF, old_abscutoff);
    if (cleanup_level >= 2)
        XP_CUMUL_CTRL(set, int, XPRS_LPITERLIMIT, old_lpiterlimit);
    if (cleanup_level >= 1)
        XP_CUMUL_CTRL(set, int, XPRS_AUTOPERTURB, old_autoperturb);

    if (rval && status)
        *status = CClp_FAILURE;

    return rval;
}

int CClp_addrows(CClp* lp, int newrows, int newnz, double* rhs, char* sense,
                 int* rmatbeg, int* rmatind, double* rmatval) {
    double* rng = (double*)NULL;
    return XPRSaddrows(lp->xp_prob, newrows, newnz, sense, rhs, rng, rmatbeg,
                       rmatind, rmatval);
}

int CClp_addcols(CClp* lp, int newcols, int newnz, double* obj, int* cmatbeg,
                 int* cmatind, double* cmatval, double* lb, double* ub) {
    return XPRSaddcols(lp->xp_prob, newcols, newnz, obj, cmatbeg, cmatind,
                       cmatval, lb, ub);
}

int CClp_delete_row(CClp* lp, int i) {
    int rval = 0;

    /* TODO(cava): Manually adjust basis */
    XP_CHK(XPRSdelrows(lp->xp_prob, 1, &i));

CLEANUP:
    return rval;
}

int CClp_delete_set_of_rows(CClp* lp, int* delstat) {
    int  rval    = 0;
    int  rcnt    = CClp_nrows(lp);
    int* dellist = (int*)NULL;
    int  delcnt  = 0;

    priv_extract_nz(delstat, rcnt, &dellist, &delcnt);
    if (delcnt == 0)
        return 0;

    /* TODO(cava): Manually adjust basis */
    XP_CHK(XPRSdelrows(lp->xp_prob, delcnt, dellist));

CLEANUP:
    CC_IFFREE(dellist, int);
    return rval;
}

int CClp_delete_column(CClp* lp, int i) {
    int rval = 0;

    /*  TODO(cava): find another column to enter the basis replacing i */
    XP_CHK(XPRSdelcols(lp->xp_prob, 1, &i));

CLEANUP:
    return rval;
}

int CClp_delete_set_of_columns(CClp* lp, int* delstat) {
    int  rval    = 0;
    int  ccnt    = CClp_ncols(lp);
    int* dellist = (int*)NULL;
    int  delcnt  = 0;

    priv_extract_nz(delstat, ccnt, &dellist, &delcnt);
    if (delcnt == 0)
        return 0;

    /* TODO(cava): Manually adjust basis */
    XP_CHK(XPRSdelcols(lp->xp_prob, delcnt, dellist));

CLEANUP:
    CC_IFFREE(dellist, int);

    return rval;
}

int CClp_setbnd(CClp* lp, int col, char lower_or_upper, double bnd) {
    return XPRSchgbounds(lp->xp_prob, 1, &col, &lower_or_upper, &bnd);
}

int CClp_get_warmstart(CClp* lp, CClp_warmstart** w) {
    int        rval  = 0;
    CClp_info* winfo = (CClp_info*)NULL;

    if (*w == (CClp_warmstart*)NULL)
        XP_CHK_CALLOC(*w, 1, CClp_warmstart);

    winfo = &(*w)->i;
    XP_CHK(priv_init_info(winfo, CClp_ncols(lp), CClp_nrows(lp)));
    XP_CHK(XPRSgetbasis(lp->xp_prob, winfo->rstat, winfo->cstat));

    return 0;

CLEANUP:
    CClp_free_warmstart(w);
    return rval;
}

int CClp_load_warmstart(CClp* lp, CClp_warmstart* w) {
    int rval = 0;

    if (w->i.cstat != (int*)NULL && w->i.rstat != (int*)NULL) {
        XP_CHK(XPRSloadbasis(lp->xp_prob, w->i.rstat, w->i.cstat));
    } else {
        printf("WARNING: No basis in call to load_warmstart\n");
        fflush(stdout);
    }

CLEANUP:
    return rval;
}

int CClp_build_warmstart(CClp_warmstart** w, CClp_info* i) {
    int        rval  = 0;
    int        j     = 0;
    CClp_info* winfo = (CClp_info*)NULL;

    if (*w == (CClp_warmstart*)NULL)
        XP_CHK_CALLOC(*w, 1, CClp_warmstart);

    winfo = &(*w)->i;
    XP_CHK(priv_init_info(winfo, i->ccount, i->rcount));

    for (j = 0; j < i->ccount; j++)
        winfo->cstat[j] = i->cstat[j];
    for (j = 0; j < i->rcount; j++)
        winfo->rstat[j] = i->rstat[j];

    return 0;

CLEANUP:
    CClp_free_warmstart(w);
    return rval;
}

void CClp_free_warmstart(CClp_warmstart** w) {

    if (w == (CClp_warmstart**)NULL || *w == (CClp_warmstart*)NULL)
        return;
    priv_free_info(&(*w)->i);
    CC_FREE(*w, CClp_warmstart);
}

int CClp_sread_warmstart(CC_SFILE* f, CClp_warmstart** w) {
    int        rval       = 0;
    char       name[5]    = {0};
    int        i          = 0;
    int        ccount     = 0;
    int        rcount     = 0;
    int        has_dnorms = 0;
    CClp_info* winfo      = (CClp_info*)NULL;

    for (i = 0; i < 4; i++)
        XP_CHK(CCutil_sread_char(f, name + i));

    if (strncmp(name, SOLVER_WARMSTART_NAME, 4) != 0) {
        fprintf(stderr, "warmstart for another solver (%s) ignored\n", name);
        return 0;
    }

    XP_CHK(CCutil_sread_int(f, &ccount));
    XP_CHK(CCutil_sread_int(f, &rcount));

    if (*w == (CClp_warmstart*)NULL)
        XP_CHK_CALLOC(*w, 1, CClp_warmstart);
    winfo = &(*w)->i;
    XP_CHK(priv_init_info(winfo, ccount, rcount));

    for (i = 0; i < ccount; i++)
        XP_CHK(CCutil_sread_bits(f, winfo->cstat + i, 2));
    for (i = 0; i < rcount; i++)
        XP_CHK(CCutil_sread_bits(f, winfo->rstat + i, 1));

    XP_CHK(CCutil_sread_int(f, &has_dnorms));
    return 0;

CLEANUP:
    CClp_free_warmstart(w);
    return rval;
}

int CClp_swrite_warmstart(CC_SFILE* f, CClp_warmstart* w) {
    int        rval = 0;
    int        i    = 0;
    CClp_info* info = &(w->i);

    for (i = 0; i < 4; i++)
        XP_CHK(CCutil_swrite_char(f, SOLVER_WARMSTART_NAME[i]));

    XP_CHK(CCutil_swrite_int(f, info->ccount));
    XP_CHK(CCutil_swrite_int(f, info->rcount));

    for (i = 0; i < info->ccount; i++)
        XP_CHK(CCutil_swrite_bits(f, info->cstat[i], 2));

    for (i = 0; i < info->rcount; i++)
        XP_CHK(CCutil_swrite_bits(f, info->rstat[i], 1));
    XP_CHK(CCutil_swrite_int(f, 0));

CLEANUP:
    return rval;
}

void CClp_free_info(CClp_info** i) {

    if (i == (CClp_info**)NULL || *i == (CClp_info*)NULL)
        return;
    priv_free_info(*i);
    CC_FREE(*i, CClp_info);
}

int CClp_get_info(CClp* lp, CClp_info** i) {
    int rval = 0;

    if (*i == (CClp_info*)NULL)
        XP_CHK_CALLOC(*i, 1, CClp_info);
    XP_CHK(priv_init_info(*i, CClp_ncols(lp), CClp_nrows(lp)));
    XP_CHK(XPRSgetbasis(lp->xp_prob, (*i)->rstat, (*i)->cstat));

    return 0;

CLEANUP:
    CClp_free_info(i);
    return rval;
}

int CClp_create_info(CClp_info** i, int rcount, int ccount) {
    int rval = 0;
    int j    = 0;

    if (*i == (CClp_info*)NULL)
        XP_CHK_CALLOC(*i, 1, CClp_info);
    XP_CHK(priv_init_info(*i, ccount, rcount));

    for (j = 0; j < ccount; j++)
        (*i)->cstat[j] = 0;
    for (j = 0; j < rcount; j++)
        (*i)->rstat[j] = 0;

    return 0;

CLEANUP:
    CClp_free_info(i);
    return rval;
}

int CClp_is_col_active(CClp_info* i, int c) {
    if (c < 0 || c >= i->ccount)
        return 0;
    return i->cstat[c] == 1 || i->cstat[c] == 2;
}

int CClp_is_row_active(CClp_info* i, int c) {
    if (c < 0 || c >= i->rcount)
        return 0;
    return i->rstat[c] == 0 || i->rstat[c] == 2;
}

void CClp_set_col_active(CClp_info* i, int c) {
    if (c >= 0 && c < i->ccount)
        i->cstat[c] = 1; /*variable is basic*/
}

void CClp_set_col_inactive(CClp_info* i, int c) {
    if (c >= 0 && c < i->ccount)
        i->cstat[c] = 0; /*variable at lower bound*/
}

void CClp_set_col_upper(CClp_info* i, int c) {
    if (c >= 0 && c < i->ccount)
        i->cstat[c] = 2; /*variable at upper bound*/
}

void CClp_set_row_active(CClp_info* i, int r) {
    if (r >= 0 && r < i->rcount)
        i->rstat[r] = 2;
    /*associated slack/surplus/artificial-var is nonbasic at value 0.0 */
}

void CClp_set_row_inactive(CClp_info* i, int r) {
    if (r >= 0 && r < i->rcount)
        i->rstat[r] = 1;
    /*associated slack/surplus/artificial-var is basic*/
}

int CClp_x(CClp* lp, double* x) {
    int rval   = 0;
    int status = 0;

    int ncols = CClp_ncols(lp);
    if (ncols == 0)
        return 0;

    XP_CHK(XPRSgetsolution(lp->xp_prob, &status, x, 0, ncols - 1));
    XP_CHK(status == XPRS_LP_INFEAS ? status : 0);

CLEANUP:
    return rval;
}

int CClp_rc(CClp* lp, double* rc) {
    int rval   = 0;
    int status = 0;

    int ncols = CClp_ncols(lp);
    if (ncols == 0)
        return 0;

    XP_CHK(XPRSgetredcosts(lp->xp_prob, &status, rc, 0, ncols - 1));

CLEANUP:
    return rval;
}

int CClp_pi(CClp* lp, double* pi) {
    int rval    = 0;
    int nrows   = 0;
    int method  = 0;
    int status  = 0;
    int has_ray = 0;

    XP_CHK(XPRSgetintcontrol(lp->xp_prob, XPRS_LPFLAGS, &method));
    XP_CHK(XPRSgetintattrib(lp->xp_prob, XPRS_LPSTATUS, &status));

    if (method == XPRS_LPFLAGS_DUAL && status == XPRS_LP_INFEAS) {
        XP_CHK(XPRSgetdualray(lp->xp_prob, pi, &has_ray));
        XP_CHK(!has_ray ? 1 : 0);
    } else {
        nrows = CClp_nrows(lp);
        XP_CHK(nrows == 0 ? 1 : 0);
        XP_CHK(XPRSgetduals(lp->xp_prob, &status, pi, 0, nrows - 1));
    }

CLEANUP:
    return rval;
}

int CClp_objval(CClp* lp, double* obj) {
    return XPRSgetdblattrib(lp->xp_prob, XPRS_LPOBJVAL, obj);
}

/* The following functions are not expected to fail */
int CClp_nrows(CClp* lp) {
    int res = 0;

    if (XPRSgetintattrib(lp->xp_prob, XPRS_ORIGINALROWS, &res)) {
        fprintf(stderr, "XPRSgetintattrib XPRS_ORIGINALROWS failed\n");
        return 0;
    }
    return res;
}

int CClp_ncols(CClp* lp) {
    int res = 0;

    if (XPRSgetintattrib(lp->xp_prob, XPRS_ORIGINALCOLS, &res)) {
        fprintf(stderr, "XPRSgetintattrib XPRS_ORIGINALCOLS failed\n");
        return 0;
    }
    return res;
}

int CClp_nnonzeros(CClp* lp) {
    int is_presolved = 0;
    int res          = 0;

    if (XPRSgetintattrib(lp->xp_prob, XPRS_PRESOLVESTATE, &is_presolved)) {
        fprintf(stderr, "XPRSgetintattrib XPRS_PRESOLVESTATE failed\n");
        return 0;
    }

    if (is_presolved)
        fprintf(stderr, "XPRS_ELEMS: returned non-zeros of presolved model\n");

    if (XPRSgetintattrib(lp->xp_prob, XPRS_ELEMS, &res)) {
        fprintf(stderr, "XPRSgetintattrib XPRS_ELEMS failed\n");
        return 0;
    }

    return res;
}

int CClp_status(CClp* lp, int* status) {
    int rval      = 0;
    int solmethod = 0;
    int solstat   = 0;

    XP_DO_CTRL(get, int, XPRS_LPFLAGS, &solmethod);
    if (solmethod == XPRS_LPFLAGS_DUAL || solmethod == XPRS_LPFLAGS_PRIMAL) {
        fprintf(stderr, "lp not solved by usual methods: %d\n", solmethod);
        *status = -2;
        return 1;
    }

    XP_CHK(XPRSgetintattrib(lp->xp_prob, XPRS_LPSTATUS, &solstat));
    if (solstat == XPRS_LP_OPTIMAL) {
        *status = 0;
    } else if (solstat == XPRS_LP_INFEAS && solmethod == XPRS_LPFLAGS_DUAL) {
        *status = 1;
    } else {
        fprintf(stderr, "lp in an unknown state: %d %d\n", solmethod, solstat);
        *status = -1;
        return 1;
    }

CLEANUP:
    return rval;
}

int CClp_getweight(CClp* lp, int nrows, int* rmatbeg, int* rmatind,
                   double* rmatval, double* weight) {

    XP_IF_DEBUG printf("CClp_getweight not supported with Xpress.\n");

    if (lp || nrows || rmatbeg || rmatind || rmatval || weight)
        return 1;
    return 0;
}

int CClp_dump_lp(CClp* lp, char const* fname) {
    return XPRSwriteprob(lp->xp_prob, fname, "lp");
}

#define OURXPRESSZERO    (1.0E-10)
#define OURXPRESS_INTTOL (0.0001)

int CClp_getgoodlist(CClp* lp, int* goodlist, int* goodlen_p, double* downpen,
                     double* uppen) {
    int     rval   = 0;
    int     ncols  = 0;
    int     i      = 0;
    int     j      = 0;
    int     status = 0;
    int*    cstat  = (int*)NULL;
    double* x      = (double*)NULL;

    /* Call XPRSlpoptimize and verify optimality */
    XP_CHK(XPRSlpoptimize(lp->xp_prob, "d"));
    ncols = CClp_ncols(lp);
    XP_CHK(ncols == 0 ? 1 : 0);

    XP_CHK_CALLOC(x, ncols, double);
    XP_CHK(XPRSgetsolution(lp->xp_prob, &status, x, 0, ncols - 1));
    XP_CHK(status == XPRS_LP_INFEAS ? status : 0);

    XP_CHK_CALLOC(cstat, ncols, int);
    XP_CHK(XPRSgetbasis(lp->xp_prob, (int*)NULL, cstat));

    *goodlen_p = 0;
    for (i = 0, j = 0; i < ncols; i++) {
        if (cstat[i] == 1 && x[i] >= OURXPRESS_INTTOL &&
            x[i] <= 1.0 - OURXPRESS_INTTOL) {
            goodlist[j] = i;
            downpen[j]  = x[i];
            uppen[j]    = 1.0 - x[i];
            j++;
        }
    }

    *goodlen_p = j;

    /* TODO(cava): compute Driebeek penalties and do the same thing that is done
     * in lpcplex8.c (can be efficiently computed by ratio tests) */

CLEANUP:
    CC_IFFREE(cstat, int);
    CC_IFFREE(x, double);
    return rval;
}

int CClp_strongbranch(CClp* lp, int* candidatelist, int ncand, double* downpen,
                      double* uppen, int iterations, double upperbound) {

    int     rval            = 0;
    double  old_upperbound  = .0;
    int     old_autoperturb = 0;
    int     cleanup_level   = 0;
    int     i               = 0;
    int     cand            = 0;
    int     status          = 0;
    int     ncols           = 0;
    double* x               = (double*)NULL;
    int*    candlist_twice  = (int*)NULL;
    char*   boundstypes     = (char*)NULL;
    double* bounds          = (double*)NULL;
    double* objs            = (double*)NULL;
    int*    statuses        = (int*)NULL;

    /* Does it limit also strong branching? */
    XP_DO_CTRL(get, dbl, XPRS_MIPABSCUTOFF, &old_upperbound);
    ++cleanup_level;
    XP_DO_CTRL(get, int, XPRS_AUTOPERTURB, &old_autoperturb);
    ++cleanup_level;

    XP_DO_CTRL(set, dbl, XPRS_MIPABSCUTOFF, upperbound);
    XP_DO_CTRL(set, int, XPRS_AUTOPERTURB, 0);

    ncols = CClp_ncols(lp);
    XP_CHK(ncols == 0 ? 1 : 0);

    XP_CHK_CALLOC(x, ncols, double);
    XP_CHK(XPRSgetsolution(lp->xp_prob, &status, x, 0, ncols - 1));
    XP_CHK(status == XPRS_LP_INFEAS ? status : 0);

    XP_CHK_CALLOC(candlist_twice, 2 * ncand, int);
    XP_CHK_CALLOC(boundstypes, 2 * ncand, char);
    XP_CHK_CALLOC(bounds, 2 * ncand, double);
    for (i = 0; i < ncand; i++) {
        cand                      = candidatelist[i];
        candlist_twice[2 * i + 0] = cand;
        candlist_twice[2 * i + 1] = cand;
        boundstypes[2 * i + 0]    = 'L';
        boundstypes[2 * i + 1]    = 'U';
        bounds[2 * i + 0]         = floor(x[cand]);
        bounds[2 * i + 1]         = ceil(x[cand]);
    }

    XP_CHK_CALLOC(objs, 2 * ncand, double);
    XP_CHK_CALLOC(statuses, 2 * ncand, int);
    XP_CHK(XPRSstrongbranch(lp->xp_prob, 2 * ncand, candlist_twice, boundstypes,
                            bounds, iterations, objs, statuses));

    for (i = 0; i < ncand; i++) {
        downpen[i] = XP_MIN(objs[i * 2 + 0], upperbound);
        uppen[i]   = XP_MIN(objs[i * 2 + 1], upperbound);
    }

CLEANUP:
    CC_IFFREE(x, double);
    CC_IFFREE(candlist_twice, int);
    CC_IFFREE(boundstypes, char);
    CC_IFFREE(bounds, double);
    CC_IFFREE(objs, double);
    CC_IFFREE(statuses, int);

    if (cleanup_level >= 2)
        XP_CUMUL_CTRL(set, int, XPRS_AUTOPERTURB, old_autoperturb);
    if (cleanup_level >= 1)
        XP_CUMUL_CTRL(set, dbl, XPRS_MIPABSCUTOFF, old_upperbound);

    return rval;
}
