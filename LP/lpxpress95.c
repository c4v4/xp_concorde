/*LD_PRELOAD=~/opt/xpressmp/lib/libxprs.so.43  gdb --args TSP/concorde
 * att48.tsp*/

#include "../INCLUDE/lp.h"
#include "../INCLUDE/machdefs.h"
#include "../INCLUDE/macrorus.h"
#include "../INCLUDE/util.h"
#include "../ext/xprs.h"
#include <math.h>

#undef CC_XPRESS_DISPLAY

#undef XP_DEBUG
#ifdef XP_DEBUG
#define XP_IF_DEBUG
#else
#define XP_IF_DEBUG if (0)
#endif

struct CClp {
  XPRSprob xp_prob;
};

struct CClp_info {
  int rcount;
  int ccount;
  int *rstat;
  int *cstat;
};

struct CClp_warmstart {
  struct CClp_info i;
};

static int CClp_env_count = 0;

/* Without macros, error handling results in probably more than 80% of the code.
 * Maybe it's quite heavy, but for sure also the repeated code has been reduced
 * a lot.*/

#define SOLVER_WARMSTART_NAME "CPL5"

#define _XP_STR(X) _XP__STR(X)
#define _XP__STR(X) #X
#define _XP_LINE_STR _XP_STR(__LINE__)

/* Current location in source file */
#define XP_LOC_STR __FILE__ ":" _XP_LINE_STR " > "

/* Check if EXPR return a non zero value, if it does go to CLEANUP */
#define XP_CHK(EXPR)                                                           \
  do {                                                                         \
    rval = (EXPR);                                                             \
    if (rval != 0) {                                                           \
      fprintf(stderr, XP_LOC_STR "'" #EXPR "' failed with status %d\n", rval); \
      goto CLEANUP;                                                            \
    }                                                                          \
  } while (0)

/* Work with Xpress controls, OP={get, set}, TYPE={int, dbl} */
#define XP_DO_CTRL(OP, TYPE, CTRL, VALUE)                                      \
  XP_CHK(XPRS##OP##TYPE##control(lp->xp_prob, CTRL, VALUE))

/* Get Xpress attribute */
#define XP_GET_ATTR(TYPE, CTRL, VALUE)                                         \
  XP_CHK(XPRSget##TYPE##attrib(lp->xp_prob, CTRL, VALUE))

/* During cleanup we cannot goto CLEANUP, so we cumulate errors */
#define XP_CUMUL_CHK(EXPR)                                                     \
  do {                                                                         \
    rval |= (EXPR);                                                            \
    if (rval != 0) {                                                           \
      fprintf(stderr,                                                          \
              XP_LOC_STR "'" #EXPR                                             \
                         "' failed with status %d, continuing anyway.\n",      \
              rval);                                                           \
    }                                                                          \
  } while (0)

/* During cleanup we cannot goto CLEANUP, so we cumulate errors */
#define XP_CUMUL_CTRL(OP, TYPE, CTRL, VALUE)                                   \
  XP_CUMUL_CHK(XPRS##OP##TYPE##control(lp->xp_prob, CTRL, VALUE))

/* Allocate, zero and check some space in memory */
#define XP_CHK_CALLOC(PTR, COUNT, TYPE)                                        \
  do {                                                                         \
    PTR = CC_SAFE_MALLOC(COUNT, TYPE);                                         \
    if (PTR == (TYPE *)NULL) {                                                 \
      fprintf(stderr,                                                          \
              XP_LOC_STR "Out of memory allocating " #PTR " = %d " #TYPE "\n", \
              COUNT);                                                          \
      rval = 1;                                                                \
      goto CLEANUP;                                                            \
    }                                                                          \
    memset(PTR, 0, COUNT * sizeof(TYPE));                                      \
  } while (0)

/* Private linkage helper functions: */

static int priv_check_lp_status(CClp *lp);
static int priv_print_stop_reason(CClp *lp);
static int priv_check_sol_status(CClp *lp, int *sol_status);
static int priv_primalopt(CClp *lp);
static int priv_dualopt(CClp *lp);
static int priv_baropt(CClp *lp);
static int priv_lpopt(CClp *lp, char pd, char pd_backup);
static void priv_messagecb(XPRSprob cbprob, void *cbdata, const char *msg,
                           int len, int msgtype);

static int priv_check_lp_status(CClp *lp) {
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

static int priv_print_stop_reason(CClp *lp) {
  int rval = 0;
  int stop_status = 0;

  fprintf(stderr, "Xpress optimization was stopped\n");
  XP_GET_ATTR(int, XPRS_STOPSTATUS, &stop_status);
  if (stop_status == XPRS_STOP_NONE)
    fprintf(stderr, "Xpress no interruption - the solve completed normally\n");
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

static int priv_check_sol_status(CClp *lp, int *sol_status) {
  int rval = 0;
  int status = 0;

  XP_GET_ATTR(int, XPRS_LPSTATUS, &status);

  if (status != XPRS_LP_OPTIMAL) {
    if (status == XPRS_LP_INFEAS) {
      fprintf(stderr, "Xpress LP: Infeasible\n");
    } else if (status == XPRS_LP_UNBOUNDED) {
      fprintf(stderr, "Xpress LP: Unbounded\n");
    } else if (status == XPRS_LP_UNSTARTED) {
      fprintf(stderr, "Xpress LP: Unstarted\n");
    } else if (status == XPRS_LP_CUTOFF) {
      fprintf(stderr, "Xpress LP: Objective worse than cutoff\n");
    } else if (status == XPRS_LP_UNFINISHED) {
      fprintf(stderr, "Xpress LP: Unfinished\n");
    } else if (status == XPRS_LP_CUTOFF_IN_DUAL) {
      fprintf(stderr, "Xpress LP: Cutoff in dual\n");
    } else if (status == XPRS_LP_UNSOLVED) {
      fprintf(stderr, "Xpress LP: Problem encountered numerical issues\n");
    } else if (status == XPRS_LP_NONCONVEX) {
      fprintf(stderr, "Xpress LP: Problem not convex\n");
    }
  }

  if (sol_status != (int *)NULL) {
    *sol_status = status;
  }

CLEANUP:
  return rval;
}

static int priv_lpopt(CClp *lp, char pd, char pd_backup) {
  int rval = 0;
  int solstatus;

#ifdef XP_DEBUG
  double obj;
  int ncols = CClp_ncols(lp);
  int nrows = CClp_nrows(lp);
#endif

  XP_CHK(XPRSlpoptimize(lp->xp_prob, &pd));
  XP_CHK(priv_check_lp_status(lp));
  XP_CHK(priv_check_sol_status(lp, &solstatus));
  if ((solstatus == XPRS_LP_INFEAS || solstatus == XPRS_LP_UNBOUNDED) &&
      pd_backup != '\0') {
    XP_CHK(priv_lpopt(lp, pd_backup, '\0'));
  }

#ifdef XP_DEBUG
  XP_CHK(XPRSgetdblattrib(lp->xp_prob, XPRS_LPOBJVAL, &obj));
  printf("CClp_opt: obj = %f (%d x %d)\n", obj, nrows, ncols);
#endif

CLEANUP:
  return rval;
}

static int priv_primalopt(CClp *lp) { return priv_lpopt(lp, 'p', '\0'); }

static int priv_dualopt(CClp *lp) { return priv_lpopt(lp, 'd', '\0'); }

static int priv_baropt(CClp *lp) { return priv_lpopt(lp, 'b', 'd'); }

static void priv_free_info(CClp_info *i) {
  if (i != (CClp_info *)NULL) {
    CC_IFFREE(i->cstat, int);
    CC_IFFREE(i->rstat, int);
  }
}

static int priv_init_info(CClp_info *i, int ccount, int rcount) {
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

static void priv_messagecb(XPRSprob cbprob, void *cbdata, const char *msg,
                           int len, int msgtype) {
  (void)cbprob; /* unused */
  (void)cbdata; /* unused */
  if (msgtype == 3 || msgtype == 4) {
    fprintf(stderr, "XPR> %*s\n", len, msg);
    fflush(stdout);
  }
#ifdef CC_XPRESS_DISPLAY
  if (msgtype == 1 || msgtype == 2) {
    fprintf(stdout, "XPR> %*s\n", len, msg);
    fflush(stdout);
  }
#endif
}

/* Public linkage exported functions: */

int CClp_init(CClp **lp_ptr) {
  XP_IF_DEBUG printf("CClp_init\n");
  int rval = 0;
  char message[512];
  CClp *lp = (CClp *)NULL;
  CClp_free(lp_ptr);

  XP_CHK_CALLOC(*lp_ptr, 1, CClp);
  lp = *lp_ptr;

  if (CClp_env_count++ == 0) {
    rval = XPRSinit("/home/cava/opt/xpressmp/bin/xpauth.xpr");
    if (rval != 0) {
      fprintf(stderr, "XPRSinit failed, return code %d\n", rval);
      XPRSgetlicerrmsg(message, sizeof(message));
      fprintf(stderr, "Licensing error: %s\n", message);
      goto CLEANUP;
    }
  }

CLEANUP:
  return rval;
}

int CClp_force_perturb(CClp *lp) {
  int rval = 0;

  XP_IF_DEBUG printf("CClp_force_perturb\n");
  XP_DO_CTRL(set, int, XPRS_AUTOPERTURB, 1);

CLEANUP:
  return 0;
}

int CClp_tune_small(CClp *lp) {
  int rval = 0;

  XP_DO_CTRL(set, int, XPRS_DUALGRADIENT, 0);
  XP_DO_CTRL(set, int, XPRS_PRICINGALG, 0);

CLEANUP:
  return rval;
}

int CClp_disable_presolve(CClp *lp) {
  int rval = 0;

  XP_DO_CTRL(set, int, XPRS_PRESOLVE, 0);

CLEANUP:
  return rval;
}

void CClp_free(CClp **lp_ptr) {
  CClp *lp = *lp_ptr;

  if (lp != (CClp *)NULL) {
    if (lp->xp_prob) {
      XPRSdestroyprob(lp->xp_prob);
    }
    lp->xp_prob = (XPRSprob)NULL;

    --CClp_env_count;
    if (CClp_env_count == 0) {
      XPRSfree();
    }

    CC_FREE(lp, CClp);
  }
}

void CClp_freelp(CClp **lp_ptr) {
  CClp *lp = *lp_ptr;

  if (lp != (CClp *)NULL) {
    if (lp->xp_prob) {
      XPRSdestroyprob(lp->xp_prob);
    }
    lp->xp_prob = (XPRSprob)NULL;
  }
}

int CClp_create(CClp *lp, const char *name) {
  int rval = 0;

  XP_CHK(XPRScreateprob(&lp->xp_prob));
  XP_CHK(XPRSsetprobname(lp->xp_prob, name));

  if (lp->xp_prob == NULL) {
    fprintf(stderr, "XPRScreateprob failed, return NULL object\n");
    rval = 2;
    goto CLEANUP;
  }

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

CLEANUP:
  return rval;
}

int CClp_loadlp(CClp *lp, const char *name, int ncols, int nrows, int objsense,
                double *obj, double *rhs, char *sense, int *matbeg, int *matcnt,
                int *matind, double *matval, double *lb, double *ub) {
  int rval = 0;
  double *rng = (double *)NULL;
  XP_IF_DEBUG printf("CClp_loadlp\n");

  XP_CHK(CClp_create(lp, name));
  XP_CHK(XPRSloadlp(lp->xp_prob, name, ncols, nrows, sense, rhs, rng, obj,
                    matbeg, matcnt, matind, matval, lb, ub));
  XP_CHK(XPRSchgobjsense(lp->xp_prob, objsense));
CLEANUP:
  return rval;
}

int CClp_new_row(CClp *lp, char sense, double rhs) {
  int *start = (int *)NULL, *colind = (int *)NULL;
  double *rng = (double *)NULL, *rowcoef = (double *)NULL;
  return XPRSaddrows(lp->xp_prob, 1, 0, &sense, &rhs, rng, start, colind,
                     rowcoef);
}

int CClp_change_sense(CClp *lp, int row, char sense) {
  return XPRSchgrowtype(lp->xp_prob, 1, &row, &sense);
}

int CClp_opt(CClp *lp, int method) {
  switch (method) {
  case CClp_METHOD_PRIMAL:
    return priv_primalopt(lp);
  case CClp_METHOD_DUAL:
    return priv_dualopt(lp);
  case CClp_METHOD_BARRIER:
    return priv_baropt(lp);
  }

  fprintf(stderr, "Nonexistent method in CClp_opt\n");
  return 1;
}

#define CC_MAX_REFACTORFREQ 150

int CClp_limited_dualopt(CClp *lp, int iterationlim, int *status,
                         double *objupperlim) {

  int rval = 0;
  int cleanup_level = 0;
  int solstatus;

  int old_autoperturb;
  int old_lpiterlimit;
  double old_abscutoff;
  int old_presolve;
  int old_invertfreq;

  XP_DO_CTRL(get, int, XPRS_AUTOPERTURB, &old_autoperturb);
  XP_DO_CTRL(get, int, XPRS_LPITERLIMIT, &old_lpiterlimit);
  XP_DO_CTRL(get, dbl, XPRS_MIPABSCUTOFF, &old_abscutoff);
  XP_DO_CTRL(get, int, XPRS_PRESOLVE, &old_presolve);
  XP_DO_CTRL(get, int, XPRS_INVERTFREQ, &old_invertfreq);

  XP_DO_CTRL(set, int, XPRS_AUTOPERTURB, 0);
  ++cleanup_level;
  XP_DO_CTRL(set, int, XPRS_LPITERLIMIT, iterationlim);
  ++cleanup_level;
  XP_DO_CTRL(set, dbl, XPRS_MIPABSCUTOFF, *objupperlim);
  ++cleanup_level;
  XP_DO_CTRL(set, int, XPRS_PRESOLVE, 0);
  ++cleanup_level;
  XP_DO_CTRL(set, int, XPRS_INVERTFREQ, iterationlim + 1);
  ++cleanup_level;

  XP_CHK(XPRSlpoptimize(lp->xp_prob, "d"));
  XP_CHK(priv_check_sol_status(lp, &solstatus));

  if (status == (int *)NULL)
    goto CLEANUP;

  if (solstatus == XPRS_LP_OPTIMAL) {
    *status = CClp_SUCCESS;
  } else if (solstatus == XPRS_LP_INFEAS) {
    *status = CClp_INFEASIBLE;
  } else if (solstatus == XPRS_LP_UNFINISHED) {
    *status = CClp_UNKNOWN;
  } else if (solstatus == XPRS_LP_UNBOUNDED) {
    *status = CClp_UNKNOWN;
  } else
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

  if (rval && status) {
    *status = CClp_FAILURE;
  }

  return rval;
}

int CClp_addrows(CClp *lp, int newrows, int newnz, double *rhs, char *sense,
                 int *rmatbeg, int *rmatind, double *rmatval) {

  double *rng = (double *)NULL;
  XP_IF_DEBUG printf("CClp_addrows: adding %d rows, %d nnz\n", newrows, newnz);
  return XPRSaddrows(lp->xp_prob, newrows, newnz, sense, rhs, rng, rmatbeg,
                     rmatind, rmatval);
}

int CClp_addcols(CClp *lp, int newcols, int newnz, double *obj, int *cmatbeg,
                 int *cmatind, double *cmatval, double *lb, double *ub) {

  XP_IF_DEBUG printf("CClp_addcols: adding %d cols, %d nnz\n", newcols, newnz);
  return XPRSaddcols(lp->xp_prob, newcols, newnz, obj, cmatbeg, cmatind,
                     cmatval, lb, ub);
}

int CClp_delete_row(CClp *lp, int i) {
  int rval = 0;
  int outlist;
  double *x = (double *)NULL;
  double objval;
  int npivots;
  int is_basic;

  /* Uncomment if you want to "help" the solver by not invalidating the basis */
  /*
    XP_CHK(XPRSgetbasisval(lp->xp_prob, i, 0, &is_basic, (int *)NULL));
    if (!is_basic) {
        XP_CHK(XPRSgetpivots(lp->xp_prob, i, &outlist, x, &objval, &npivots,
    1)); XP_CHK(XPRSpivot(lp->xp_prob, i, outlist));
    }
  */

  XP_CHK(XPRSdelrows(lp->xp_prob, 1, &i));

CLEANUP:
  return rval;
}

int CClp_delete_set_of_rows(CClp *lp, int *delstat) {
  int rval = 0;
  int *dellist = (int *)NULL;
  int delcnt = 0;
  int i, j;
  int rcnt = CClp_nrows(lp);

  for (i = 0; i < rcnt; i++) {
    if (delstat[i])
      delcnt++;
  }
  if (delcnt == 0) {
    fprintf(stderr, "delete_set_of_rows with no deleted rows\n");
    goto CLEANUP;
  }

  XP_CHK_CALLOC(dellist, delcnt, int);
  for (i = 0, j = 0; i < rcnt; i++) {
    if (delstat[i]) {
      dellist[j++] = i;
    }
  }

#ifdef XP_DEBUG
  printf("CClp_delete_set_of_rows: ");
  for (i = 0; i < delcnt; i++)
    printf("%d ", dellist[i]);
  printf("\n");
#endif

  XP_CHK(XPRSdelrows(lp->xp_prob, delcnt, dellist));

CLEANUP:
  CC_FREE(dellist, int);
  return rval;
}
int CClp_delete_column(CClp *lp, int i) {
  int rval = 0;
  /*
    double bd = 0.0;
    int is_basic;

    XP_CUMUL_CHK(XPRSchgbounds(lp->xp_prob, 1, &i, "B", &bd));
    XP_CUMUL_CHK(XPRSlpoptimize(lp->xp_prob, "d"));

    Here I should find another column to enter the basis replacing i
  */

  XP_CHK(XPRSdelcols(lp->xp_prob, 1, &i));

CLEANUP:
  return rval;
}

int CClp_delete_set_of_columns(CClp *lp, int *delstat) {
  int rval = 0;
  int *dellist = (int *)NULL;
  char *lu = (char *)NULL;
  double *bd = (double *)NULL;
  int delcnt = 0;
  int i;
  int j;
  int ccnt = CClp_ncols(lp);

  for (i = 0; i < ccnt; i++) {
    if (delstat[i])
      delcnt++;
  }

  if (delcnt == 0) {
    fprintf(stderr, XP_LOC_STR "delete_set_of_columns no columns to delete\n");
    return 0;
  }

  XP_CHK_CALLOC(dellist, delcnt, int);
  /* XP_CHK_CALLOC(lu, delcnt, char);
  XP_CHK_CALLOC(bd, delcnt, double); */
  for (i = 0, j = 0; i < ccnt; i++) {
    if (delstat[i]) {
      /*  lu[j] = 'B';
       bd[j] = 0.0; */
      dellist[j++] = i;
    }
  }

  /*
    XP_CUMUL_CHK(XPRSchgbounds(lp->xp_prob, delcnt, dellist, lu, bd));
    XP_CUMUL_CHK(XPRSlpoptimize(lp->xp_prob, "d"));
  */

  XP_CHK(XPRSdelcols(lp->xp_prob, delcnt, dellist));

CLEANUP:
  CC_FREE(dellist, int);
  CC_FREE(lu, char);
  CC_FREE(bd, double);

  return rval;
}

int CClp_setbnd(CClp *lp, int col, char lower_or_upper, double bnd) {
  XP_IF_DEBUG printf("CClp_setbnd: %d %c %f\n", col, lower_or_upper, bnd);
  return XPRSchgbounds(lp->xp_prob, 1, &col, &lower_or_upper, &bnd);
}

int CClp_get_warmstart(CClp *lp, CClp_warmstart **w) {
  int rval = 0;
  CClp_info *winfo = (CClp_info *)NULL;

  if (*w == (CClp_warmstart *)NULL)
    XP_CHK_CALLOC(*w, 1, CClp_warmstart);

  winfo = &(*w)->i;
  XP_CHK(priv_init_info(winfo, CClp_ncols(lp), CClp_nrows(lp)));
  XP_CHK(XPRSgetbasis(lp->xp_prob, winfo->rstat, winfo->cstat));

  return 0;

CLEANUP:

  CClp_free_warmstart(w);
  return rval;
}

int CClp_load_warmstart(CClp *lp, CClp_warmstart *w) {
  int rval = 0;
  XP_IF_DEBUG printf("CClp_load_warmstart\n");

  if (w->i.cstat != (int *)NULL && w->i.rstat != (int *)NULL) {
    XP_CHK(XPRSloadbasis(lp->xp_prob, w->i.rstat, w->i.cstat));
  } else {
    printf("WARNING: No basis in call to load_warmstart\n");
    fflush(stdout);
  }

CLEANUP:

  return rval;
}

int CClp_build_warmstart(CClp_warmstart **w, CClp_info *i) {
  XP_IF_DEBUG printf("CClp_build_warmstart\n");
  int rval = 0;
  CClp_info *winfo = (CClp_info *)NULL;
  int j;

  if (*w == (CClp_warmstart *)NULL)
    XP_CHK_CALLOC(*w, 1, CClp_warmstart);

  winfo = &(*w)->i;
  XP_CHK(priv_init_info(winfo, i->ccount, i->rcount));

  for (j = 0; j < i->ccount; j++) {
    winfo->cstat[j] = i->cstat[j];
  }
  for (j = 0; j < i->rcount; j++) {
    winfo->rstat[j] = i->rstat[j];
  }

  return 0;

CLEANUP:
  CClp_free_warmstart(w);
  return rval;
}

void CClp_free_warmstart(CClp_warmstart **w) {
  if (w == (CClp_warmstart **)NULL || *w == (CClp_warmstart *)NULL)
    return;
  priv_free_info(&(*w)->i);
  CC_FREE(*w, CClp_warmstart);
}

int CClp_sread_warmstart(CC_SFILE *f, CClp_warmstart **w) {
  int rval = 0;
  char name[5];
  int i;
  int ccount;
  int rcount;
  int has_dnorms;
  CClp_info *winfo = (CClp_info *)NULL;

  for (i = 0; i < 4; i++) {
    XP_CHK(CCutil_sread_char(f, name + i));
  }
  name[4] = '\0';

  if (strncmp(name, SOLVER_WARMSTART_NAME, 4)) {
    fprintf(stderr, "warmstart for another solver (%s) ignored\n", name);
    return 0;
  }

  XP_CHK(CCutil_sread_int(f, &ccount));
  XP_CHK(CCutil_sread_int(f, &rcount));

  if (*w == (CClp_warmstart *)NULL)
    XP_CHK_CALLOC(*w, 1, CClp_warmstart);

  winfo = &(*w)->i;
  XP_CHK(priv_init_info(winfo, ccount, rcount));

  for (i = 0; i < ccount; i++) {
    XP_CHK(CCutil_sread_bits(f, winfo->cstat + i, 2));
  }
  for (i = 0; i < rcount; i++) {
    XP_CHK(CCutil_sread_bits(f, winfo->rstat + i, 1));
  }

  XP_CHK(CCutil_sread_int(f, &has_dnorms));
  /* Still have to find a way to get and set norms with Xpress */

  return 0;

CLEANUP:

  CClp_free_warmstart(w);
  return rval;
}

int CClp_swrite_warmstart(CC_SFILE *f, CClp_warmstart *w) {
  int rval = 0;
  int i;
  const char *name = SOLVER_WARMSTART_NAME;
  CClp_info *info = &(w->i);

  for (i = 0; i < 4; i++) {
    XP_CHK(CCutil_swrite_char(f, name[i]));
  }

  XP_CHK(CCutil_swrite_int(f, info->ccount));
  XP_CHK(CCutil_swrite_int(f, info->rcount));

  for (i = 0; i < info->ccount; i++) {
    XP_CHK(CCutil_swrite_bits(f, info->cstat[i], 2));
  }

  for (i = 0; i < info->rcount; i++) {
    XP_CHK(CCutil_swrite_bits(f, info->rstat[i], 1));
  }
  XP_CHK(CCutil_swrite_int(f, 0));

CLEANUP:
  return rval;
}

void CClp_free_info(CClp_info **i) {
  if (i == (CClp_info **)NULL || *i == (CClp_info *)NULL)
    return;
  priv_free_info(*i);
  CC_FREE(*i, CClp_info);
}

int CClp_get_info(CClp *lp, CClp_info **i) {
  int rval = 0;

  if (*i == (CClp_info *)NULL)
    XP_CHK_CALLOC(*i, 1, CClp_info);
  XP_CHK(priv_init_info(*i, CClp_ncols(lp), CClp_nrows(lp)));
  XP_CHK(XPRSgetbasis(lp->xp_prob, (*i)->rstat, (*i)->cstat));

  XP_IF_DEBUG printf(
      "CClp_get_info: %d cols, %d rows, c %d %d %d..., r %d %d %d...\n",
      (*i)->ccount, (*i)->rcount, (*i)->cstat[0], (*i)->cstat[1],
      (*i)->cstat[2], (*i)->rstat[0], (*i)->rstat[1], (*i)->rstat[2]);

  return 0;

CLEANUP:

  CClp_free_info(i);
  return rval;
}

int CClp_create_info(CClp_info **i, int rcount, int ccount) {
  XP_IF_DEBUG printf("CClp_create_info\n");
  int rval = 0;
  int j;

  if (*i == (CClp_info *)NULL)
    XP_CHK_CALLOC(*i, 1, CClp_info);
  XP_CHK(priv_init_info(*i, ccount, rcount));

  for (j = 0; j < ccount; j++) {
    (*i)->cstat[j] = 0;
  }
  for (j = 0; j < rcount; j++) {
    (*i)->rstat[j] = 0;
  }

  return 0;

CLEANUP:
  CClp_free_info(i);
  return rval;
}

int CClp_is_col_active(CClp_info *i, int c) {
  if (c < 0 || c >= i->ccount)
    return 0;
  return i->cstat[c] == 1 || i->cstat[c] == 2;
}

int CClp_is_row_active(CClp_info *i, int r) {
  if (r < 0 || r >= i->rcount)
    return 0;
  return i->rstat[r] == 0 || i->rstat[r] == 2;
}

void CClp_set_col_active(CClp_info *i, int c) {
  if (c >= 0 && c < i->ccount)
    i->cstat[c] = 1; /*variable is basic*/
}

void CClp_set_col_inactive(CClp_info *i, int c) {
  if (c >= 0 && c < i->ccount)
    i->cstat[c] = 0; /*variable at lower bound*/
}

void CClp_set_col_upper(CClp_info *i, int c) {
  if (c >= 0 && c < i->ccount)
    i->cstat[c] = 2; /*variable at upper bound*/
}

void CClp_set_row_active(CClp_info *i, int r) {
  if (r >= 0 && r < i->rcount)
    i->rstat[r] = 2;
  /*associated slack/surplus/artificial-var is nonbasic at value 0.0 */
}

void CClp_set_row_inactive(CClp_info *i, int r) {
  if (r >= 0 && r < i->rcount)
    i->rstat[r] = 1;
  /*associated slack/surplus/artificial-var is basic*/
}

int CClp_x(CClp *lp, double *x) {
  int rval = 0;
  int status;

  int ncols = CClp_ncols(lp);
  if (ncols == 0) {
    fprintf(stderr, "No columns in LP\n");
    goto CLEANUP;
  }

  XP_CHK(XPRSgetsolution(lp->xp_prob, &status, x, 0, ncols - 1));
  if (status == XPRS_LP_INFEAS) {
    fprintf(stderr, "Xpress solution status is %d\n", status);
    rval = 1;
    goto CLEANUP;
  }
  XP_IF_DEBUG printf("CClp_x: %.2f %.2f %.2f %.2f...\n", x[0], x[1], x[2],
                     x[3]);

CLEANUP:
  return rval;
}

int CClp_rc(CClp *lp, double *rc) {
  int rval = 0;
  int status;

  int ncols = CClp_ncols(lp);
  if (ncols == 0) {
    fprintf(stderr, "No columns in LP\n");
    goto CLEANUP;
  }

  XP_CHK(XPRSgetredcosts(lp->xp_prob, &status, rc, 0, ncols - 1));
  XP_IF_DEBUG printf("CClp_rc: %.2f %.2f %.2f %.2f...\n", rc[0], rc[1], rc[2],
                     rc[3]);
CLEANUP:
  return rval;
}

int CClp_pi(CClp *lp, double *pi) {
  int rval = 0;
  int nrows;
  int method;
  int status;

  XP_CHK(XPRSgetintcontrol(lp->xp_prob, XPRS_LPFLAGS, &method));
  XP_CHK(XPRSgetintattrib(lp->xp_prob, XPRS_LPSTATUS, &status));

  if (method == XPRS_LPFLAGS_DUAL && status == XPRS_LP_INFEAS) {
    int has_ray = 0;
    XP_CHK(XPRSgetdualray(lp->xp_prob, pi, &has_ray));
    if (!has_ray) {
      fprintf(stderr, "XPRSgetdualray failed\n");
      rval = 1;
    }
    XP_IF_DEBUG printf("CClp_pi: infeas %.2f %.2f %.2f %.2f...\n", pi[0], pi[1],
                       pi[2], pi[3]);

  } else {
    nrows = CClp_nrows(lp);
    if (nrows == 0) {
      fprintf(stderr, "No rows in LP\n");
      rval = 1;
      goto CLEANUP;
    }

    XP_CHK(XPRSgetduals(lp->xp_prob, &status, pi, 0, nrows - 1));
    XP_IF_DEBUG printf("CClp_pi: feas %.2f %.2f %.2f %.2f...\n", pi[0], pi[1],
                       pi[2], pi[3]);
  }

CLEANUP:
  return rval;
}

int CClp_objval(CClp *lp, double *obj) {
  int rval = XPRSgetdblattrib(lp->xp_prob, XPRS_LPOBJVAL, obj);
  XP_IF_DEBUG printf("CClp_objval: %f\n", *obj);
  return rval;
}

/* The following functions are not expected to fail */
int CClp_nrows(CClp *lp) {
  int res;
  if (XPRSgetintattrib(lp->xp_prob, XPRS_ORIGINALROWS, &res)) {
    fprintf(stderr, "XPRSgetintattrib XPRS_ORIGINALROWS failed\n");
    return 0;
  }
  return res;
}

int CClp_ncols(CClp *lp) {
  int res;
  if (XPRSgetintattrib(lp->xp_prob, XPRS_ORIGINALCOLS, &res)) {
    fprintf(stderr, "XPRSgetintattrib XPRS_ORIGINALCOLS failed\n");
    return 0;
  }
  return res;
}

int CClp_nnonzeros(CClp *lp) {
  int is_presolved;
  int res;

  int rval = XPRSgetintattrib(lp->xp_prob, XPRS_PRESOLVESTATE, &is_presolved);
  if (rval) {
    fprintf(stderr, "XPRSgetintattrib XPRS_PRESOLVESTATE failed\n");
    return 0;
  }

  if (is_presolved) {
    fprintf(stderr, "XPRS_ELEMS: returned non-zeros of presolved model\n");
  }

  rval = XPRSgetintattrib(lp->xp_prob, XPRS_ELEMS, &res);
  if (rval) {
    fprintf(stderr, "XPRSgetintattrib XPRS_ELEMS failed\n");
    return 0;
  }

  return res;
}

int CClp_status(CClp *lp, int *status) {
  int rval = 0;
  int solmethod, solstat;

  XP_DO_CTRL(get, int, XPRS_LPFLAGS, &solmethod);
  if ((solmethod & (XPRS_LPFLAGS_DUAL | XPRS_LPFLAGS_PRIMAL)) == 0) {
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

int CClp_getweight(CClp *lp, int nrows, int *rmatbeg, int *rmatind,
                   double *rmatval, double *weight) {
  XP_IF_DEBUG printf("CClp_getweight ...\n");
  fflush(stdout);

  if (lp || nrows || rmatbeg || rmatind || rmatval || weight) {
    return 1;
  } else {
    return 0;
  }
}

int CClp_dump_lp(CClp *lp, const char *fname) {
  return XPRSwriteprob(lp->xp_prob, fname, "lp");
}

#define OURXPRESSZERO (1.0E-10)
#define OURXPRESS_INTTOL (0.0001)

int CClp_getgoodlist(CClp *lp, int *goodlist, int *goodlen_p, double *downpen,
                     double *uppen) {
  int rval = 0;
  int ncols, i, j;
  int status;
  int *cstat = (int *)NULL;
  double *x = (double *)NULL;

  /* Call XPRSlpoptimize and verify optimality */

  XP_CHK(XPRSlpoptimize(lp->xp_prob, "d"));
  ncols = CClp_ncols(lp);
  if (ncols == 0) {
    fprintf(stderr, "No columns in LP\n");
    rval = 1;
    goto CLEANUP;
  }

  XP_CHK_CALLOC(x, ncols, double);
  XP_CHK(XPRSgetsolution(lp->xp_prob, &status, x, 0, ncols - 1));
  if (status == XPRS_LP_INFEAS) {
    fprintf(stderr, "Xpress solution status is %d\n", status);
    rval = 1;
    goto CLEANUP;
  }

  XP_CHK_CALLOC(cstat, ncols, int);
  XP_CHK(XPRSgetbasis(lp->xp_prob, (int *)NULL, cstat));

  *goodlen_p = 0;
  for (i = 0, j = 0; i < ncols; i++) {
    if (cstat[i] == 1 && x[i] >= OURXPRESS_INTTOL &&
        x[i] <= 1.0 - OURXPRESS_INTTOL) {
      goodlist[j] = i;
      downpen[j] = x[i];
      uppen[j] = 1.0 - x[i];
      j++;
    }
  }

  *goodlen_p = j;

  /* TODO(cava): compute Driebeek penalties and do the same thing that is done
   * in lpcplex8.c (can be efficiently computed by ratio test) s*/

CLEANUP:

  CC_IFFREE(cstat, int);
  CC_IFFREE(x, double);
  return rval;
}

int CClp_strongbranch(CClp *lp, int *candidatelist, int ncand, double *downpen,
                      double *uppen, int iterations, double upperbound) {
  int rval = 0;
  double old_upperbound;
  int old_autoperturb;
  int cleanup_level = 0;
  int i;
  int status;
  int ncols;
  double *x = (double *)NULL;
  int *candlist_twice = (int *)NULL;
  char *boundstypes = (char *)NULL;
  double *bounds = (double *)NULL;
  double *objs = (double *)NULL;
  int *statuses = (int *)NULL;

  XP_IF_DEBUG printf("CClp_strongbranch\n");

  /* Does it limit also strong branching? */
  XP_DO_CTRL(get, dbl, XPRS_MIPABSCUTOFF, &old_upperbound);
  ++cleanup_level;
  XP_DO_CTRL(get, int, XPRS_AUTOPERTURB, &old_autoperturb);
  ++cleanup_level;

  XP_DO_CTRL(set, dbl, XPRS_MIPABSCUTOFF, upperbound);
  XP_DO_CTRL(set, int, XPRS_AUTOPERTURB, 0);

  ncols = CClp_ncols(lp);
  if (ncols == 0) {
    fprintf(stderr, "No columns in LP\n");
    return 1;
  }

  XP_CHK_CALLOC(x, ncols, double);
  XP_CHK(XPRSgetsolution(lp->xp_prob, &status, x, 0, ncols - 1));
  if (status == XPRS_LP_INFEAS) {
    fprintf(stderr, "Xpress solution status is %d\n", status);
    rval = 1;
    goto CLEANUP;
  }

  XP_CHK_CALLOC(candlist_twice, 2 * ncand, int);
  XP_CHK_CALLOC(boundstypes, 2 * ncand, char);
  XP_CHK_CALLOC(bounds, 2 * ncand, double);
  for (i = 0; i < ncand; i++) {
    candlist_twice[2 * i + 0] = candidatelist[i];
    candlist_twice[2 * i + 1] = candidatelist[i];
    boundstypes[2 * i + 0] = 'L';
    boundstypes[2 * i + 1] = 'U';
    bounds[2 * i + 0] = floor(x[candidatelist[i]]);
    bounds[2 * i + 1] = ceil(x[candidatelist[i]]);
  }

  XP_CHK_CALLOC(objs, 2 * ncand, double);
  XP_CHK_CALLOC(statuses, 2 * ncand, int);
  XP_CHK(XPRSstrongbranch(lp->xp_prob, 2 * ncand, candlist_twice, boundstypes,
                          bounds, iterations, objs, statuses));

  for (i = 0; i < ncand; i++) {
    downpen[i] = objs[i * 2] < upperbound ? objs[i * 2] : upperbound;
    uppen[i] = objs[i * 2 + 1] < upperbound ? objs[i * 2 + 1] : upperbound;
  }

CLEANUP:
  CC_IFFREE(candlist_twice, int);
  CC_IFFREE(boundstypes, char);
  CC_IFFREE(bounds, double);
  CC_IFFREE(objs, double);
  CC_IFFREE(x, double);

  if (cleanup_level >= 2)
    XP_CUMUL_CTRL(set, int, XPRS_AUTOPERTURB, old_autoperturb);
  if (cleanup_level >= 1)
    XP_CUMUL_CTRL(set, dbl, XPRS_MIPABSCUTOFF, old_upperbound);

  return rval;
}
