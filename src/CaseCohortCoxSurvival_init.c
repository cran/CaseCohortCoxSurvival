#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME:
   Check these declarations against the C/Fortran source code.
*/


/* .C calls */
extern void C_getS0t(void *, void *, void *, void *, void *, void *);
extern void C_getS1t(void *, void *, void *, void *, void *, void *, void *, void *);
extern void C_getSumAAwgt(void *, void *, void *, void *, void *);
extern void C_getdNtColSums(void *, void *, void *, void *, void *);
extern void C_getdNtWgtColSums(void *, void *, void *, void *, void *, void *);
extern void C_get_drond_U_eta(void *, void *, void *, void *, void *, void *, void *, void *,
                              void *, void *, void *, void *, void *, void *, void *, void *);
extern void C_getBetaScore(void *, void *, void *, void *, void *, void *, void *,
                           void *, void *, void *, void *, void *);
extern void C_infl_lambda0_tau12(void *, void *, void *, void *, void *, void *, void *,
                                 void *, void *, void *, void *, void *, void *, void *,
                                 void *, void *, void *, void *, void *, void *, void *);
extern void C_infl_lambda0_tau12_noCalib(void *, void *, void *, void *, void *, void *, void *,
                                         void *, void *, void *, void *, void *, void *, void *);
extern void C_getPhase2Var(void *, void *, void *, void *, void *, void *, void *,
                           void *, void *, void *);
extern void C_getS0GammaCasetimes(void *, void *, void *, void *, void *, void *, void *, void *);
extern void  C_infl2_lambda0t(void *, void *, void *, void *, void *, void *, void *,
                              void *, void *, void *, void *, void *, void *, void *);
extern void  C_infl3_lambda0t(void *, void *, void *, void *, void *, void *, void *,
                              void *, void *, void *, void *, void *, void *, void *, void *);
extern void  C_infl2_lambda0t_noEst(void *, void *, void *, void *, void *, void *, void *);
extern void  C_infl3_lambda0t_noEst(void *, void *, void *, void *, void *, void *, void *,
                                    void *, void *, void *, void *, void *);
extern void  C_phase23VarEstF(void *, void *, void *, void *, void *, void *, void *, void *,
                              void *, void *, void *);
extern void  C_phase23VarEstT(void *, void *, void *, void *, void *, void *, void *, void *,
                              void *, void *, void *, void *, void *);
extern void  C_getRiskMatCol(void *, void *, void *, void *, void *);


static const R_CMethodDef CEntries[] = {
    {"C_getS0t",                     (DL_FUNC) &C_getS0t,                      6},
    {"C_getS1t",                     (DL_FUNC) &C_getS1t,                      8},
    {"C_getSumAAwgt",                (DL_FUNC) &C_getSumAAwgt,                 5},
    {"C_getdNtColSums",              (DL_FUNC) &C_getdNtColSums,               5},
    {"C_getdNtWgtColSums",           (DL_FUNC) &C_getdNtWgtColSums,            6},
    {"C_get_drond_U_eta",            (DL_FUNC) &C_get_drond_U_eta,            16},
    {"C_getBetaScore",               (DL_FUNC) &C_getBetaScore,               12},
    {"C_infl_lambda0_tau12",         (DL_FUNC) &C_infl_lambda0_tau12,         21},
    {"C_infl_lambda0_tau12_noCalib", (DL_FUNC) &C_infl_lambda0_tau12_noCalib, 14},
    {"C_getPhase2Var",               (DL_FUNC) &C_getPhase2Var,               10},
    {"C_getS0GammaCasetimes",        (DL_FUNC) &C_getS0GammaCasetimes,         8},
    {"C_infl2_lambda0t",             (DL_FUNC) &C_infl2_lambda0t,             14},
    {"C_infl3_lambda0t",             (DL_FUNC) &C_infl3_lambda0t,             15},
    {"C_infl2_lambda0t_noEst",       (DL_FUNC) &C_infl2_lambda0t_noEst,        7},
    {"C_infl3_lambda0t_noEst",       (DL_FUNC) &C_infl3_lambda0t_noEst,       12},
    {"C_phase23VarEstF",             (DL_FUNC) &C_phase23VarEstF,             11},
    {"C_phase23VarEstT",             (DL_FUNC) &C_phase23VarEstT,             13},
    {"C_getRiskMatCol",              (DL_FUNC) &C_getRiskMatCol,               5},
    {NULL, NULL, 0}
};

void R_init_CaseCohortCoxSurvival(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
