#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME:
Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP _kentr_countKmers(SEXP);
extern SEXP _kentr_find_mismatch(SEXP, SEXP);
extern SEXP _kentr_get_hamming(SEXP, SEXP);
extern SEXP _kentr_getKmers(SEXP, SEXP);
extern SEXP _kentr_getSeq(SEXP, SEXP);
extern SEXP _kentr_read_bam(SEXP);
extern SEXP _kentr_revComp(SEXP);

static const R_CallMethodDef CallEntries[] = {
  {"_kentr_countKmers",    (DL_FUNC) &_kentr_countKmers,    1},
  {"_kentr_find_mismatch", (DL_FUNC) &_kentr_find_mismatch, 2},
  {"_kentr_get_hamming",   (DL_FUNC) &_kentr_get_hamming,   2},
  {"_kentr_getKmers",      (DL_FUNC) &_kentr_getKmers,      2},
  {"_kentr_getSeq",        (DL_FUNC) &_kentr_getSeq,        2},
  {"_kentr_read_bam",      (DL_FUNC) &_kentr_read_bam,      1},
  {"_kentr_revComp",       (DL_FUNC) &_kentr_revComp,       1},
  {NULL, NULL, 0}
};

void R_init_kentr(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
