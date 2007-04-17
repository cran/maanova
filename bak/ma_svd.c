#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Lapack.h>

//#include "ma_svd.h"
//#include "Defn.h"
#include "config.h"

//static int initialized = 0;

static void ma_Init(void)
{
    int res = R_moduleCdynload("lapack", 1, 1);
//    initialized = -1;
    if(!res) return;
/*    if(!ptr->svd)
//    	error("lapack routines cannot be accessed in module");
	error(_("lapack routines cannot be accessed in module"));
    initialized = 1;*/
    return;
}
/*
SEXP La_svd(SEXP jobu, SEXP jobv, SEXP x, SEXP s, SEXP u, SEXP v, SEXP method)
{
    if(!initialized) La_Init();
    if(initialized > 0)
	return ma_La_svd(jobu, jobv, x, s, u, v, method);
    else {
		error("lapack rountines cannot be loaded");
	//error(_("lapack routines cannot be loaded"));
	return R_NilValue;
    }
}*/

SEXP ma_La_svd(SEXP jobu, SEXP jobv, SEXP x, SEXP s, SEXP u, SEXP v,
		      SEXP method)
{
	ma_Init();
    int *xdims, n, p, lwork, info = 0;
    double *work, *xvals, tmp;
    SEXP val, nm;
    char *meth;

    if (!(isString(jobu) && isString(jobv)))
	error("'jobu' and 'jobv' must be character strings");
    if (!isString(method))
	error("'method' must be a character string");
    meth = CHAR(STRING_ELT(method, 0));
/*#ifndef IEEE_754
    if (strcmp(meth, "dgesdd") == 0)
	error("method = \"dgesdd\" requires IEEE 754 arithmetic");
#endif*/
    xdims = INTEGER(coerceVector(getAttrib(x, R_DimSymbol), INTSXP));
    n = xdims[0]; p = xdims[1];
    xvals = (double *) R_alloc(n * p, sizeof(double));
    /* work on a copy of x */
    Memcpy(xvals, REAL(x), (size_t) (n * p));

    if(strcmp(meth, "dgesdd")) {
	/* ask for optimal size of work array */
	lwork = -1;
	F77_CALL(dgesvd)(CHAR(STRING_ELT(jobu, 0)), CHAR(STRING_ELT(jobv, 0)),
			 &n, &p, xvals, &n, REAL(s),
			 REAL(u), INTEGER(getAttrib(u, R_DimSymbol)),
			 REAL(v), INTEGER(getAttrib(v, R_DimSymbol)),
			 &tmp, &lwork, &info);
	if (info != 0)
		error("error code from Lapack routine dgesvd");
	    //error(_("error code %d from Lapack routine '%s'"), info, "dgesvd");
	lwork = (int) tmp;

	work = (double *) R_alloc(lwork, sizeof(double));
	F77_CALL(dgesvd)(CHAR(STRING_ELT(jobu, 0)), CHAR(STRING_ELT(jobv, 0)),
			 &n, &p, xvals, &n, REAL(s),
			 REAL(u), INTEGER(getAttrib(u, R_DimSymbol)),
			 REAL(v), INTEGER(getAttrib(v, R_DimSymbol)),
			 work, &lwork, &info);
	if (info != 0)
		error("error code from Lapack routine dgesvd");
	    //error(_("error code %d from Lapack routine '%s'"), info, "dgesvd");
    }else {
	int ldu = INTEGER(getAttrib(u, R_DimSymbol))[0],
	    ldvt = INTEGER(getAttrib(v, R_DimSymbol))[0];
	int *iwork= (int *) R_alloc(8*(n<p ? n : p), sizeof(int));

	/* ask for optimal size of work array */
	lwork = -1;
	F77_CALL(dgesdd)(CHAR(STRING_ELT(jobu, 0)),
			 &n, &p, xvals, &n, REAL(s),
			 REAL(u), &ldu,
			 REAL(v), &ldvt,
			 &tmp, &lwork, iwork, &info);
	if (info != 0)
		error("error code from Lapack routine dgesdd");
	    //error(_("error code %d from Lapack routine '%s'"), info, "dgesdd");
	lwork = (int) tmp;
	work = (double *) R_alloc(lwork, sizeof(double));
	F77_CALL(dgesdd)(CHAR(STRING_ELT(jobu, 0)),
			 &n, &p, xvals, &n, REAL(s),
			 REAL(u), &ldu,
			 REAL(v), &ldvt,
			 work, &lwork, iwork, &info);
	if (info != 0)
		error("error code from Lapack routine dgesdd");
	    //error(_("error code %d from Lapack routine '%s'"), info, "dgesdd");
    }

    val = PROTECT(allocVector(VECSXP, 3));
    nm = PROTECT(allocVector(STRSXP, 3));
    SET_STRING_ELT(nm, 0, mkChar("d"));
    SET_STRING_ELT(nm, 1, mkChar("u"));
    SET_STRING_ELT(nm, 2, mkChar("vt"));
    setAttrib(val, R_NamesSymbol, nm);
    SET_VECTOR_ELT(val, 0, s);
    SET_VECTOR_ELT(val, 1, u);
    SET_VECTOR_ELT(val, 2, v);
    UNPROTECT(2);
    return val;
}
