#include <math.h>
#include <stdlib.h>
#include <ctype.h>
#include <R.h>
#include <R_ext/Memory.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

#define CHECK_MEM(obj) if (obj == NULL) {Rprintf("ERROR: allocating memory \n"); error("1");}
#define MISSTIMEIND -2


/*
static double print_dVec(v, n)
double *v;
int n;
{
  int i;
  for (i=0; i<n; i++) Rprintf("%g ", v[i]);
  Rprintf("\n");
}

static double print_dMat(m, n1, n2)
double **m;
int n1, n2;
{
  int i, j;
  for (i=0; i<n1; i++) {
    for (j=0; j<n2; j++) Rprintf("%g ", m[i][j]);
    Rprintf("\n");
  }
}
*/

static double * dVec_alloc(int n, int initFlag, double initVal)
{
  int i;
  double *ret, *p;

  if (n < 1) error("n < 1 in dVec_alloc");
  ret = (double *) R_Calloc(n, double);
  CHECK_MEM(ret);
  if (initFlag) {
    for (i=0, p=ret; i<n; i++, p++) *p = initVal;
  }

  return(ret);

} 

static double ** dMat_alloc(int nr, int nc, int initFlag, double initVal)
{
  int i, j;
  double **ret, *p;

  ret = (double **) R_Calloc(nr, double *);
  CHECK_MEM(ret);
  for (i=0; i<nr; i++) {
    ret[i] = (double *) R_Calloc(nc, double);
    CHECK_MEM(ret[i]);
    if (initFlag) {
      p = ret[i];
      for (j=0; j<nc; j++) p[j] = initVal; 
    }
  }
  return(ret);

} 

static double *** dArr_alloc(int n1, int n2, int n3)
{
  int i, j;
  double ***ret;

  ret = (double ***) R_Calloc(n1, double **);
  CHECK_MEM(ret);
  for (i=0; i<n1; i++) {
    ret[i] = (double **) R_Calloc(n2, double *);
    CHECK_MEM(ret[i]);
    for (j=0; j<n2; j++) {
      ret[i][j] = (double *) R_Calloc(n3, double);
      CHECK_MEM(ret[i][j]);
    }
  }
  return(ret);

} 

static void dMat_free(double **x, int nr)
{
  int i;
  
  for (i=0; i<nr; i++) R_Free(x[i]);
  R_Free(x);

}

static void dArr_free(double ***x, int n1, int n2)
{
  int i, j;
  
  for (i=0; i<n1; i++) {
    for (j=0; j<n2; j++) R_Free(x[i][j]); 
    R_Free(x[i]);
  }
  R_Free(x);

}


static void putMatInVecByRow(double **mat, int nr, int nc, double *ret)
{
  int i, j, ind;
  
  ind = 0;
  for (i=0; i<nr; i++) {
    for (j=0; j<nc; j++) ret[ind++] = mat[i][j];
  }

}

static double dotprod(double *x1, double *x2, int n)
{
  double ret=0.0;
  int i;
  for (i=0; i<n; i++) ret += x1[i]*x2[i];
  return(ret);
}

static double dotprodV1plusV2(double *x, double *v1, double *v2, int n)
{
  double ret=0.0;
  int i;
  for (i=0; i<n; i++) ret += x[i]*(v1[i] + v2[i]);
  return(ret);
}

static void vecprod(double *x1, double *x2, int n, double *ret)
{
  int i;
  for (i=0; i<n; i++) ret[i] = x1[i]*x2[i];
}

static double vecsum(double *x, int n)
{
  double ret=0.0;
  int i;
  for (i=0; i<n; i++) ret += x[i];
  return(ret);
}

static double vecsumprod(double *x1, double *x2, int n)
{
  double ret=0.0;
  int i;
  for (i=0; i<n; i++) ret += x1[i]*x2[i];
  return(ret);
}

static void tcrossprod(double *x1, double *x2, int n1, int n2, double **ret)
{
  /* ret has dim n1 x n2 */
  int i, j;
  double x1i;

  for (i=0; i<n1; i++) {
    x1i = x1[i];
    for (j=0; j<n2; j++) ret[i][j] = x1i*x2[j];
  }  
}

static void getColFromMatByRow(double *matByRow, int nr, int nc, int col, double *ret)
{
  /* Ret has length nr */
  int i, ind=col;
  
  for (i=0; i<nr; i++) {
    ret[i] = matByRow[ind];
    ind += nc;
  }
}

/*
static void getMatMultRow(mat1row, mat2ByCol, n1, nc2, ret)
double *mat1row, *mat2ByCol, *ret;
int n1, nc2;
{
   ret has length nc2 
  int i;
  double *p2;

  p2 = mat2ByCol;
  for (i=0; i<nc2; i++) {
    ret[i] = dotprod(mat1row, p2, n1);
    p2 += n1;
  }

}
*/

static void colSumsColMatTimesVec(double *colMat, int nr, int nc, double *vec, double *ret)
{
  /* colMat is column major 
     ret must have length nc */

  int i;
  double *pmat;

  pmat = colMat;
  for (i=0; i<nc; i++) {
    ret[i] = vecsumprod(pmat, vec, nr);
    pmat += nr;
  }
}

static void colSumsRowMatTimesVec(double *rowMat, int nr, int nc, double *vec, double *tmpvec, double *ret)
{
  /* rowMat is row major 
     ret must have length nc
     tmpvec has length nr 
  */

  int i;

  for (i=0; i<nc; i++) {
    getColFromMatByRow(rowMat, nr, nc, i, tmpvec);
    ret[i] = vecsumprod(tmpvec, vec, nr);
  }
}

static void initVec(double *vec, int n, double val)
{
  /* vec is overwritten */

  int i;
  for (i=0; i<n; i++) vec[i] = val;
}

static void initMat(double **mat, int nr, int nc, double val)
{
  /* Mat is overwritten */

  int i, j;
  for (i=0; i<nr; i++) {
    for (j=0; j<nc; j++) mat[i][j] = val;
  }
}

static void applySweep3sum(double ***arr, int n1, int n2, int n3, double *vec, double **ret)
{
  /* ret must be n1 x n2
     vec must have length n3
     arr has dimension n1 x n2 x n3
  */
  int i, j, k;
  double val;

  initMat(ret, n1, n2, 0.0);

  for (k=0; k<n3; k++) {
    val = vec[k];
    for (i=0; i<n1; i++) {
      for (j=0; j<n2; j++) ret[i][j] += val*arr[i][j][k];
    } 
  }

}

static void getRiskMatCol2(int col, int *timeInd, int *time0Ind, int n, double *ret)
{
  /* time0Ind indicator for start time */
  int i, tind, t0ind;
  for (i=0; i<n; i++) {
    tind  = timeInd[i];
    t0ind = time0Ind[i];
    if ((t0ind < col) && (col <= tind)) {
      ret[i] = 1.0;
    } else {
      ret[i] = 0.0;
    }
  }
} 

static void getRiskMatCol(int col, int *timeInd, int *time0Ind, int n, double *ret)
{
  int i, tind;
  if (time0Ind[0] < MISSTIMEIND) {
    /* 1 time var */
    for (i=0; i<n; i++) {
      tind = timeInd[i];
      if ((tind >= 0) && (col <= tind)) {
        ret[i] = 1.0;
      } else {
        ret[i] = 0.0;
      }
    }
  } else {
    getRiskMatCol2(col, timeInd, time0Ind, n, ret);
  }
} 

void C_getRiskMatCol(int *col, int *timeInd, int *time0Ind, int *n, double *ret) {

  /* Use time0Ind[0] < -9999 for NULL */
  getRiskMatCol(*col, timeInd, time0Ind, *n, ret);
  return;

}

void C_getRiskMatCol2(int *col, int *timeInd, int *time0Ind, int *n, double *ret) {

  getRiskMatCol2(*col, timeInd, time0Ind, *n, ret);
  return;

}

static void getS0t(int *timeInd, int *time0Ind, double *expXwgt, int nsub, int ntimes, double *tmpvec, double *ret)
{
  /* ret should be length ntimes
     tmpvec must have length nsub
 */
  int i;
  for (i=0; i<ntimes; i++) {
    getRiskMatCol(i, timeInd, time0Ind, nsub, tmpvec);
    ret[i] = dotprod(tmpvec, expXwgt, nsub);
  }

}

void C_getS0t(int *timeInd, int *time0Ind, double *expXwgt, int *nsub, int *ntimes, double *ret) {

  double *tmpvec;

  tmpvec = dVec_alloc(*nsub, 0, 0.0);
  getS0t(timeInd, time0Ind, expXwgt, *nsub, *ntimes, tmpvec, ret);
  R_Free(tmpvec);

  return;
}


static void getS1t(int *timeInd, int *time0Ind, double *expXwgt, double *Xvec, int nsub, int ntimes, int p, double *tmpvec, double *tmpvec2, double **ret)
{
  /* X is nsub x p, stored as a vector by columns
     tmpvec, tmpvec2 length nsub
     ret has dim ntimes x p
  */
  int i, j;
  double *pX;

  for (i=0; i<ntimes; i++) {
    getRiskMatCol(i, timeInd, time0Ind, nsub, tmpvec);
    pX = Xvec;
    for (j=0; j<p; j++) {
      vecprod(pX, expXwgt, nsub, tmpvec2);
      ret[i][j] = dotprod(tmpvec, tmpvec2, nsub);
      pX += nsub;
    }
  }

}

void C_getS1t(int *timeInd, int *time0Ind, double *expXwgt, double *Xvec, int *nsub, int *ntimes, 
              int *p, double *ret)
{
  /* Xvec is nsub x p, stored as a vector by columns
     tmpvec, tmpvec2 length nsub
     ret has dim ntimes x p
  */
  double *tmpvec, *tmpvec2, **tmpmat;

  tmpvec  = dVec_alloc(*nsub, 0, 0.0);
  tmpvec2 = dVec_alloc(*nsub, 0, 0.0);
  tmpmat  = dMat_alloc(*ntimes, *p, 0, 0.0);

  getS1t(timeInd, time0Ind, expXwgt, Xvec, *nsub, *ntimes, *p, tmpvec, tmpvec2, tmpmat);
  putMatInVecByRow(tmpmat, *ntimes, *p, ret); /* By row */

  R_Free(tmpvec);
  R_Free(tmpvec2);
  dMat_free(tmpmat, *ntimes);

  return;
}

static void addInCrossProd(double *vec1, double *vec2, int n1, int n2, double **ret) 
{
  /* ret must have dimension n1 x n2 */
  int i, j;
  double vec1i;
  for (i=0; i<n1; i++) {
    vec1i = vec1[i];
    for (j=0; j<n2; j++) ret[i][j] += vec1i*vec2[j];
  }

}

static void getSumAAwgt(double *Aphase2ByRow, double *AwgtByRow, int nsub, int q, double **ret)
{
  /* ret mat dim q x q and already initialized to 0 */
  int i;
  double *p1, *p2;

  p1 = Aphase2ByRow;
  p2 = AwgtByRow;
  for (i=0; i<nsub; i++) {
    addInCrossProd(p1, p2, q, q, ret);
    p1 += q;
    p2 += q;
  }

}

void C_getSumAAwgt(double *Aphase2ByRow, double *AwgtByRow, int *nsub, int *q, double *ret) {

  double **tmpmat;

  tmpmat  = dMat_alloc(*q, *q, 1, 0.0); /* Initialize to 0 */
  getSumAAwgt(Aphase2ByRow, AwgtByRow, *nsub, *q, tmpmat);
  putMatInVecByRow(tmpmat, *q, *q, ret); /* By row */
  dMat_free(tmpmat, *q);

  return;
}

static void getdNtCol(int col, int *timeInd, int nsub, double *y, double *ret)
{
  int i;
  for (i=0; i<nsub; i++) {
    if (col == timeInd[i]) {
      ret[i] = y[i];
    } else {
      ret[i] = 0.0;
    }
  }
}

static void getdNtWgtCol(int col, int *timeInd, int nsub, double *y, double *wgt, double *ret)
{
  int i;
  for (i=0; i<nsub; i++) {
    if (col == timeInd[i]) {
      ret[i] = y[i]*wgt[i];
    } else {
      ret[i] = 0.0;
    }
  }
}

static void getYexpXWgtCol(int col, int *timeInd, int *time0Ind, int nsub, double *expXwgt, double *tmpvec, double *ret)
{

  getRiskMatCol(col, timeInd, time0Ind, nsub, tmpvec);
  vecprod(tmpvec, expXwgt, nsub, ret);

}

void C_getYexpXWgtCol(int *col, int *timeInd, int *time0Ind, int *nsub, double *expXwgt, double *ret) {

  double *tmpvec;
  tmpvec  = dVec_alloc(*nsub, 0, 0.0);

  getYexpXWgtCol(*col, timeInd, time0Ind, *nsub, expXwgt, tmpvec, ret);
  R_Free(tmpvec);

  return;

}

static void getdNtColSums(int *timeInd, int nsub, int ncol, double *y, double *tmpvec, double *ret)
{
  int i;
  for (i=0; i<ncol; i++) {
    getdNtCol(i, timeInd, nsub, y, tmpvec);
    ret[i] = vecsum(tmpvec, nsub);
  }
}

void C_getdNtColSums(int *timeInd, int *nsub, int *ncol, 
                       double *y, double *ret) {

  double *tmpvec;
  tmpvec  = dVec_alloc(*nsub, 0, 0.0);

  getdNtColSums(timeInd, *nsub, *ncol, y, tmpvec, ret);

  R_Free(tmpvec);
  return;
}

static void getdNtWgtColSums(int *timeInd, int nsub, int ncol, double *y, double *wgt, double *tmpvec, double *ret)
{
  int i;
  for (i=0; i<ncol; i++) {
    getdNtWgtCol(i, timeInd, nsub, y, wgt, tmpvec);
    ret[i] = vecsum(tmpvec, nsub);
  }
}

void C_getdNtWgtColSums(int *timeInd, int *nsub, int *ncol, 
                       double *y, double *wgt, double *ret) {

  double *tmpvec;
  tmpvec  = dVec_alloc(*nsub, 0, 0.0);

  getdNtWgtColSums(timeInd, *nsub, *ncol, y, wgt, tmpvec, ret);

  R_Free(tmpvec);
  return;
}

static void getXA(double *XbyRow, double *AbyRow, int p, int q, int n, double ***ret)
{
  /* ret must have dimension p x q x n 
     XA <- array(NA, dim = c(p, q, n))
     for(i in 1:n) {
        XA[,, i] <- tcrossprod(X[i,], A.phase2[i,])
     }
  */
  int i, j, k;
  double *pX, *pA, Xi;

  pX = XbyRow;
  pA = AbyRow;
  for (k=0; k<n; k++) {
    for (i=0; i<p; i++) {
      Xi = pX[i];
      for (j=0; j<q; j++) ret[i][j][k] = Xi*pA[j];
    }
    pX += p;
    pA += q; 
  }

}

static void getXexpectRow(int row, double *S1tByRow, int p, double *S0t, double *ret)
{
  /* S1tByRow ntimes x p, ntimes is not used
     ret has length p
  */
  int i;
  double *pS1, val0;
 
  pS1  = &S1tByRow[row*p];
  val0 = S0t[row];
  for (i=0; i<p; i++) ret[i] = pS1[i]/val0;

}

static void get_drond_U_eta(int DEBUG, double *XbyRow, double *AbyRow, int p, int q, int n, int ntimes, int *timeInd, int *time0Ind, double *y, 
                            double *S1tByRow, double *S0t, double *expXwgt, double *weights, double *dS0tEtaByCol, double *dUtEtaByRow)
{
  double ***XA, *tmpvecn, *tmp2vecn, **drondG1t, *drondG0t, **drondS1t, *drondS0t;
  double *csumsDntWgt, **tmpmatpq, **tmp2matpq, tmp1, tmp2, S0t_t, csum_t, *S1t_t;
  int t, i, j, ind;

  /* Allocate */
  if (DEBUG) Rprintf("Allocate memory\n");
  XA          = dArr_alloc(p, q, n);
  tmpvecn     = dVec_alloc(n, 0, 0.0);
  tmp2vecn    = dVec_alloc(n, 0, 0.0);
  drondG1t    = dMat_alloc(p, q, 0, 0.0);
  drondG0t    = dVec_alloc(q, 0, 0.0);
  drondS1t    = dMat_alloc(p, q, 0, 0.0);
  csumsDntWgt = dVec_alloc(ntimes, 0, 0.0);
  tmpmatpq    = dMat_alloc(p, q, 0, 0.0);
  tmp2matpq   = dMat_alloc(p, q, 0, 0.0);

  /* Initialize dUEtaByRow to 0 */
  if (DEBUG) Rprintf("init\n");
  initVec(dUtEtaByRow, p*q, 0.0);

  /* Compute colSums(dNt.weighted),  dNt.weighted is n x ntimes */
  if (DEBUG) Rprintf("colsums\n");
  getdNtWgtColSums(timeInd, n, ntimes, y, weights, tmpvecn, csumsDntWgt);

  /* Get XA array array(NA, dim = c(p, q, n)) */
  if (DEBUG) Rprintf("getXA\n");
  getXA(XbyRow, AbyRow, p, q, n, XA); 

  for (t=0; t<ntimes; t++) {
    if (DEBUG) Rprintf("t = %d, ntimes=%d\n", t, ntimes);

    /* Compute drond.G1t.eta[,, t] */
    if (DEBUG) Rprintf("getdNtWgtCol\n");
    getdNtWgtCol(t, timeInd, n, y, weights, tmpvecn);
    if (DEBUG) Rprintf("applySweep3sum\n");
    applySweep3sum(XA, p, q, n, tmpvecn, drondG1t);

    /* Compute drond.G0t.eta[, t] */
    if (DEBUG) Rprintf("colSumsRowMatTimesVec\n");
    colSumsRowMatTimesVec(AbyRow, n, q, tmpvecn, tmp2vecn, drondG0t);

    /* Compute drond.S1t.eta[,, t] */
    if (DEBUG) Rprintf("getYexpWgtCol\n");
    getYexpXWgtCol(t, timeInd, time0Ind, n, expXwgt, tmpvecn, tmp2vecn);
    if (DEBUG) Rprintf("applySweep3sum\n");
    applySweep3sum(XA, p, q, n, tmp2vecn, drondS1t);

    /* Compute drond.S0t.eta[, t] dim q x ntimes */
    if (DEBUG) Rprintf("drondS0t\n");
    drondS0t = &dS0tEtaByCol[t*q]; 

    if (DEBUG) Rprintf("colSumsRowMatTimesVec\n");
    colSumsRowMatTimesVec(AbyRow, n, q, tmp2vecn, tmpvecn, drondS0t);

    if (DEBUG) Rprintf("getXepectRow\n");
    getXexpectRow(t, S1tByRow, p, S0t, tmpvecn); /* X.expect[t, ] */

    if (DEBUG) Rprintf("tcrossprod\n");
    tcrossprod(tmpvecn, drondG0t, p, q, tmpmatpq); /* tcrossprod(X.expect[t, ], drond.G0t.eta[, t]) */

    S0t_t  = S0t[t];
    csum_t = csumsDntWgt[t];
    tmp1   = csum_t/S0t_t;
    tmp2   = csum_t/(S0t_t*S0t_t);
    S1t_t  = &S1tByRow[t*p];
    if (DEBUG) Rprintf("tcrossprod\n");
    tcrossprod(S1t_t, drondS0t, p, q, tmp2matpq); /* tcrossprod(S1t[t, ], drond.S0t.eta[, t]) */
    ind = 0;
    for (i=0; i<p; i++) {
      for (j=0; j<q; j++) {
        dUtEtaByRow[ind] += drondG1t[i][j] - tmpmatpq[i][j] - tmp1*drondS1t[i][j] +
                           tmp2*tmp2matpq[i][j];
        ind++;
      }
    }
  }

  if (DEBUG) Rprintf("FREE XA\n");
  dArr_free(XA, p, q);
  if (DEBUG) Rprintf("FREE tempvecn\n");
  R_Free(tmpvecn);
  if (DEBUG) Rprintf("FREE tempvec2n\n");
  R_Free(tmp2vecn);
  if (DEBUG) Rprintf("FREE drondG1t\n");
  dMat_free(drondG1t, p);
  if (DEBUG) Rprintf("FREE drondG0t\n");
  R_Free(drondG0t);
  if (DEBUG) Rprintf("FREE drondS1t\n");
  dMat_free(drondS1t, p);
  if (DEBUG) Rprintf("FREE csumsDntWgt\n");
  R_Free(csumsDntWgt);
  if (DEBUG) Rprintf("FREE tmpmatpq\n");
  dMat_free(tmpmatpq, p);
  if (DEBUG) Rprintf("FREE tmp2matpq\n");
  dMat_free(tmp2matpq, p);

}

void C_get_drond_U_eta(int *DEBUG, double *XbyRow, double *AbyRow, int *p, int *q, int *n, int *ntimes, 
                       int *timeInd, int *time0Ind, double *y, double *S1tByRow, double *S0t, double *expXwgt, 
                       double *weights, double *ret_dS0tEtaByCol, double *ret_dUtEtaByRow) {

  get_drond_U_eta(*DEBUG, XbyRow, AbyRow, *p, *q, *n, *ntimes, timeInd, time0Ind, y, 
                  S1tByRow, S0t, expXwgt, weights, ret_dS0tEtaByCol, ret_dUtEtaByRow);
  return;
}

static void initScoreForDim3(int dim3, int nsub, int p, double *XeByRow, double *XbyRow, double **ret)
{
  /* ret must be dim nsub x p 
     0 <= dim3 <= ntimes-1  
     score <- array(NA, dim = c(n, p, number.times))
     for(i in 1:n) score[i,,] <- - sweep(t(X.expect), 1, X[i,], "-")

     Return the nxp matrix for the third dimension at dim3
  */
  int i, j;
  double *pXe, *pX;

  pXe = XeByRow + dim3*p;
  pX  = XbyRow;
  for (i=0; i<nsub; i++) {
    for (j=0; j<p; j++) ret[i][j] = -(pXe[j] - pX[j]);
    pX += p; 
  }
}

static void getBetaScore(int nsub, int p, double *XeByRow, double *XbyRow, int *timeInd, int *time0Ind, int ntimes,
                            double *y, double *weights, double *expXwgt, double *csumS0, double *retByRow)
{
  /* retByRow must be length nsub x p
  */

  double *csumsDntWgt, **score, *tmpvecn, *tmp2vecn, dMt, csumS0t;
  int t, i, j, ind;

  initVec(retByRow, nsub*p, 0.0);

  score       = dMat_alloc(nsub, p, 0, 0.0);
  csumsDntWgt = dVec_alloc(ntimes, 0, 0.0);
  tmpvecn     = dVec_alloc(nsub, 0, 0.0);
  tmp2vecn    = dVec_alloc(nsub, 0, 0.0);

  /* Compute colSums(dNt.weighted) */
  getdNtWgtColSums(timeInd, nsub, ntimes, y, weights, tmpvecn, csumsDntWgt);

  for (t=0; t<ntimes; t++) { 
    /* Initialize score to - sweep(t(X.expect), 1, X[i,], "-") */
    initScoreForDim3(t, nsub, p, XeByRow, XbyRow, score);

    getdNtWgtCol(t, timeInd, nsub, y, weights, tmpvecn);
    getRiskMatCol(t, timeInd, time0Ind, nsub, tmp2vecn);
  
    csumS0t = csumS0[t];
    for (i=0; i<nsub; i++) {
      dMt = tmpvecn[i] - tmp2vecn[i]*expXwgt[i]*csumS0t;
      for (j=0; j<p; j++) score[i][j] = score[i][j]*dMt;
    }
    /* Sum up over time points */
    ind = 0;
    for (i=0; i<nsub; i++) {
      for (j=0; j<p; j++) retByRow[ind++] += score[i][j];
    }
  }
 
  dMat_free(score, nsub);
  R_Free(csumsDntWgt);
  R_Free(tmpvecn);
  R_Free(tmp2vecn);

}

void C_getBetaScore(int *nsub, int *p, double *XeByRow, double *XbyRow, 
                           int *timeInd, int *time0Ind, int *ntimes, double *y, double *weights, 
                           double *expXwgt, double *csumS0, double *retByRow) {

  getBetaScore(*nsub, *p, XeByRow, XbyRow, timeInd, time0Ind, *ntimes,
                            y, weights, expXwgt, csumS0, retByRow);
  return;
}

static void infl_lambda0_tau12(double *infl1BetaByRow, double *S1tByRow, double *infl1EtaByRow, double *dS0tByCol,
                          double *infl2BetaByRow, double *infl2EtaByRow,
                          double *lambda0, double *S0, int *timeInd, int *time0Ind, double *y, double *expXwgt, int nsub, int nsubA, int ntimes, int p, int q,
                          int *Tau12Times, int nTau12Times, double *ret1, double *ret2) 
{
  /* ret1 has length nsubA, ret2 has length nsub 
     Loop over Tau12 columns to compute the rowsums for infl1 and infl2
     Tau12Times must have 1 subtracted already (C vs R indices)
  */
  int t, col, i;  
  double *dNt, *YexpXwgt, *tmpvecn, *pBeta, *pS1, *pEta, *pdS0;
  double lambda0t, S0t, val1, val2, infl1, infl2, tmp;

  initVec(ret1, nsubA, 0.0);
  initVec(ret2, nsub, 0.0);

  dNt      =  dVec_alloc(nsub, 0, 0.0);
  YexpXwgt =  dVec_alloc(nsub, 0, 0.0);
  tmpvecn  =  dVec_alloc(nsub, 0, 0.0);

  for (t=0; t<nTau12Times; t++) {
    col       = Tau12Times[t];
    lambda0t  = lambda0[col];
    S0t       = S0[col];
    tmp       = lambda0t/S0t;
    pBeta     = infl1BetaByRow;
    pEta      = infl1EtaByRow;
    pS1       = S1tByRow + col*p;
    pdS0      = dS0tByCol + col*q;

    for (i=0; i<nsubA; i++) {
      val1    = dotprod(pBeta, pS1, p);
      val2    = dotprod(pEta, pdS0, q);
      infl1   = -(val1 + val2)*tmp;
      ret1[i] = ret1[i] + infl1;
      pBeta  += p;
      pEta   += q;
    }

    getdNtCol(col, timeInd, nsub, y, dNt);
    getYexpXWgtCol(col, timeInd, time0Ind, nsub, expXwgt, tmpvecn, YexpXwgt);
    pBeta     = infl2BetaByRow;
    pEta      = infl2EtaByRow;
    pS1       = S1tByRow + col*p;
    pdS0      = dS0tByCol + col*q;
    for (i=0; i<nsub; i++) {
      val1    = dotprod(pBeta, pS1, p);
      val2    = dotprod(pEta, pdS0, q);
      infl2   = (dNt[i] - (YexpXwgt[i] + val1 + val2)*lambda0t)/S0t;
      ret2[i] = ret2[i] + infl2;
      pBeta  += p;
      pEta   += q;
    }
  }

  R_Free(dNt);
  R_Free(YexpXwgt);
  R_Free(tmpvecn);

}

void C_infl_lambda0_tau12(double *infl1BetaByRow, double *S1tByRow, double *infl1EtaByRow, 
                          double *dS0tByCol, double *infl2BetaByRow, double *infl2EtaByRow,
                          double *lambda0, double *S0, int *timeInd, int *time0Ind, double *y, double *expXwgt, 
                          int *nsub, int *nsubA, int *ntimes, int *p, int *q,
                          int *Tau12Times, int *nTau12Times, 
                          double *ret1, double *ret2) {

  infl_lambda0_tau12(infl1BetaByRow, S1tByRow, infl1EtaByRow, dS0tByCol,
                          infl2BetaByRow, infl2EtaByRow,
                          lambda0, S0, timeInd, time0Ind, y, expXwgt, *nsub, *nsubA, *ntimes, *p, *q,
                          Tau12Times, *nTau12Times, ret1, ret2);
  return;
} 

static void infl_lambda0_tau12_noCalib(double *inflBetaByRow, double *S1tByRow, 
                          double *lambda0, double *S0, int *timeInd, int *time0Ind, double *y, double *expXwgt, int nsub, int ntimes, int p,
                          int *Tau12Times, int nTau12Times, double *ret)
{
  /* ret has length nsub 
     Loop over Tau12 columns to compute the rowsums for infl
     Tau12Times must have 1 subtracted already (C vs R indices)
  */
  int t, col, i;  
  double *dNt, *YexpXwgt, *tmpvecn, *pBeta, *pS1;
  double lambda0t, S0t, val, infl;

  initVec(ret, nsub, 0.0);

  dNt      =  dVec_alloc(nsub, 0, 0.0);
  YexpXwgt =  dVec_alloc(nsub, 0, 0.0);
  tmpvecn  =  dVec_alloc(nsub, 0, 0.0);

  for (t=0; t<nTau12Times; t++) {
    col       = Tau12Times[t];
    lambda0t  = lambda0[col];
    S0t       = S0[col];
    
    getdNtCol(col, timeInd, nsub, y, dNt);
    getYexpXWgtCol(col, timeInd, time0Ind, nsub, expXwgt, tmpvecn, YexpXwgt);
    pBeta     = inflBetaByRow;
    pS1       = S1tByRow + col*p;
    for (i=0; i<nsub; i++) {
      val    = dotprod(pBeta, pS1, p);
      infl   = (dNt[i] - (YexpXwgt[i] + val)*lambda0t)/S0t;
      ret[i] = ret[i] + infl;
      pBeta  += p;
    }
  }

  R_Free(dNt);
  R_Free(YexpXwgt);
  R_Free(tmpvecn);

}

void C_infl_lambda0_tau12_noCalib(double *inflBetaByRow, double *S1tByRow, 
                          double *lambda0, double *S0, int *timeInd, int *time0Ind, double *y, 
                          double *expXwgt, int *nsub, int *ntimes, int *p,
                          int *Tau12Times, int *nTau12Times, double *ret) {

  infl_lambda0_tau12_noCalib(inflBetaByRow, S1tByRow, 
                          lambda0, S0, timeInd, time0Ind, y, expXwgt, *nsub, *ntimes, *p,
                          Tau12Times, *nTau12Times, ret);
  return;

} 

static void getProdCovarWgtRow_strat(int *W, int row, double *omegaConst, double *covarConst, double *omegaDiag, double *covarDiag,
                               int n, double *ret)
{
  /* ret must be length n 
     W must be coded 0, 1, 2, ..., m
     omegaConst, etc coded for each stratum (in order)
  */

  int i, Wrow;
  double omega, covar;

  Wrow = W[row];
  for (i=0; i<n; i++) {
    if (i == row) { 
      /* Diag elements */
      omega = omegaDiag[Wrow];
      covar = covarDiag[Wrow];
    } else if (W[i] == Wrow) {
      omega = omegaConst[Wrow];
      covar = covarConst[Wrow];
    } else {
      omega = 0.0;
      covar = 0.0;
    }
    if (omega < -0.5) omega = 0.0;
    if (covar < -0.5) covar = 0.0;
    ret[i] = omega*covar; 
  }

}

static void getProdCovarWgtRow_unstrat(int row, double omegaConst, double covarConst, double omegaDiag, double covarDiag,
                               int n, double *ret)
{
  /* ret must be length n 
  */

  int i;
  double omega, covar;

  for (i=0; i<n; i++) {
    if (i != row) { 
      omega = omegaConst;
      covar = covarConst; 
    } else {
      omega = omegaDiag;
      covar = covarDiag;
    }
    ret[i] = omega*covar; 
  }

}

static void getPhase2Var(int stratFlag, int *W, double *omegaConst, double *covarConst, double *omegaDiag, double *covarDiag,
                         int n, double *infl2ByCol, int p, double *ret)
{
  /* ret has length p
     If stratFlag = 0, then omegaConst, etc will be scalars
     Loop over rows of prod.covar.weight to compute 
                                prod.covar.weight %*% 
                                infl2[(status == 0) & (subcohort == TRUE),])
     prod.covar.weight is an nxn matrix
  */

  double **mat, *vecn, *pinfl2, tmp;
  int i, j;

  mat  =  dMat_alloc(n, p, 0, 0.0);
  vecn =  dVec_alloc(n, 0, 0.0);

  for (i=0; i<n; i++) {
    /* Get a row of prod.covar.weight */
    if (stratFlag) {
      getProdCovarWgtRow_strat(W, i, omegaConst, covarConst, omegaDiag, covarDiag, n, vecn);
    } else {
      getProdCovarWgtRow_unstrat(i, *omegaConst, *covarConst, *omegaDiag, *covarDiag, n, vecn);
    }

    pinfl2 = infl2ByCol;
    for (j=0; j<p; j++) {
      mat[i][j] = dotprod(vecn, pinfl2, n);
      pinfl2 += n;
    }
  }
  R_Free(vecn);
  
  /* Now compute the quadratic form to get the diagonal elements */
  pinfl2 = infl2ByCol;
  for (j=0; j<p; j++) {
    tmp = 0.0;
    for (i=0; i<n; i++) tmp += pinfl2[i]*mat[i][j];
    ret[j] = tmp;
    pinfl2 += n;
  }

  dMat_free(mat, n);
  
}

void C_getPhase2Var(int *stratFlag, int *W, double *omegaConst, double *covarConst, 
                    double *omegaDiag, double *covarDiag, int *n, double *infl2ByCol, 
                    int *p, double *ret) {

  getPhase2Var(*stratFlag, W, omegaConst, covarConst, omegaDiag, covarDiag,
               *n, infl2ByCol, *p, ret);
  return;

}

/******************************************************************
  Missing data
  phase3 is a subset of phase2
  phase2 risk matrix is from the whole cohort fit
*******************************************************************/

static void getS0GammaCasetimes(double *B3ByCol, int *obstimes3, int *obstimes3_0, int neventTimes3, double *expXwgt,
                                int n3, int J3, double *retByCol)
{
  /*
  ret must have length J3 x neventTimes3
  drond.S0t.gamma.casestimes[, t] <- colSums(B.phase3 * matrix(Y.exp.X.weighted.casestimes[, t], 
                                                                   nrow = n.phase3, 
                                                                   ncol = J3, 
                                                                   byrow = FALSE))
  */

  int t;
  double *tmpvec, *yexpXwgt, *pret;

  tmpvec   =  dVec_alloc(n3, 0, 0.0);
  yexpXwgt =  dVec_alloc(n3, 0, 0.0);

  /* ret must be initialized to 0 */
  initVec(retByCol, neventTimes3*J3, 0.0);

  /* Loop over each event time */
  pret = retByCol;
  for (t=0; t<neventTimes3; t++) {
    /* Get Y.exp.X.weighted.casestimes[, t] = sweep(riskmat.phase3, 1, exp.X.weighted, "*")[, t]
    */
    getYexpXWgtCol(t, obstimes3, obstimes3_0, n3, expXwgt, tmpvec, yexpXwgt);
    colSumsColMatTimesVec(B3ByCol, n3, J3, yexpXwgt, pret);
    pret += J3;
  }

  R_Free(yexpXwgt);
  R_Free(tmpvec);
}

void C_getS0GammaCasetimes(double *B3ByCol, int *obstimes3, int *obstimes3_0, int *neventTimes3, double *expXwgt,
                           int *n3, int *J3, double *retByCol) {

  getS0GammaCasetimes(B3ByCol, obstimes3, obstimes3_0, *neventTimes3, expXwgt,
                      *n3, *J3, retByCol);
  return;

}

static void infl2_lambda0t(int *obstimes2, double *y2, double *infl2BetaByRow, double *S1tCaseByRow,
                           double *infl2GammaByRow, double *dS0tGammaCaseByCol, double *S0tCase,
                           double *lambda0, int *Tau12Times, int nTau12Times, int n2, int p, int J3, double *ret)
{
  /* ret has length n2
     Loop over Tau12 columns to compute the rowsums for infl2
     Tau12Times must have 1 subtracted already (C vs R indices)
  */
  int t, col, i;  
  double *dNt2, *pBeta, *pS1, *pGamma, *pdS0;
  double lambda0t, S0t, val1, val2, infl2;

  initVec(ret, n2, 0.0);

  dNt2 =  dVec_alloc(n2, 0, 0.0);

  for (t=0; t<nTau12Times; t++) {
    col       = Tau12Times[t];
    lambda0t  = lambda0[col];
    S0t       = S0tCase[col];

    getdNtCol(col, obstimes2, n2, y2, dNt2);
    pBeta     = infl2BetaByRow;
    pGamma    = infl2GammaByRow;
    pS1       = S1tCaseByRow + col*p;
    pdS0      = dS0tGammaCaseByCol + col*J3;
    for (i=0; i<n2; i++) {
      val1    = dotprod(pBeta, pS1, p);
      val2    = dotprod(pGamma, pdS0, J3);
      infl2   = (dNt2[i] - (val1 + val2)*lambda0t)/S0t;
      ret[i]  = ret[i] + infl2;
      pBeta  += p;
      pGamma += J3;
    }    
  }

  R_Free(dNt2);

}

void C_infl2_lambda0t(int *obstimes2, double *y2, double *infl2BetaByRow, double *S1tCaseByRow,
                      double *infl2GammaByRow, double *dS0tGammaCaseByCol, double *S0tCase,
                      double *lambda0, int *Tau12Times, int *nTau12Times, int *n2, int *p, int *J3, 
                      double *ret) {

  infl2_lambda0t(obstimes2, y2, infl2BetaByRow, S1tCaseByRow,
                 infl2GammaByRow, dS0tGammaCaseByCol, S0tCase,
                 lambda0, Tau12Times, *nTau12Times, *n2, *p, *J3, 
                 ret);

  return;
}

static void infl3_lambda0t(int *obstimes3, int *obstimes3_0, double *expXwgt, double *infl3BetaByRow, double *S1tCaseByRow,
                           double *infl3GammaByRow, double *dS0tGammaCaseByCol, double *S0tCase,
                           double *lambda0, int *Tau12Times, int nTau12Times, int n3, int p, int J3, double *ret)
{
  /* ret has length n3
     Loop over Tau12 columns to compute the rowsums for infl3
     Tau12Times must have 1 subtracted already (C vs R indices)
  */
  int t, col, i;  
  double *yexpwgt, *tmpvec, *pBeta, *pS1, *pGamma, *pdS0;
  double lambda0t, S0t, val1, val2, infl3, tmp;

  initVec(ret, n3, 0.0);

  yexpwgt =  dVec_alloc(n3, 0, 0.0);
  tmpvec  =  dVec_alloc(n3, 0, 0.0);

  for (t=0; t<nTau12Times; t++) {
    col       = Tau12Times[t];
    lambda0t  = lambda0[col];
    S0t       = S0tCase[col];
    tmp       = lambda0t/S0t;

    /* Get column col of matrix Y.exp.X.weighted.casestimes */
    getYexpXWgtCol(col, obstimes3, obstimes3_0, n3, expXwgt, tmpvec, yexpwgt);

    pBeta     = infl3BetaByRow;
    pGamma    = infl3GammaByRow;
    pS1       = S1tCaseByRow + col*p;
    pdS0      = dS0tGammaCaseByCol + col*J3;
    for (i=0; i<n3; i++) {
      val1    = dotprod(pBeta, pS1, p);
      val2    = dotprod(pGamma, pdS0, J3);
      infl3   = -(yexpwgt[i] + val1 + val2)*tmp;
      ret[i]  = ret[i] + infl3;
      pBeta  += p;
      pGamma += J3;
    }    
  }

  R_Free(yexpwgt);
  R_Free(tmpvec);

}

void C_infl3_lambda0t(int *obstimes3, int *obstimes3_0, double *expXwgt, double *infl3BetaByRow, double *S1tCaseByRow,
                    double *infl3GammaByRow, double *dS0tGammaCaseByCol, double *S0tCase,
                    double *lambda0, int *Tau12Times, int *nTau12Times, int *n3, 
                    int *p, int *J3, double *ret) {

  infl3_lambda0t(obstimes3, obstimes3_0, expXwgt, infl3BetaByRow, S1tCaseByRow,
                           infl3GammaByRow, dS0tGammaCaseByCol, S0tCase,
                           lambda0, Tau12Times, *nTau12Times, *n3, *p, *J3, ret);
  return;

}

static void infl2_lambda0t_noEst(int *obstimes2, double *y2, double *S0tCase, int *Tau12Times, int nTau12Times, int n2, double *ret)
{
  /* ret has length n2
     Loop over Tau12 columns to compute the rowsums for infl2
     Tau12Times must have 1 subtracted already (C vs R indices)
  */
  int t, col, i;  
  double *dNt2, S0t;

  initVec(ret, n2, 0.0);

  dNt2 =  dVec_alloc(n2, 0, 0.0);

  for (t=0; t<nTau12Times; t++) {
    col = Tau12Times[t];
    S0t = S0tCase[col];

    getdNtCol(col, obstimes2, n2, y2, dNt2);

    for (i=0; i<n2; i++) ret[i] = ret[i] + dNt2[i]/S0t;
  }

  R_Free(dNt2);

}
void C_infl2_lambda0t_noEst(int *obstimes2, double *y2, double *S0tCase, 
                     int *Tau12Times, int *nTau12Times, int *n2, double *ret) {

  infl2_lambda0t_noEst(obstimes2, y2, S0tCase, Tau12Times, *nTau12Times, *n2, ret);

  return;
}

static void infl3_lambda0t_noEst(int *obstimes3, int *obstimes3_0, double *expXwgt, double *infl3BetaByRow, double *S1tCaseByRow,
                           double *S0tCase, double *lambda0, int *Tau12Times, int nTau12Times, int n3, int p, double *ret)
{
  /* ret has length n3
     Loop over Tau12 columns to compute the rowsums for infl3
     Tau12Times must have 1 subtracted already (C vs R indices)
  */
  int t, col, i;  
  double *yexpwgt, *tmpvec, *pBeta, *pS1;
  double lambda0t, S0t, val1, infl3, tmp;

  initVec(ret, n3, 0.0);

  yexpwgt =  dVec_alloc(n3, 0, 0.0);
  tmpvec  =  dVec_alloc(n3, 0, 0.0);

  for (t=0; t<nTau12Times; t++) {
    col       = Tau12Times[t];
    lambda0t  = lambda0[col];
    S0t       = S0tCase[col];
    tmp       = lambda0t/S0t;

    /* Get column col of matrix Y.exp.X.weighted.casestimes */
    getYexpXWgtCol(col, obstimes3, obstimes3_0, n3, expXwgt, tmpvec, yexpwgt);

    pBeta     = infl3BetaByRow;
    pS1       = S1tCaseByRow + col*p;
    for (i=0; i<n3; i++) {
      val1    = dotprod(pBeta, pS1, p);
      infl3   = -(yexpwgt[i] + val1)*tmp;
      ret[i]  = ret[i] + infl3;
      pBeta  += p;
    }    
  }

  R_Free(yexpwgt);
  R_Free(tmpvec);

}

void C_infl3_lambda0t_noEst(int *obstimes3, int *obstimes3_0, double *expXwgt, double *infl3BetaByRow, 
                            double *S1tCaseByRow, double *S0tCase, double *lambda0, 
                            int *Tau12Times, int *nTau12Times, int *n3, int *p, double *ret) {

  infl3_lambda0t_noEst(obstimes3, obstimes3_0, expXwgt, infl3BetaByRow, S1tCaseByRow,
                           S0tCase, lambda0, Tau12Times, *nTau12Times, *n3, *p, ret);
  return;
}

static void prodCovWgtStratT_row(int *W, int *W3, int *y, int row, int n, double *valVec, double *transWgt, double *ret)
{
  /* ret must have length n (row row of nxn matrix)
  */
  int i, Wrow, W3row;

  Wrow  = W[row];
  W3row = W3[row];
  for (i=0; i<n; i++) {
    ret[i] = 0.0;
    if (i != row) {
      if (!y[i] && (W[i] == Wrow) && (W3[i] == W3row)) ret[i] = valVec[Wrow]; 
    } else {
      ret[i] = transWgt[i];
    }
  }

} 

void C_prodCovWgtStratT_row(int *W, int *W3, int *y, int *row, int *n, 
          double *valVec, double *transWgt, double *ret) {

  prodCovWgtStratT_row(W, W3, y, *row, *n, valVec, transWgt, ret);
  return;

}

static void prodCovWgtStratF_row(int *W3, int *y, int row, int n, double val, double *transWgt, double *ret)
{
  /* ret must have length n (row row of nxn matrix)
     unstratified, 1 value
  */
  int i, W3row;

  W3row = W3[row];
  for (i=0; i<n; i++) {
    ret[i] = 0.0;
    if (i != row) {
      if (!y[i] && (W3[i] == W3row)) ret[i] = val; 
    } else {
      ret[i] = transWgt[i];
    }
  }

} 

void C_prodCovWgtStratF_row(int *W3, int *y, int *row, int *n, 
                          double *val, double *transWgt, double *ret) {

  prodCovWgtStratF_row(W3, y, *row, *n, *val, transWgt, ret);
  return;

}

static void phase23VarEstF(int DEBUG, int *W, int *W3, int *y, int n, double *transWgt, double *valVec, double *infl3ByCol, int p, int stratFlag, double *ret)
{
  /* Note that y is integer in this function !!!

  */
  double **mat, *vecn, *pinfl3, tmp, val0;
  int i, j;

  mat  =  dMat_alloc(n, p, 0, 0.0);
  vecn =  dVec_alloc(n, 0, 0.0);
  val0 = valVec[0];

  for (i=0; i<n; i++) {
    /* Get a row of prod.covar.weight */
    if (stratFlag) {
      if (DEBUG) Rprintf("prodCovWgtStratT_row, i=%d\n", i); 
      prodCovWgtStratT_row(W, W3, y, i, n, valVec, transWgt, vecn);
    } else {
      if (DEBUG) Rprintf("prodCovWgtStratF_row, i=%d\n", i);
      prodCovWgtStratF_row(W3, y, i, n, val0, transWgt, vecn);
    }

    /* Multiply row by column of infl3 */
    pinfl3 = infl3ByCol;
    for (j=0; j<p; j++) {
      if (DEBUG) Rprintf("dotprod, j=%d\n", j);
      mat[i][j] = dotprod(vecn, pinfl3, n);
      pinfl3 += n;
    }
  }
  if (DEBUG) Rprintf("FREE vecn\n");
  R_Free(vecn);
  
  /* Now compute the quadratic form to get the diagonal elements */
  pinfl3 = infl3ByCol;
  for (j=0; j<p; j++) {
    if (DEBUG) Rprintf("quadform, j=%d\n", j);
    tmp = 0.0;
    for (i=0; i<n; i++) tmp += pinfl3[i]*mat[i][j];
    ret[j] = tmp;
    pinfl3 += n;
  }
  if (DEBUG) Rprintf("FREE mat\n");
  dMat_free(mat, n);

}

void C_phase23VarEstF(int *DEBUG, int *W, int *W3, int *y, int *n, double *transWgt, double *valVec, 
                     double *infl3ByCol, int *p, int *stratFlag, 
                     double *ret) {

  phase23VarEstF(*DEBUG, W, W3, y, *n, transWgt, valVec, infl3ByCol, *p, *stratFlag, ret);
  return;

}

static void phase23VarEstT(int DEBUG, int *W, int *W3, int *y, int n, double *transWgt, double *transWgt2, double *valVec, 
                           double *infl2ByCol, double *infl3ByCol, int p, int stratFlag, double *ret)
{
  /* Note that y is integer in this function !!!
     transWgt  is 1 - 1 / weights.phase2 
     transWgt2 is 1 / weights.p2.phase2 * (1 - 1 / weights.p3.phase2)
  */
  double **mat, *vecn, *pinfl2, *pinfl3, quad1, quad2, val0;
  int i, j;

  mat  =  dMat_alloc(n, p, 0, 0.0);
  vecn =  dVec_alloc(n, 0, 0.0);
  val0 = valVec[0];

  for (i=0; i<n; i++) {
    /* Get a row of prod.covar.weight */
    if (stratFlag) {
      prodCovWgtStratT_row(W, W3, y, i, n, valVec, transWgt, vecn);
    } else {
      prodCovWgtStratF_row(W3, y, i, n, val0, transWgt, vecn);
    }

    /* Multiply row by column of infl2 + infl3 */
    pinfl2 = infl2ByCol;
    pinfl3 = infl3ByCol;
    for (j=0; j<p; j++) {
      mat[i][j] = dotprodV1plusV2(vecn, pinfl2, pinfl3, n);
      pinfl2 += n;
      pinfl3 += n;
    }
  }
  R_Free(vecn);
  
  /* Now compute the quadratic form to get the diagonal elements */
  pinfl2 = infl2ByCol;
  pinfl3 = infl3ByCol;
  for (j=0; j<p; j++) {
    quad1 = 0.0;
    for (i=0; i<n; i++) quad1 += (pinfl2[i] + pinfl3[i])*mat[i][j];
    quad2 = 0.0;
    for (i=0; i<n; i++) quad2 += pinfl2[i]*transWgt2[i]*(pinfl2[i] + 2.0*pinfl3[i]);
    ret[j] = quad1 - quad2;
    pinfl2 += n;
    pinfl3 += n;
  }

  dMat_free(mat, n);

}

void C_phase23VarEstT(int *DEBUG, int *W, int *W3, int *y, int *n, double *transWgt, double *transWgt2,
                      double *valVec, double *infl2ByCol, double *infl3ByCol, int *p, 
                      int *stratFlag, double *ret) {

  phase23VarEstT(*DEBUG, W, W3, y, *n, transWgt, transWgt2, valVec, infl2ByCol, infl3ByCol, 
                 *p, *stratFlag, ret);
  return;

}


