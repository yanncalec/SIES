#include "mex.h"
#include <math.h>
#include <omp.h>

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
  int px, py, ax, ay, ordmin, ordmax;
  size_t Ntx, Nty, nrows, ncols;
  double x0, x1, y0, y1, dhx, dhy;
  double *Tphi, *Tpsi, *toto, *Z;

  /* check for proper number of arguments */
  if(nrhs!=6) {
    mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nrhs","Seven inputs required.");
  }
  if(nlhs!=1) {
    mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nlhs","One output required.");
  }
    
  Tphi = mxGetPr(prhs[0]); /* Table of phi */  
  toto = mxGetPr(prhs[1]); 
  x0 = toto[0]; x1 = toto[1]; /* support range of phi */
  Ntx = mxGetNumberOfElements(prhs[0]); /* Number of table values */
  dhx = (x1-x0)/(double)Ntx; /* Sampling step */

  Tpsi = mxGetPr(prhs[2]); /* Table of psi */  
  toto = mxGetPr(prhs[3]); 
  y0 = toto[0]; y1 = toto[1]; /* support range of psi */
  Nty = mxGetNumberOfElements(prhs[2]); /* Number of table values */
  dhy = (y1-y0)/(double)Nty; /* Sampling step */

  toto = mxGetPr(prhs[4]); ordmin = (int)(*toto);
  toto = mxGetPr(prhs[5]); ordmax = (int)(*toto);


  /* mexPrintf("x0,x1,y0,y1:%f,%f,%f,%f\n",x0,x1,y0,y1); */
  /* mexPrintf("Sampling step:%f\n",dh); */
  /* mexPrintf("ordmin, ordmax, Nt:%d,%d,%d\n",ordmin, ordmax, Nt); */

  nrows=ordmax-ordmin+1;
  ncols=ordmax-ordmin+1;

  plhs[0] = mxCreateDoubleMatrix(nrows, ncols, mxREAL);

  double *A = mxGetPr(plhs[0]);
  size_t idx;

  /* Loop for all position indexes in the range of rx X ry */
#pragma omp parallel for shared(A,Tphi,Tpsi,x0,y0) private(ax, ay, px,py,idx)
  for(ay=ordmin; ay<=ordmax; ay++){ /* column iteration */
    for(ax=ordmin; ax<=ordmax; ax++){ /* row iteration */

      idx = (ay-ordmin) * nrows + (ax-ordmin);
      A[idx] = 0;

      /* Loop inside the support of the wavelet for the inner product*/
      for(px=0; px<Ntx; px++){
	for(py=0; py<Nty; py++){
	  A[idx] += Tphi[px]*Tpsi[py] * pow(px*dhx+x0, ax)*pow(py*dhy+y0, ay);
	}
      }
      A[idx] *= (dhx*dhy);
    }
  }
}
