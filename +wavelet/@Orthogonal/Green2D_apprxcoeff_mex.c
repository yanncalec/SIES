#include "mex.h"
#include <math.h>
/* #include <omp.h> */

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
  int rx0, rx1, ry0, ry1, Jmax, row, col, px, py, nx, ny;
  size_t Nt, nrows, ncols;
  double t0, t1, dh, xx, yy;
  double *Tphi, *toto, *Z, *out;

  /* check for proper number of arguments */
  if(nrhs!=6) {
    mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nrhs","Seven inputs required.");
  }
  if(nlhs!=1) {
    mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nlhs","One output required.");
  }
    
  Tphi = mxGetPr(prhs[0]); /* Table of phi */  
  Nt = mxGetNumberOfElements(prhs[0]); /* Number of table values */

  toto = mxGetPr(prhs[1]); 
  t0 = toto[0]; t1 = toto[1]; /* support range of the scaling function phi */
  dh = (t1-t0)/(double)Nt; /* Sampling step */

  toto = mxGetPr(prhs[2]); Jmax = (int)(*toto); 

  toto = mxGetPr(prhs[3]); 
  rx0 = (int)(toto[0]); rx1 = (int)(toto[1]); /* range in x-axis of the position index */

  toto = mxGetPr(prhs[4]); 
  ry0 = (int)(toto[0]); ry1 = (int)(toto[1]); /* range in x-axis of the position index */

  Z = mxGetPr(prhs[5]); /* singularity of the Green function */

  nrows = ry1-ry0+1; ncols = rx1-rx0+1;

  /* mexPrintf("t0,t1:%f,%f\n",t0,t1); */
  /* mexPrintf("rx, ry, Jmax:%d,%d,%d,%d,%d\n",rx0, rx1, ry0, ry1, Jmax); */
  /* mexPrintf("Sampling step:%f\n",dh); */
  /* mexPrintf("nrows, ncols, Nt:%d,%d,%d\n",nrows,ncols,Nt); */

  plhs[0] = mxCreateDoubleMatrix((mwSize)nrows,(mwSize)ncols, mxREAL);
  double *A = mxGetPr(plhs[0]);
  double cst = pow(2, Jmax);
  size_t idx;

  /* mexPrintf("val=%e\n",(cst*dh*dh/M_PI/4)); */

  /* Loop for all position indexes in the range of rx X ry */
#pragma omp parallel for shared(A,Tphi) private(nx,ny,px,py,idx,xx,yy)
  for(nx=rx0; nx<=rx1; nx++){
    for(ny=ry0; ny<=ry1; ny++){
      idx = (nx-rx0) * nrows + (ny-ry0);

      A[idx] = 0;

      /* Loop inside the support of Phi for the inner product*/
      for(px=0; px<Nt; px++){
	xx = cst*(px*dh+nx) - Z[0];

	for(py=0; py<Nt; py++){
	  yy = cst*(py*dh+ny) - Z[1];
	  A[idx] += Tphi[px]*Tphi[py] * log(xx*xx+yy*yy);
	}
      }
      A[idx] *= (cst*dh*dh/M_PI/4);
      /* mexPrintf("A[%d]=%f\n",idx,A[idx]); */
    }
  }
}
