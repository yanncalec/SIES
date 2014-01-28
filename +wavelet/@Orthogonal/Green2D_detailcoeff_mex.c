#include "mex.h"
#include <math.h>
#include <omp.h>

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
  int rx0, rx1, ry0, ry1, J, row, col, px, py, nx, ny;
  size_t Ntx, Nty, nrows, ncols;
  double x0, x1, y0, y1, dhx, dhy, xx, yy;
  double *Tphi, *Tpsi, *toto, *Z;

  /* check for proper number of arguments */
  if(nrhs!=8) {
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

  toto = mxGetPr(prhs[4]); J = (int)(*toto); 

  toto = mxGetPr(prhs[5]); 
  rx0 = (int)(toto[0]); rx1 = (int)(toto[1]); /* range in x-axis of the position index */

  toto = mxGetPr(prhs[6]); 
  ry0 = (int)(toto[0]); ry1 = (int)(toto[1]); /* range in x-axis of the position index */

  Z = mxGetPr(prhs[7]); /* singularity of the Green function */

  nrows = ry1-ry0+1; ncols = rx1-rx0+1;

  /* mexPrintf("t0,t1:%f,%f\n",t0,t1); */
  /* mexPrintf("rx, ry, Jmax:%d,%d,%d,%d,%d\n",rx0, rx1, ry0, ry1, Jmax); */
  /* mexPrintf("Sampling step:%f\n",dh); */
  /* mexPrintf("nrows, ncols, Nt:%d,%d,%d\n",nrows,ncols,Nt); */

  plhs[0] = mxCreateDoubleMatrix((mwSize)nrows, (mwSize)ncols, mxREAL);
  double *D = mxGetPr(plhs[0]);
  double cst = pow(2, J);
  size_t idx;

  /* Loop for all position indexes in the range of rx X ry */
#pragma omp parallel for shared(D,Tphi,Tpsi,x0,y0) private(nx,ny,px,py,idx,xx,yy)
  for(nx=rx0; nx<=rx1; nx++){
    for(ny=ry0; ny<=ry1; ny++){
      idx = (nx-rx0) * nrows + (ny-ry0);

      D[idx] = 0;

      /* Loop inside the support of Phi for the inner product*/
      for(px=0; px<Ntx; px++){
	xx = cst*(px*dhx+nx) - Z[0];

	for(py=0; py<Nty; py++){
	  yy = cst*(py*dhy+ny) - Z[1];
	  D[idx] += Tphi[px]*Tpsi[py] * log(xx*xx+yy*yy);
	}
      }
      D[idx] *= (cst*dhx*dhy/M_PI/4);
      /* mexPrintf("A[%d]=%f\n",idx,A[idx]); */
    }
  }
}
