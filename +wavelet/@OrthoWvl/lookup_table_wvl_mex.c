# include "mex.h"
# include <math.h>
/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
  int c, r;
  size_t ncols, nrows, nzmax;
  int Kx;
  double sl, sr, dx, percent_sparse, Mx;
  double *X, *N, *tt, *Tx, *dTx, *sr_V, *sr_dV;
  mwIndex *irs_V, *jcs_V, *irs_dV, *jcs_dV, k;

  /* check for proper number of arguments */
  if(nrhs!=7) {
    mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nrhs","Seven inputs required.");
  }
  if(nlhs!=2) {
    mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nlhs","Two outputs required.");
  }
    
  /* create pointers to the real data in the input matrices  */
  X = mxGetPr(prhs[0]);
  N = mxGetPr(prhs[1]);
  tt = mxGetPr(prhs[2]); sl = *tt; 
  /* mexPrintf("Lower bound:%f\n",sl); */
  tt = mxGetPr(prhs[3]); sr = *tt; 
  /* mexPrintf("Upper bound:%f\n",sr); */
  tt = mxGetPr(prhs[4]); dx = *tt; 
  /* mexPrintf("Sampling step:%f\n",dx); */

  Tx = mxGetPr(prhs[5]);
  dTx = mxGetPr(prhs[6]);
  
  /* get dimensions of inputs */
  nrows = mxGetNumberOfElements(prhs[0]);
  ncols = mxGetNumberOfElements(prhs[1]);
  size_t nTx = mxGetNumberOfElements(prhs[5]);
  size_t ndTx = mxGetNumberOfElements(prhs[6]);

  percent_sparse = 0.1;
  nzmax=(mwSize)ceil((double)ncols*(double)nrows*percent_sparse);

  plhs[0] = mxCreateSparse(nrows,ncols,nzmax,mxREAL);
  plhs[1] = mxCreateSparse(nrows,ncols,nzmax,mxREAL);

  /* mexPrintf("nTx=%d, ndTx=%d, nzmax=%d\n", nTx, ndTx, nzmax); */

  sr_V  = mxGetPr(plhs[0]);
  irs_V = mxGetIr(plhs[0]);
  jcs_V = mxGetJc(plhs[0]);

  sr_dV  = mxGetPr(plhs[1]);
  irs_dV = mxGetIr(plhs[1]);
  jcs_dV = mxGetJc(plhs[1]);

  k = 0;

  for(c=0; c<ncols; c++){
    jcs_V[c] = k;
    jcs_dV[c] = k;

    /* mexPrintf("nrows=%d, ncols=%d, nzmax=%d, col=%d, k=%d\n", nrows, ncols, nzmax, c, k); */

    for(r=0; r<nrows; r++){
      if (k>=nzmax){
	/* mexPrintf("Reallocating memory\n"); */

	mwSize oldnzmax = nzmax;
	percent_sparse += 0.1;
	nzmax = (mwSize)ceil((double)ncols*(double)nrows*percent_sparse);
                    
	/* make sure nzmax increases atleast by 1 */
	if (oldnzmax == nzmax)
	  nzmax++;
                    
	mxSetNzmax(plhs[0], nzmax);
	mxSetPr(plhs[0], mxRealloc(sr_V, nzmax*sizeof(double)));
	mxSetPr(plhs[1], mxRealloc(sr_dV, nzmax*sizeof(double)));
	mxSetIr(plhs[0], mxRealloc(irs_V, nzmax*sizeof(mwIndex)));
	mxSetIr(plhs[1], mxRealloc(irs_dV, nzmax*sizeof(mwIndex)));
                    
	sr_V  = mxGetPr(plhs[0]);
	irs_V = mxGetIr(plhs[0]);
	sr_dV  = mxGetPr(plhs[1]);
	irs_dV = mxGetIr(plhs[1]);
      }

      Mx = X[r] + N[c];
      /* mexPrintf("Mx=%f\n", Mx); */
  
      if((Mx >= sl) && (Mx <= sr)){

	Kx = (int)(floor((Mx - sl)/ dx) + 1);
	/* Kx = min(max((size_t)(floor((Mx - sl)/ dx) + 1), 0), nTx); */
	if ((Kx>=0) && (Kx<nTx)) {
	  /* if ((Kx>=nTx) || (Kx<0)) */
	  /*   mexPrintf("Index range error\n"); */

	  sr_V[k] = Tx[Kx];
	  irs_V[k] = r;

	  sr_dV[k] = dTx[Kx];
	  irs_dV[k] = r;

	  k++;
	}
      }
    }
  }
  jcs_V[ncols] = k;
  jcs_dV[ncols] = k;
}
