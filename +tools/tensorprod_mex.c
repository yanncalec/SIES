# include "mex.h"

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
  int c, r;
  size_t ncols, nrows;
  double *X, *Y, *XY;

  /* check for proper number of arguments */
  if(nrhs!=2) {
    mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nrhs","Two inputs required.");
  }
  if(nlhs!=1) {
    mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nlhs","One output required.");
  }
  /* make sure input arguments are arrys */
  if( !(mxIsDouble(prhs[0]) || 
	mxIsComplex(prhs[0]))) {
    mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notScalar","Inputs must be a numerical arrayrs.");
  }
    
  /* create a pointer to the real data in the input matrix  */
  X = mxGetPr(prhs[0]);
  Y = mxGetPr(prhs[1]);

  /* get dimensions of inputs */
  nrows = mxGetNumberOfElements(prhs[0]);
  ncols = mxGetNumberOfElements(prhs[1]);

  /* create the output matrix */
  if(mxIsComplex(prhs[0]) || mxIsComplex(prhs[1]))
    plhs[0] = mxCreateDoubleMatrix((mwSize)nrows, (mwSize)ncols, mxREAL);
  else
    plhs[0] = mxCreateDoubleMatrix((mwSize)nrows, (mwSize)ncols, mxCOMPLEX);    

  /* get a pointer to the real data in the output matrix */
  XY = mxGetPr(plhs[0]);

  for(c=0; c<ncols; c++)
    for(r=0; r<nrows; r++)
      XY[c*nrows + r] = X[r] * Y[c];
    
  /*printf("Mex file\n");*/
}
