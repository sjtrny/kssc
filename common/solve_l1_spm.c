#include <math.h>
#include <matrix.h>
#include <mex.h>
#include <string.h>

#define sign(a) ( ( (a) < 0 )  ?  -1   : ( (a) > 0 ) )
#define MAX(x, y) (((x) > (y)) ? (x) : (y))

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    double lambda = mxGetScalar(prhs[1]);
    
    mwSize num_vals = mxGetNzmax(prhs[0]);
    
    mwSize m = mxGetM(prhs[0]);
    mwSize n = mxGetN(prhs[0]);
    
    double *pr_A = mxGetPr(prhs[0]);
    mwIndex *in_ir = mxGetIr(prhs[0]);
    mwIndex *in_jc = mxGetJc(prhs[0]);
    
    plhs[0] = mxCreateSparse(m,n,num_vals,mxREAL);
    double *pr = mxGetPr(plhs[0]);
    mwIndex *out_ir = mxGetIr(plhs[0]);
    mwIndex *out_jc = mxGetJc(plhs[0]);
    
    for (int i = 0; i < num_vals; i++)
    {
        pr[i] = sign(pr_A[i]) * MAX(fabs(pr_A[i]) - lambda, 0);
    }
    
    memcpy(out_ir, in_ir, sizeof(mwIndex) * num_vals);
    memcpy(out_jc, in_jc, sizeof(mwIndex) * (n+1));
}