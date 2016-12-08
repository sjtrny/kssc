#include <math.h>
#include <matrix.h>
#include <mex.h>
#include <string.h>

#define sign(a) ( ( (a) < 0 )  ?  -1   : ( (a) > 0 ) )
#define MAX(x, y) (((x) > (y)) ? (x) : (y))

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    double lambda = mxGetScalar(prhs[1]);
    
    mwSize num_vals = mxGetM(prhs[0]);
    
    double *pr_A = mxGetPr(prhs[0]);
//     mwIndex *ir_A = mxGetIr(prhs[0]);
//     mwIndex *jc_A = mxGetJc(prhs[0]);

//     plhs[0] = mxCreateSparse(mxGetM(prhs[0]), mxGetN(prhs[0]), num_vals, mxREAL);
// 
//     double *pr = mxGetPr(plhs[0]);
// //     mwIndex *ir = mxGetIr(plhs[0]);
// //     mwIndex *jc = mxGetJc(plhs[0]);
//     
//     
// //     double *pr = mxMalloc(num_vals * sizeof(double));
// //     memcpy(pr, pr_A, num_vals * sizeof(double));
// //     mxSetPr(plhs[0], pr);
//     
//     mwIndex *ir = mxMalloc(num_vals * sizeof(mwIndex));
//     memcpy(ir, ir_A, num_vals * sizeof(mwIndex));
//     mxSetIr(plhs[0], ir);
//     
//     mwIndex *jc = mxMalloc(num_vals * sizeof(mwIndex));
//     memcpy(jc, jc_A, num_vals * sizeof(mwIndex));
//     mxSetJc(plhs[0], jc);
//     
//     /* For each element we wish to calculate */
//     for (int i = 0; i < num_vals; i++)
//     {
//         pr[i] = sign(pr_A[i]) * MAX(fabs(pr_A[i]) - lambda, 0);
//     }
    
    plhs[0] = mxCreateDoubleMatrix(num_vals, 1, mxREAL);
    double *pr = mxGetPr(plhs[0]);
    for (int i = 0; i < num_vals; i++)
    {
        pr[i] = sign(pr_A[i]) * MAX(fabs(pr_A[i]) - lambda, 0);
    }
    
    
}