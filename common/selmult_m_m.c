#include <math.h>
#include <matrix.h>
#include <mex.h>
#include <string.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    double * A = mxGetPr(prhs[0]);
    double * B = mxGetPr(prhs[1]);
    double * rows = mxGetPr(prhs[2]);
    double * cols = mxGetPr(prhs[3]);
    
    mwSize n_inds = mxGetM(prhs[2]);
    mwSize m_len_A = mxGetM(prhs[0]);
    mwSize m_len_B = mxGetM(prhs[1]);

    plhs[0] = mxCreateDoubleMatrix(n_inds, 1, mxREAL);
    double *pr = mxGetPr(plhs[0]);
    
    /* For each element we wish to calculate */
    for (int i = 0; i < n_inds; i++)
    {
        mwSize cur_row = rows[i];
        mwSize cur_col = cols[i];

        for(int k = 0; k < m_len_B; k++)
        {
            pr[i] = pr[i] + (A[k*m_len_A + cur_row-1] * B[(cur_col-1)*m_len_B + k]);
        }
    } 
    
}