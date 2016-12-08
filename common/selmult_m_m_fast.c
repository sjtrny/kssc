#include <math.h>
#include <matrix.h>
#include <mex.h>
#include <string.h>
#include <blas.h>

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
    
    mwSize cur_row;
    mwSize cur_col;
    
    mwSignedIndex one = 1;
    mwSignedIndex row_gap = (mwSignedIndex)m_len_A;
    
    mwSignedIndex dot_length = (mwSignedIndex)m_len_B;
    
    /* For each element we wish to calculate */
    for (int i = 0; i < n_inds; i++)
    {
        cur_row = rows[i];
        cur_col = cols[i];

        pr[i] = ddot(&dot_length, &A[cur_row-1], &row_gap, &B[m_len_B * (cur_col-1)], &one);
    } 
    
}