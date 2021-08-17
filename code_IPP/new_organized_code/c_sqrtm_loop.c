#include <math.h>
#include "mex.h"

void c_sqrt_loop(int n1, double *T_re, double *T_im, double *sqrt_T_diag_re, double *sqrt_T_diag_im, double *R_re, double *R_im)
{
    int i, j, k;
    double s_re, s_im, a, b, c, d, e;
    for (i = 0; i < n1*n1; i ++)
    {   
		R_re[i] = R_im[i] = 0.0;
		
	}
		
    for (j = 0; j < n1; j ++)
    {   
	    //R[j][j] = csqrt(T[j][j]);
		R_re[j*n1+j]=sqrt_T_diag_re[j];
		R_im[j*n1+j]=sqrt_T_diag_im[j];
		
        for (i = j-1; i >= 0; i--)
        {   s_re = 0.0;
 		    s_im = 0.0;
		    if (i+1<=j-1) 
		    { 
               for (k=i+1; k<=j-1; k++)			
		        { //s += R[i][k]*R[k][j]; 
			      s_re += R_re[k*n1+i]*R_re[j*n1+k]-R_im[k*n1+i]*R_im[j*n1+k];
				  s_im += R_re[k*n1+i]*R_im[j*n1+k]+R_im[k*n1+i]*R_re[j*n1+k];
			    }
			} 
			
		    //R[i][j] = (T[i][j] - s)/(R[i][i]+R[j][j]);
			a = T_re[j*n1 + i] - s_re;
			b = T_im[j*n1 + i] - s_im;
			
			c = R_re[i*n1 + i] + R_re[j*n1 + j];
			d = R_im[i*n1 + i] + R_im[j*n1 + j];
			e = c*c+d*d;
			R_re[j*n1+i] = (a*c+b*d)/e;
			R_im[j*n1+i] = (b*c-a*d)/e;
        }		
    }
}


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	double *T_re, *T_im, *sqrt_T_diag_re, *sqrt_T_diag_im, *R_re, *R_im;
    int n1, n2, n3, n4;
    char msg[1000];
    
    if (nrhs != 2) mexErrMsgTxt("Two inputs required.");
    if (nlhs != 1) mexErrMsgTxt("One output required.");
       
    T_re = mxGetPr(prhs[0]); // get the real part.
	T_im = mxGetPi(prhs[0]); // get the imaginery part.	
    n1 = (int)mxGetM(prhs[0]);  // get the num of row.
    n2 = (int)mxGetN(prhs[0]);  // get the num of cols.    
    if (n1 != n2) mexErrMsgTxt("Input T must be a square matrix.");
	
	sqrt_T_diag_re = mxGetPr(prhs[1]); 
	sqrt_T_diag_im = mxGetPi(prhs[1]); // get the imaginery part.	
    n3 = (int)mxGetM(prhs[1]);  // get the num of row.
    n4 = (int)mxGetN(prhs[1]);  // get the num of cols.    
    if (n3 != n1 || n4 != 1) mexErrMsgTxt("Input sqrt_T_diag must be a vector of length n.");
    
    mexPrintf("Check inputs:");
    sprintf(msg, "dim of T = %d,\n", n1);
    mexPrintf(msg);
    
    plhs[0] = mxCreateDoubleMatrix(n1, n1, mxCOMPLEX);
    R_re = mxGetPr(plhs[0]);
	R_im = mxGetPi(plhs[0]);
        
    c_sqrt_loop(n1, T_re, T_im, sqrt_T_diag_re, sqrt_T_diag_im, R_re, R_im);
	// Matlab call: R=c_sqrtm_loop(T,sqrt_T_diag);
}

/*
compile in matlab:
mex c_sqrtm_loop.c

*/