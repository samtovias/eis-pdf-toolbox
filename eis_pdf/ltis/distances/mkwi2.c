/*       
    C MEX Function.
    MKWI2 Minkowski unnormalized distance of order 2.  
        MD2 = MKWI2(X1,X2) computes the Minkowski unnormalized distance between 
        two single vectors. X1 and X2 are single vectors of 1-by-D size. MD2 
        is the Minkowski unnormalized distance of order 2. 
        
        Compile MEX function: 
        --------------------
        mex -O mkwi2.c -largeArrayDims; % Build an executable file for Matlab engine
        
        Example:
        -------
        load concentric3.mat;           % Load a dataset
        x1 = X(1,:);                    % Extract the first instance 
        x2 = X(2,:);                    % Extract the second instance
        md2 = mkwi2(x1,x2);             % Minkowski unnormalized distance
        
        References:
        ----------
        Yang, R., Jiang, Y., Mathews, S., Housworth, E. A., Hahn, M. W., & 
        Radivojac, P. (2019). A new class of metrics for learning on 
        real-valued and structured data. Data Mining and Knowledge Discovery, 
        33(4), 995–1016. https://doi.org/10.1007/s10618-019-00622-6
        
        ------------------------------------------------------------------------
        Cinvestav-IPN (Mexico)
        MKWI2 Version 1.0 (Matlab R2022a)
        June 2022
        Copyright (c) 2022, Wilfrido Gomez Flores 
        ------------------------------------------------------------------------
*/      
    
// Libraries 
#include "math.h"
#include "mex.h"
    
// Function definition 
void dist_xy(double* x, double* y, double* dist, size_t rows, size_t cols){
    mwSize i, j;
    double t;
    double xj, yij;
    for (i=0; i<rows; i++) {
        t = 0;
        for (j=0; j<cols; j++) {
            xj  = x[j];
            yij = y[i+j*rows];
            t += pow(fabs(xj-yij),2);
        }
        dist[i] = sqrt(t); 
    }
}   
    
// Main function
void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]){
    
    #if MX_HAS_INTERLEAVED_COMPLEX
        mxDouble *Xptr, *Yptr, *Dptr;
    #else
        double *Xptr, *Yptr, *Dptr;
    #endif
    
    size_t rowsX, rowsY, colsX, colsY;
    mwSize i, j;
    
    if (nrhs != 2) {
        mexErrMsgIdAndTxt("MATLAB:d1Nmetric:invalidNumInputs", "Two input arguments required.");
    } else if (nlhs > 1) {
        mexErrMsgIdAndTxt("MATLAB:d1Nmetric:maxlhs", "Too many output arguments.");
    }
    
    /* Check the dimensions of x and y */
    rowsX = mxGetM(prhs[0]);
    colsX = mxGetN(prhs[0]);
    rowsY = mxGetM(prhs[1]);
    colsY = mxGetN(prhs[1]);
    
    if (colsX != colsY)
        mexErrMsgIdAndTxt("MATLAB:d1Nmetric:invalidD", "x and y must have the same dimension.");
    if (rowsX > 1)
        mexErrMsgIdAndTxt("MATLAB:d1Nmetric:invalidX", "x must be the a single vector.");
    
    /* Create a matrix for the return argument */
    plhs[0] = mxCreateDoubleMatrix((mwSize)rowsY, (mwSize)rowsX, mxREAL);
    
    /* Assign pointers */
    #if MX_HAS_INTERLEAVED_COMPLEX
        Dptr = mxGetDoubles(plhs[0]);
        Xptr = mxGetDoubles(prhs[0]);
        Yptr = mxGetDoubles(prhs[1]);
    #else
        Dptr  = mxGetPr(plhs[0]);
        Xptr  = mxGetPr(prhs[0]);
        Yptr  = mxGetPr(prhs[1]);
    #endif
    
    /* Calculate distance */
    dist_xy(Xptr,Yptr,Dptr,rowsY,colsY);
    
}   