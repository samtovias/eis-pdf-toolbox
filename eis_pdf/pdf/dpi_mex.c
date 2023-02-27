/*      
    C MEX Function.
    DPI_MEX Direct plug-in rule.  
        H = DPI_MEX(XI) computes the bandwidth with the direct plug-in rule 
        for PDF estimation using KDE. XI is an univariate continuous variable 
        of 1-by-N size. H is the bandwidth. 
        
        Compile MEX function: 
        --------------------
        mex -O dpi_mex.c -largeArrayDims; % Build an executable file for Matlab engine 
        
        Example:
        -------
        load concentric3.mat;             % Load a dataset
        xi = X(:,1);                      % Extract values of one feature
        xi = minmaxnorm(xi);              % Normalize the univariate variable
        xi = xi';                         % Transpose xi to obtain a vector of 1-by-N size. 
        h = dpi_mex(xi);                  % Compute the bandwidth with the DPI rule 
        xh = linspace(-1.5,1.5,100);      % Generate a linearly space vector 
        ph = kde(xi,xh,h);                % Compute the density estimation 
        plot(xh,ph);                      % Plot the estimation 
        t = 'KDE with DPI rule, h =';
        t = strcat(t,{' '},num2str(h)); 
        title(t);
        xlabel('x_1'); 
        ylabel('PDF'); 
        legend('PDF of x_1');
        
        References:
        ----------
        Wand, M. P., & Jones, M. C. (1994). Kernel Smoothing. Chapman & Hall/CRC. 
        pp. 72
        
        Sheather, S.J. and Jones, Chris (1991). A reliable data-based bandwidth 
        selection method for kernel density estimation. Journal of the Royal 
        Statistical Society: Series B (Statistical Methodology), 53(3) pp. 683–690.
        
        ------------------------------------------------------------------------
        Cinvestav-IPN (Mexico)
        DPI_MEX Version 1.0 (Matlab R2022a)
        June 2022
        Copyright (c) 2022, Samuel Omar Tovias Alanis 
        ------------------------------------------------------------------------
*/      
        
// Libraries
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mex.h"
#include "matrix.h" 
    
// Macros 
// Gaussian kernel K(x)
// Integral of K(x)^2 evaluated at [-inf,inf]: R(K) 
#define RK 0.282094791773878
// Integral of (x^2)*K(x) evaluated at [-inf,inf]: mu2(K) 
#define mu2K 1.0
// Sixth order derivative of K(x) evaluated at zero: K6
#define K6 -5.984134206021495
// Fourth order derivative of K(x) evaluated at zero: K4 
#define K4 1.196826841204299
// Square root of PI: SQRTPI
#define SQRTPI 1.77245385090552
    
// Function definitions
// Quicksort 
void quicksort(double array[], mwSize n) {
	double pivot, temp;
	int i, j; 
	if (n < 2) return;
	pivot = array[n / 2];
	for (i = 0, j = n - 1; ; i++, j--) {
	    while (array[i] < pivot) i++;
	    while (array[j] > pivot) j--;
	    if (i >= j) break;
	    temp = array[i];
	    array[i] = array[j];
	    array[j] = temp;
	}	
	quicksort(array, i);
	quicksort(array + i, n - i);
}   
// Gaussian kernel K(x)
double phi(double x){
    // 0.398... is the inverse of the square root of 2*PI
    return 0.398942280401433*exp(-0.5*(x*x));  
}   
// Sexth derivative of K(x) 
double D6K(double x){
	return 0.398942280401433*exp(-0.5*(x*x))*(x*x*x*x*x*x - 15.0*x*x*x*x + 45.0*x*x - 15.0);
}   
// Fourth derivative of K(x)
double D4K(double x){
	return 0.398942280401433*exp(-0.5*(x*x))*(x*x*x*x - 6.0*x*x + 3.0);
}   
// Median
double median(double array[], mwSize n){
	double med; 
	quicksort(array, n);
	return med = (n % 2 == 0) 
	? (array[n/2] + array[n/2 - 1]) / 2.0 
	: array[n/2];
}   
// Median absolute deviation (MAD)
double mad(double x[], double auxiliar[], mwSize n){
	int i; 
	double x_median; 
	x_median = median(x, n);
	for (i = 0; i < n; i++){
		auxiliar[i] = fabs(x[i] - x_median);
	}
	return median(auxiliar, n);
}   
    
// Main function 
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
    
    // Variable declaration  
    double *variable, *output, *auxiliar;       
    double madval, psi8, psi6, psi4, g1, g2, invg1, invg2, s, h;
    mwSize M, N, i, j;
    
    // Input checks  
    // Verify number of nrhs (number of right-hand-side) arguments (inputs)      
	if(nrhs != 1)
        mexErrMsgIdAndTxt( "MATLAB:dpi_mex:invalidNumInputs",
        "One input array of 1-by-N size is required.");
    // Dimensions of the input array
    M = mxGetM(prhs[0]);
    N = mxGetN(prhs[0]);
    //Check the number of rows of the input array  
    if(M != 1)
        mexErrMsgIdAndTxt( "MATLAB:dpi_mex:invalidSizeOfInputs",
                "The input array must be of 1-by-N size.");
    // MEX-pointers of the input array
    variable = mxGetPr(prhs[0]);
    
    // Allocates memory 
    // Auxiliar array of 1-by-N size
    auxiliar = (double *)malloc(sizeof(double)*N);
    
    // Step 1: Estimate psi8 
   	// MAD
  	madval = mad(variable, auxiliar, N);
    // MAD to the 9th power
  	madval = madval*madval*madval*madval*madval*madval*madval*madval*madval;
  	// psi8
    psi8 = 105.0 * (1.0 / (32.0 * SQRTPI * madval));
    
  	// Step 2: Estimate psi6
  	// Computes g1 and its inverse value 
  	g1 = pow(-2.0*K6*(1.0/(mu2K*psi8*N)),0.1111111111111111); 
    invg1 = 1 / g1; 
    s = 0;
    // Evaluation of D6K (avoiding redundancy) 
    for (i = 0; i < N-1; i++){
    	for (j = 0; j < N-i-1; j++){
    		s += D6K((variable[i+j+1] - variable[i]) * invg1); 
    	}
    }
    // Using the symmetry properties 
    s = (2.0*s) + K6*N;
    // Computes psi6
    psi6 = (1.0/(N*(N-1)))*(1.0/(g1*g1*g1*g1*g1*g1*g1))*s;
    
  	// Step 3: Estimate psi4
    // Computes g2 and its inverse value
    g2 = pow(- 2.0 * K4 * (1.0 / ( mu2K * psi6 * N) ) , 0.142857142857143 );
    invg2 = 1 / g2; 
    s = 0;
    // Evaluation of D4K (avoiding redundancy)
    for (i = 0; i < N-1; i++){
    	for (j = 0; j < N-i-1; j++){
    		s += D4K((variable[i+j+1] - variable[i]) * invg2);
    	}
    }
    // Using the symmetry properties
    s = (2.0*s) + K4*N;
    // Computes psi4
    psi4 = ( 1.0 / (N * (N-1)) ) * ( 1.0 / (g2*g2*g2*g2*g2) ) * s;
    
    // Step 4: Compute bandwidth 
    h = pow(RK*(1.0/(mu2K*mu2K*psi4*N)),0.2);
    
    // Allocate memory to the output value (bandwidth)
    plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
    output = mxGetPr(plhs[0]);
    *output = h;
    
    // Free memory 
    // MATLAB memory manager
	/* Use mxFree: for pointer variables created using mxCalloc, mxMalloc, mxRealloc,
       ...mxArrayToString, mxArrayToUTF8String, etc.
       Use mxDestroyArray: for mxArray variables */
    // Native C/C++ memory manager 
    free(auxiliar);
    
}   