% COMPILE_METRICS Compile the MEX functions of the distance metrics. 
%   COMPILE_METRICS Builds the executable files for Matlab engine of the distance
%   metrics used to build linkage trees and compute distances for the LTIS algorithm. 
%    
%   Example: 
%   --------
%   compile_metrics; % Compile the MEX functions of the distance metrics
%   
%   See also LTIS MOISLT
    
% ------------------------------------------------------------------------
%   Cinvestav-IPN (Mexico)
%   COMPILE_METRICS Version 1.0 (Matlab R2022a)
%   June 2022
%   Copyright (c) 2022, Wilfrido Gomez Flores 
% ------------------------------------------------------------------------
    
function compile_metrics() 
    mex -O d05metric.c -largeArrayDims 
    mex -O d1metric.c -largeArrayDims
    mex -O d2metric.c -largeArrayDims
    mex -O dINFmetric.c -largeArrayDims
    mex -O mkwi05.c -largeArrayDims
    mex -O mkwi2.c -largeArrayDims
    mex -O mkwiINF.c -largeArrayDims 
end 