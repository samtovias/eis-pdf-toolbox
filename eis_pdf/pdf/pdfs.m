% PDFS Estimate the univariate PDFs of a dataset using the KDE method.
%   PDF = PDFS(X,Y,PARAMS) estimates the univariate PDFs per class and dimension.
%   X is a dataset of N-by-D size (N instances and D features) with C classes. 
%   Y is the class labels vector of N-by-1 size. PARAMS is a structure that 
%   contains the input parameters of the optimization problem. PDF is a cell 
%   array of C-by-D size and each cell contains a vector of 1-by-100 size that 
%   defines a probability density estimation.
%   
%   Example:
%   --------
%   load concentric3.mat;               % Load a dataset
%   X = minmaxnorm(X);                  % Normalize the dataset  
%   params.d  = size(X,2);	            % Dimensions of the dataset
%   params.c = max(Y);                  % Number of classes
%   params.xh = linspace(-1.5,1.5,100); % Linearly spaced vector                 
%   params.hopt = 'Silverman';          % Method to compute bandwidth
%   params.H = bandwidths(X,Y,params);  % Computes the bandwidths 
%   PDF = pdfs(X,Y,params);             % Estimate the univariate PDFs
%   plot_pdfs(PDF,params);              % Plot PDFs
%   
%   See also KDE PLOT_PDFS  
    
% ------------------------------------------------------------------------
%   Cinvestav-IPN (Mexico)
%   PDFS Version 1.0 (Matlab R2022a)
%   June 2022
%   Copyright (c) 2022, Samuel Omar Tovias Alanis 
% ------------------------------------------------------------------------
    
function PDF = pdfs(X,Y,params)
PDF  = cell(params.c,params.d);
for i = 1:params.c
    for j = 1:params.d 
        PDF{i,j} = kde(X(Y==i,j),params.xh,params.H(i,j));
    end
end 
end 