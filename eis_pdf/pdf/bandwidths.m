% BANDWIDTHS Compute a matrix of bandwidths for PDF estimation using KDE. 
%   H = BANDWIDTHS(X,Y,PARAMS) computes the bandwidths required for the PDF 
%   estimation using the KDE algorithm. X is a dataset of N-by-D size (N instances 
%   and D features). Y is a class label vector of N-by-1 size. PARAMS is a 
%   structure that contains the input parameters of the optimization problem.
%   H is a matrix of bandwidths of C-by-D size (C classes and D features).
%   
%   Example: 
%   --------
%   load concentric3.mat                    % Load a dataset 
%   [data,params] = setup_moislt(X,Y);      % Setup with the default values
%   X = data.X; mn = data.mn; mx = data.mx; % Normalized dataset
%   H = bandwidths(X,Y,params);             % Computes the bandwidths
%   
%   See also SILVERMAN DPI  
    
% ------------------------------------------------------------------------
%   Cinvestav-IPN (Mexico)
%   BANDWIDTHS Version 1.0 (Matlab R2022a)
%   June 2022
%   Copyright (c) 2022, Samuel Omar Tovias Alanis 
% ------------------------------------------------------------------------
    
function H = bandwidths(X,Y,params)
H = zeros(params.c,params.d);
for i = 1:params.c
    Xi = X(Y == i,:);
    for j = 1:params.d
        xi = Xi(:,j)';
        if strcmp(params.hopt,'Silverman')  
            H(i,j) = silverman(xi);
        elseif strcmp(params.hopt,'DPI')
            H(i,j) = dpi_mex(xi);
            %H(i,j) = dpi(xi);
        end
        if H(i,j) == 0 || isnan(H(i,j)) || isinf(H(i,j)) 
            H(i,j) = 0.1; 
        end
    end
end 
end 