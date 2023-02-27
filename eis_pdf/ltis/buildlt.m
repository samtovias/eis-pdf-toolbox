% BUILDLT Build linkage trees per class. 
%   LT = BUILDTL(X,Y,PARAMS) builds the linkage trees per class. X is a dataset 
%   of N-by-D size (N instances and D features) with C classes. Y is a class 
%   labels vector of N-by-1 size. PARAMS is a structure that contains the input
%   parameters of the optimization problem. LT is a cell array of 1-by-C size. 
%   Each element of the cell is a matrix of (N-1)-by-3 size that defines a 
%   tree of hierarchical clusters.
%   
%   Example:
%   --------
%   load concentric3.mat               % Load a dataset
%   [data,params] = setup_moislt(X,Y); % Setup with the default values 
%   X = data.X;                        % Normalized dataset
%   LT = buildlt(X,Y,params);          % Creates the linkage trees per class. 
%   
%   See also SETUP_MOISLT CREATE_LINKAGE
%   
%   
%   References:
%   ---------
%   Xu, R., & Wunsch, D. (2008). Clustering. Wiley-IEEE Press.
    
% ------------------------------------------------------------------------
%   Cinvestav-IPN (Mexico)
%   BUILDLT Version 1.0 (Matlab R2022a)
%   June 2022
%   Copyright (c) 2022, Samuel Omar Tovias Alanis 
% ------------------------------------------------------------------------
    
function LT = buildlt(X,Y,params)
LT = cell(1,params.c);           
for i = 1:params.c       
    Xi = X(Y==i,:);                 
    Zi = create_linkage(Xi,params);
    id = Zi(:,3) == 0;                             
    % Avoid cutoff values equal to zero
    Zi(id,3) = 0.1*min(Zi(~id,3)); 
    LT{i} = Zi;
end 
end 