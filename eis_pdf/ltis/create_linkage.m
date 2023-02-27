% CREATE_LINKAGE Creates the linkage tree of a dataset. 
%   Z = CREATE_LINKAGE(X,PARAMS) creates the linkage tree of a given dataset. 
%   X is a dataset of N-by-D size (N instances and D features). PARAMS is a 
%   structure that contains the input parameters of the optimization problem.
%   Z is a matrix of (N-1)-by-3 size that defines a tree of hierarchical clusters.
%   
%   Clusters are created with the complete linkage algorithm using different 
%   distances metrics. A distance metric was setted in PARAMS. The distance 
%   metric options to build linkage trees include MINKOWSKI and the metric families 
%   proposed by YANG et al. Also, the valid values of the distance metrics order 
%   are {'05', '2', 'INF'} for MINKOWSKI, and {'05', '1', '2', 'INF'} for YANG distance.  
%   
%   Example: 
%   -------
%   load concentric3.mat               % Load a dataset
%   [data,params] = setup_moislt(X,Y); % Setup with the default values
%   X = data.X;                        % Normalized dataset
%   Xi = X(Y==1,:);                    % Extract the instances of the first class 
%   Z = create_linkage(Xi,params);     % Creates the linkage tree of the dataset Xi 
%    
%   See also COMPUTE_DISTANCE BUILDLT
%   
%   
%   References:
%   ----------
%   Xu, R., & Wunsch, D. (2008). Clustering. Wiley-IEEE Press.
%   
%   Yang, R., Jiang, Y., Mathews, S., Housworth, E. A., Hahn, M. W., & 
%   Radivojac, P. (2019). A new class of metrics for learning on 
%   real-valued and structured data. Data Mining and Knowledge Discovery, 
%   33(4), 995–1016. https://doi.org/10.1007/s10618-019-00622-6
    
% ------------------------------------------------------------------------
%   Cinvestav-IPN (Mexico)
%   CREATE_LINKAGE Version 1.0 (Matlab R2022a)
%   June 2022
%   Copyright (c) 2022, Samuel Omar Tovias Alanis 
% ------------------------------------------------------------------------
    
function Z = create_linkage(X,params)
if strcmp(params.distance,'Minkowski')
    switch params.p
        case '05' 
           Z = linkage(X,'complete',@mkwi05);   
        case '2'
           Z = linkage(X,'complete',@mkwi2);    
        case 'INF'
           Z = linkage(X,'complete',@mkwiINF); 
        otherwise
           error('%s\n%s\n','Invalid value of p tag for Minkowski metric.',...
                 'Valid values include: ''05'', ''2'', ''INF''.'); 
    end
elseif strcmp(params.distance,'Yang')
    switch params.p
        case '05' 
           Z = linkage(X,'complete',@d05metric);
        case '1'
           Z = linkage(X,'complete',@d1metric);
        case '2'
           Z = linkage(X,'complete',@d2metric);    
        case 'INF'
           Z = linkage(X,'complete',@dINFmetric); 
        otherwise
            error('%s\n%s\n','Invalid value of p tag for Yang et al. metric.',...
                  'Valid values include: ''05'', ''1'', ''2'', ''INF''.');
    end
else    
    error('Distance metric tag must be ''Minkowski'' or ''Yang''.');
end 
end 