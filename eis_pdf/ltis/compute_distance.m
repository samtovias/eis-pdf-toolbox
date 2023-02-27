% COMPUTE_DISTANCE Obtain the nearest instance to a centroid.
%   ID = COMPUTE_DISTANCE(XM,MN,PARAMS) obtains the closest instance to the 
%   cluster centroid. XM is a set of instances that form a cluster. MN is the 
%   centroid of XM. PARAMS is a structure that contains the input parameters 
%   of the optimization problem. ID corresponds to the index of the nearest 
%   instance in XM to MN.
%   
%   Distances are computed using different metrics. A distance metric was setted 
%   in PARAMS. The distance metric options include MINKOWSKI and the metric 
%   families proposed by YANG et al. Also, the valid values of the distance 
%   metrics order are {'05', '2', 'INF'} for MINKOWSKI, and {'05', '1', '2', 'INF'} 
%   for YANG distance.
%   
%   Example: 
%   --------
%   load concentric3.mat                % Load a dataset  
%   [data,params] = setup_moislt(X,Y);  % Setup with the default values
%   X = data.X;                         % Normalized dataset
%   mn = mean(X,1);                     % Compute the centroid of X 
%   id = compute_distance(X,mn,params); % Obtains the nearest instance to the centroid  
%   
%   See also CREATE_LINKAGE SETUP_MOISLT
%   
%   
%   Reference: 
%   ---------
%   Yang, R., Jiang, Y., Mathews, S., Housworth, E. A., Hahn, M. W., & 
%   Radivojac, P. (2019). A new class of metrics for learning on 
%   real-valued and structured data. Data Mining and Knowledge Discovery, 
%   33(4), 995–1016. https://doi.org/10.1007/s10618-019-00622-6
    
% ------------------------------------------------------------------------
%   Cinvestav-IPN (Mexico)
%   COMPUTE_DISTANCE Version 1.0 (Matlab R2022a)
%   June 2022
%   Copyright (c) 2022, Samuel Omar Tovias Alanis 
% ------------------------------------------------------------------------
    
function id = compute_distance(Xm,mn,params)
if strcmp(params.distance,'Minkowski')
    switch params.p
        case '05' 
           [~,id] = pdist2(Xm,mn,@mkwi05,'Smallest',1);
        case '2'
           [~,id] = pdist2(Xm,mn,@mkwi2,'Smallest',1);    
        case 'INF'
           [~,id] = pdist2(Xm,mn,@mkwiINF,'Smallest',1); 
        otherwise
           error('%s\n%s\n','Invalid value of p tag for Minkowski metric.',...
                 'Valid values include: ''05'', ''2'', ''INF''.'); 
    end
elseif strcmp(params.distance,'Yang')
    switch params.p
        case '05' 
           [~,id] = pdist2(Xm,mn,@d05metric,'Smallest',1);
        case '1'
           [~,id] = pdist2(Xm,mn,@d1metric,'Smallest',1);
        case '2'
           [~,id] = pdist2(Xm,mn,@d2metric,'Smallest',1);   
        case 'INF'
           [~,id] = pdist2(Xm,mn,@dINFmetric,'Smallest',1);
        otherwise
            error('%s\n%s\n','Invalid value of p tag for Yang et al. metric.',...
                  'Valid values include: ''05'', ''1'', ''2'', ''INF''.');
    end
else    
     error('Distance metric tag must be ''Minkowski'' or ''Yang''.');
end 
end 