% LTIS Instance selection based on linkage trees. 
%   [XS,YS] = ltis(X,Y,XSOL,PARAMS) selects a subset of instances using C cut-off 
%   levels and a linkage tree per class, where C is the number of classes. X 
%   is the original dataset of N-by-D size (N instances and D features), Y is 
%   the class labels vector of N-by-1 size, XSOL is the cut-off levels vector 
%   of C-by-1 or 1-by-C size. PARAMS is a structure that contains the input 
%   parameters of the optimization problem. XS is the selected subset of Ns-by-D 
%   size (Ns selected instances). YS is the class labels vector of Ns-by-1 size.   
%   
%   Example: 
%   -------- 
%   load concentric3.mat                % Load a dataset 
%   out = moislt(X,Y,10,10);            % Find a set of cut-off levels using MOISLT   
%   params = out.params;                % Parameters of the optimization problem  
%   nvar = params.nvar;                 % Number of variables of the optimization problem  
%   nobj = params.nobj;                 % Number of objectives of the optimization problem
%   fpop = out.pop(:,1:nobj);           % Fitness values of the population 
%   [~,ind]= min(dist([0 0],fpop'));    % Index of the nearest solution to the origin in the objective space (NSOL solution)
%   xpop = out.pop(ind,:);              % NSOL solution 
%   xsol = xpop(nobj+1:nobj+nvar);      % Cut-off levels of the NSOL solution
%   Xn = minmaxnorm(X);                 % Dataset normalization
%   [XSn,YS] = ltis(Xn,Y,xsol,params);  % Instance selection based on linkage trees (selected subset normalized) 
%   stats = out.minmaxnorm_stats;       % Minimum and maximum values of each dimension of the original dataset
%   XS = minmaxunnorm(XSn,stats);       % Unnormalization of the selected subset
%   plot_datasets(X,Y,XS,YS,'MOISLT');  % Plot the original dataset and the selected subset
%   
%   See also MOISLT BUILDLT 
%   
%   
%   References:
%   ---------
%   S. O. Tovias-Alanis, W. Gómez-Flores and G. Toscano-Pulido, "Instance 
%   Selection Based on Linkage Trees," 2021 18th International Conference 
%   on Electrical Engineering, Computing Science and Automatic Control (CCE), 
%   2021, pp. 1-6, doi: 10.1109/CCE53527.2021.9633116.
    
% ------------------------------------------------------------------------
%   Cinvestav-IPN (Mexico)
%   LTIS Version 1.0 (Matlab R2022a)
%   June 2022
%   Copyright (c) 2022, Samuel Omar Tovias Alanis 
% ------------------------------------------------------------------------
    
function [XS,YS] = ltis(X,Y,xsol,params)
XS = cell(1,params.c); 
YS = cell(1,params.c);
for i = 1:params.c
    Z = params.LT{i};
    th = Z(xsol(i),3);
    T = cluster(Z,'Cutoff',th,'Criterion','distance');
    mx = max(T);
    XSi = zeros(mx,params.d);
    YSi = zeros(mx,1);
    Xi = X(Y==i,:);
    for m = 1:mx
        Xm = Xi(T==m,:);
        mn = mean(Xm,1);
        id = compute_distance(Xm,mn,params);
        XSi(m,:) = Xm(id,:);
        YSi(m) = i;
    end 
    XS{i} = XSi; 
    YS{i} = YSi; 
end
XS = cat(1,XS{:});
YS = cat(1,YS{:});
end 