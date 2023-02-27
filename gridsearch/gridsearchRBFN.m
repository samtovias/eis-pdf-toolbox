% GRIDSEARCHRBFN Grid search to fit the hyperparameters of an RBFN classification model.
function [hn,grid] = gridsearchRBFN(X,Y,params)
    if nargin < 3 
        K = 10;      
        name = '';
    else
        K = params.K; 
        name = params.name; 
    end
    % Cross-validation and grid-search for RBFN tuning
    kf = crossvalind('KFold',Y,K); 
    n = size(X,1);
    hn = 3:3*round(sqrt(n));
    grid = zeros(numel(hn),K);
    for k = 1:K
        grid(:,k) = grid_search(X,Y,k,K,kf,hn,name);
    end
    m = mean(grid,2);  % Mean of the folds
    [~,ind] = max(m);  % Best mean performance value         
    hn = hn(ind);      % Number of hidden nodes   
end 
%**************************************************************************
% Grid search for the kth fold.
function grid = grid_search(X,Y,k,K,kf,hn,name)
    % Split data into training and validation 
    idx = kf == k;   
    Xtr = X(~idx,:); 
    Ytr = Y(~idx,:);
    Xvd = X(idx,:);  
    Yvd = Y(idx,:); 
    % Grid-search
    grid = zeros(numel(hn),1); 
    for i = 1:numel(hn)
        h = hn(i); 
        [model,exectime] = trainRBFN(Xtr,Ytr,h);  % Training  
        Yp = predictRBFN(model,Xvd);              % Validation 
        acc = sum(Yp==Yvd)/numel(Yvd);            % Accuracy
        grid(i) = acc;                               
        display_status(name,k,K,h,hn(end),exectime);    
    end  
end 
%**************************************************************************
% Display grid search status. 
function display_status(name,k,K,h,hn,exectime) 
    fprintf(['Grid search: RBFN | Dataset: %s | Fold: %d/%d | Hidden nodes:' ...
             ' %d/%d | Time: %.4f s.\n'],name,k,K,h,hn,exectime); 
end 