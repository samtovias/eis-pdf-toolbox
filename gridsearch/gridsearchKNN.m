% GRIDSEARCHKNN Grid search to fit the k value of an KNN classification model.
function [k,grid] = gridsearchKNN(X,Y,params)
    if nargin < 3 
        K = 10;      
        name = '';
    else
        K = params.K; 
        name = params.name; 
    end
    % Cross-validation and grid-search for SVM tuning
    kf = crossvalind('KFold',Y,K); 
    knum = [1 3 5 7 9];
    grid = zeros(numel(knum),K);
    for k = 1:K
        grid(:,k) = grid_search(X,Y,k,K,kf,knum,name);
    end
    m = mean(grid,2);  % Mean of the folds
    [~,ind] = max(m);  % Best mean performance value  
    k = knum(ind);     % Best k-value    
end 
%**************************************************************************
% Grid search for the kth fold.
function grid = grid_search(X,Y,k,K,kf,knum,name)
    % Split data into training and validation 
    idx = kf == k;   
    Xtr = X(~idx,:); 
    Ytr = Y(~idx,:);
    Xvd = X(idx,:);  
    Yvd = Y(idx,:); 
    % Grid-search
    grid = zeros(numel(knum),1); 
    for i = 1:numel(knum)
        knn = knum(i); 
        [model,exectime] = trainKNN(Xtr,Ytr,knn);  % Training  
        Yp = predictKNN(model,Xvd);                % Validation 
        acc = sum(Yp==Yvd)/numel(Yvd);             % Accuracy
        grid(i) = acc;                               
        display_status(name,k,K,knn,exectime); 
    end  
end 
%**************************************************************************
% Display grid search status.
function display_status(name,k,K,knn,exectime)
    fprintf(['Grid search: KNN | Dataset: %s | Fold: %d/%d | k-value:' ...
             ' %d | Time: %.4f s.\n'],name,k,K,knn,exectime); 
end 