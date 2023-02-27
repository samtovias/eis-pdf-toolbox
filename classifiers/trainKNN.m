% TRAINKNN Fit a kd-tree model for k-nearest neighbors classification. 
function [model,exectime] = trainKNN(Xtr,Ytr,k) 
    if nargin < 3 
        k = 3; 
    end 
    model = struct();
    tic; 
    kdtree = createns(Xtr,'nsmethod','kdtree'); 
    exectime = toc;
    model.kdtree = kdtree; 
    model.Ytr = Ytr; 
    model.k = k; 
end 