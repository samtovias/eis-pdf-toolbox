% TRAINRF Training a random forest classification model.
function [model,exectime] = trainRF(Xtr,Ytr,B)
    if nargin < 3
        B = 100;
    end
    tic; 
    [Data,Params] = split_bootstraps(Xtr,Ytr,B);
    B  = Params.B;
    % Training with bootstrap samples
    Forest = cell(1,B);
    for b = 1:B
        Forest{b} = train_trees(b,Data);
    end
    exectime = toc; 
    model.Forest = Forest;
    model.Params = Params;
end 
%***********************************************************************
% Train the trees of the random forest model 
function tree = train_trees(b,Data)
    X  = Data.Train.X;
    Y  = Data.Train.Y;
    Ti = Data.Boost(:,b);
    nfea = Data.nFeats;
    Xtr = X(Ti,:);
    Ytr = Y(Ti);
    tree = fitctree(Xtr,Ytr,'MinLeaf',1,'Prune','off','Prior','uniform','NVarToSample',nfea);
end
%***********************************************************************
% Split bootstraps 
function [Data,Params] = split_bootstraps(Xtr,Ytr,B)
    % Parameters
    C = max(Ytr);       % Number of classes
    [N,D] = size(Xtr);  % Number of patterns and features
    % Parameters structure
    Params.C = C;       % Number of classes
    Params.D = D;       % Number of features
    Params.N = N;       % Number of patterns
    Params.B = B;       % Number of bootstraps
    % Bootstrap indices
    indBtr = randi(N,N,B);
    % Number of random features
    Feats = ceil(sqrt(D)); 
    % Data to train trees 
    Data.Train.X = Xtr;
    Data.Train.Y = Ytr;
    Data.Boost  = indBtr;
    Data.nFeats = Feats;
end 