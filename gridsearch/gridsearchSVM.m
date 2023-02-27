% GRIDSEARCHSVM Grid search to fit the hyperparameters of an SVM classification model.
function [C,g,grid] = gridsearchSVM(X,Y,params)
    if nargin < 3 
        K = 10;      
        name = '';
    else
        K = params.K; 
        name = params.name; 
    end
    % Cross-validation and grid-search for SVM tuning
    kf = crossvalind('KFold',Y,K); 
    Cr = 2.^(-5:15);
    Gr = 2.^(-15:3); 
    nC = numel(Cr);
    nG = numel(Gr);
    grid = zeros(nC,nG,K);
    for k = 1:K
        grid(:,:,k) = grid_search(X,Y,k,kf,Cr,Gr,K,name);
    end
    m = mean(grid,3);              % Mean of the folds 
    [~,k] = max(m(:));             % Best mean performance value  
    [k1,k2] = ind2sub([nC nG],k);  % Indices of results 
    C = Cr(k1);                    % C 
    g = Gr(k2);                    % Gamma 
end
%**************************************************************************
% Grid search for the kth fold.
function grid = grid_search(X,Y,k,kf,Cr,Gr,K,name)
    % Split data into training and validation 
    idx = kf == k;   
    Xtr = X(~idx,:); 
    Ytr = Y(~idx,:);
    Xvd = X(idx,:);  
    Yvd = Y(idx,:); 
    % Grid-search
    nC = numel(Cr);
    nG = numel(Gr);
    grid = zeros(nC,nG);
    cont = 0; 
    for i = 1:nC 
        for j = 1:nG  
            cont = cont + 1; 
            str = ['-s 0 -t 2 -c ' num2str(Cr(i)) ' -g ' num2str(Gr(j)) ' -q'];
            tic;                                             
            model = libsvmtrain(Ytr,Xtr,str);                % Training 
            exectime = toc;                                  
            [~,perf,~] = libsvmpredict(Yvd,Xvd,model,'-q');  % Validation 
            grid(i,j) = perf(1);                             % Accuracy 
            display_status(name,k,K,cont,nC,nG,exectime);    
        end
    end    
end
%**************************************************************************
% Display grid search status
function display_status(name,k,K,cont,nC,nG,exectime)
    fprintf(['Grid search: SVM | Dataset: %s | Fold: %d/%d | Combination:' ...
             ' %d/%d | Time: %.4f s.\n'],name,k,K,cont,nC*nG,exectime); 
end