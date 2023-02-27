% TRAINRBFN Training a radial basis function network classification model.
function [model,exectime] = trainRBFN(Xtr,Ytr,hn)
    if nargin < 3 
        hn = round(sqrt(size(Xtr,1)));
    end
    % Avoid NaN's values
    Xtr(isnan(Xtr))=0;
    c = numel(unique((Ytr)));
    [ntr,~] = size(Xtr);
    tic; 
    % Unsupervised step to get centers and radii
    [C,~] = mkmeans(Xtr,hn);
    p = 2;
    if p == ntr 
        p = 1;
    end    
    DX = dist(Xtr,C');
    DC = dist(C,C');
    DC = sort(DC,2,'ascend');
    SO = mean(DC(:,2:p+1),2)';
    if sum(SO==0) ~= 0
        SO(SO==0) = eps;
    end    
    S = SO(ones(ntr,1),:);
    % Responses of hidden nodes 
    phi = exp(-(DX.^2)./(2*S.^2)); 
    phi = cat(2,ones(ntr,1),phi);
    % Compute weights 
    base = meshgrid(1:c,1:ntr);
    Trg  = Ytr(:,ones(1,c))== base;
    W = pinv(phi)*Trg;
    exectime = toc; 
    % Get model: weights, centers and radii 
    model.W = W;
    model.C = C; 
    model.SO = SO; 
end
%*******************************************************************
% Modified kmeans algorithm to avoid empty clusters
function  [C,idx] = mkmeans(X,k)
    % Take randomly k points from X
    s = size(X,1);
    k = min(s,k);
    p = randperm(s); 
    zk = X(p(1:k),:);
    % Small threshold
    val = 1e-6; 
    % Main loop
    j = 0; band = 0;
    while ~band && j < 1000
        % Compute square euclidean distance ||xj-zk||^2
        djk = dist(zk,X').^2;
        zkold = zk;
        % Index of points by cluster
        [~,idx1] = sort(djk,1,'ascend'); 
        idx2 = idx1(1,:); 
        % Points by cluster
        n = accumarray(idx2',ones(numel(idx2),1),[k,1]);
        % Update centers 
        for i=1:k
            xjck = sum(X(idx2==i,:),1);
            zk(i,:) = (1/(n(i)+1))*(xjck + zkold(i,:));
        end   
        % Check centers stability
        band = sum(sum((zk-zkold) < val)) == numel(zk);
        j = j + 1;
    end    
    % Outputs: centers and indices 
    C = zk; 
    idx = idx2;  
end