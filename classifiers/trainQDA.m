% TRAINQDA Training a quadratic discriminant analysis classification model.
function [model,exectime] = trainQDA(Xtr,Ytr)
    % Fix zero variance problem (one single element in a class)
    [Xtr,Ytr] = fixzerovarianceQDA(Xtr,Ytr);
    % Training
    tic;
    model = fitcdiscr(Xtr,Ytr,'DiscrimType','pseudoquadratic');
    exectime = toc;
end
%****************************************************** 
% Fix issue in fitcdiscr when the data have zero variance: 
% One single element in a class 
function [X,Y] = fixzerovarianceQDA(X,Y)
    C = numel(unique(Y)); 
    uniqueY = unique(Y);
    for j = 1:C
        if size(X(uniqueY(j)==Y,:),1) == 1
            sample = X(uniqueY(j)==Y,:);
            X = [X;sample]; 
            Y = [Y;uniqueY(j)];
        end    
    end
end