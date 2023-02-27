% TRAINLDA Training a linear discriminant analysis classification model.
function [model,exectime] = trainLDA(Xtr,Ytr)
    % Fix zero variance problem (one single element in a class)
    % Specific error: X must have more observations than the number of classes
    [Xtr,Ytr] = fixzerovarianceLDA(Xtr,Ytr);
    % Training 
    tic;
    model = fitcdiscr(Xtr,Ytr,'DiscrimType','pseudoLinear');
    exectime = toc;
end 
%****************************************************** 
% Fix issue in fitcdiscr when the data have zero variance: 
% One single element in a class 
function [X,Y] = fixzerovarianceLDA(X,Y)
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