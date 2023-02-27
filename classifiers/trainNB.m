% TRAINNB Train a Naive-Bayes classification model. 
function [model,exectime] = trainNB(Xtr,Ytr)
    % Training 
    [Xtr,Ytr] = fixzerovariance(Xtr,Ytr);
    tic;
    model = fitcnb(Xtr,Ytr);
    exectime = toc;
end
%****************************************************** 
% Fix zero variance problems in Naive Bayes
function [X,Y] = fixzerovariance(X,Y)
    % Fix zero variance problem type 1 (one single element in a class)
    [X,Y] = fixzerovarianceNB_1(X,Y);
    % Fix zero variance problem type 2 (a dimension with zero variance)
    X = fixzerovarianceforNB_2(X,Y);
end
%****************************************************** 
% Fix issue in fitcnb when the data have zero variance: 
% Problem type 1 (one single element in a class)        
function [X,Y] = fixzerovarianceNB_1(X,Y)               
    C = numel(unique(Y));                               
    uniqueY = unique(Y);                                 
    for j = 1:C                                         
        if size(X(uniqueY(j)==Y,:),1) == 1                      
            sample = X(uniqueY(j)==Y,:) + eps;
            X = [X;sample]; 
            Y = [Y;uniqueY(j)];
        end    
    end
end
%*******************************************************
% Fix issue in fitcnb when the data have zero variance: 
% Problem type 2 (a dimension with zero variance) 
function X = fixzerovarianceforNB_2(X,Y)
    [~,D] = size(X);
    C = numel(unique(Y)); 
    uniqueY = unique(Y);
    for i = 1:D
        for j = 1:C
            if var(X(uniqueY(j)==Y,i)) == 0
                datax = X(uniqueY(j)==Y,i);  
                datax(1) = datax(1) + eps; 
                X(uniqueY(j)==Y,i) = datax; 
            end    
        end
    end
end 