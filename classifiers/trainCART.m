% TRAINCART Train a classification and regression tree model. 
function [model,exectime] = trainCART(Xtr,Ytr)
    tic;
    model = fitctree(Xtr,Ytr);
    exectime = toc; 
end 