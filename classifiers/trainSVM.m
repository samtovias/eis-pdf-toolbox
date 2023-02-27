% TRAINSVM Training a support vector machine classification model 
function [model,exectime] = trainSVM(Xtr,Ytr,c,gamma)
    % LIBSVM settings 
    % -q : quiet mode (no outputs)
    % -s 0 : C-SVC (multi-class classification)
    % -t 2 : radial basis function: exp(-gamma*|u-v|^2)
    str = ['-q -s 0 -t 2 -c ' num2str(c) ' -g ' num2str(gamma)];
    % Training 
    tic;
    model = libsvmtrain(Ytr,Xtr,str);
    exectime = toc;
end 