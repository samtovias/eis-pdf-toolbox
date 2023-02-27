% TRAINMLP Training a multilayer-perceptron network classification model.
function [model,exectime] = trainMLP(Xtr,Ytr,h)
    % Training size 
    N = size(Xtr,1);
    % Hidden nodes is square root of N 
    if nargin < 3
        h = max(3,round(sqrt(N)));
    end
    hiddenSizes = h; 
    % Training function: scaled conjugate gradient backpropagation 
    % Optimization function: cross-entropy
    trainFcn = 'trainscg';
    performFcn = 'crossentropy';
    net = patternnet(hiddenSizes,trainFcn,performFcn);
    net.trainParam.epochs = 1000;
    net.trainParam.showWindow = false;
    Ttr = zeros(max(Ytr),size(Xtr,1)); 
    for ii = 1:size(Xtr,1)
        Ttr(Ytr(ii),ii) = 1;
    end
    % Training  
    tic;
    model = train(net,Xtr',Ttr); 
    exectime = toc;
end 