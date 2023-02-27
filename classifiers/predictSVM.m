% PREDICTSVM Predict the class label using a SVM model.
function Yp = predictSVM(model,Xtt)
    l = ones(size(Xtt,1),1); % Not used to predict
    Yp = libsvmpredict(l,Xtt,model,'-q');
end 