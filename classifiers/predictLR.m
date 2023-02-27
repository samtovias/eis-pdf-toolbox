% PREDICTLR Predict the class label using a LR model.
function Yp = predictLR(model,Xtt)
    ph = mnrval(model,Xtt);
    [~,Yp] = max(ph,[],2);
end 