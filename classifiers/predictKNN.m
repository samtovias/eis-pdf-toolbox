% PREDICTKNN Predict the class label using a KNN model. 
function Yp = predictKNN(model,Xtt) 
    % Search KNN of Xtt 
    [idx,~] = knnsearch(model.kdtree,Xtt,'k',model.k);
    if size(idx,1) > 1 
        Yp = mode(model.Ytr(idx),2);
    else 
        Yp = mode(model.Ytr(idx));
    end
end 