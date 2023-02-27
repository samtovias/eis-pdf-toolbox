% PREDICTMLP Predict the class label using a MLP model.
function Yp = predictMLP(model,Xtt)
    % Prediction 
    out = model(Xtt');  
    Yp = vec2ind(out);
    Yp = Yp'; 
end 