% PREDICTRBFN Predict the class label using a RBFN model. 
function Yp = predictRBFN(model,Xtt)
    % Avoid NaN's values 
    Xtt(isnan(Xtt)) = 0;
    % Decode model
    W = model.W; 
    C = model.C;
    SO = model.SO;
    % Predict class of test data 
    [nvl,~] = size(Xtt);
    DV  = dist(Xtt,C');
    S = SO(ones(nvl,1),:);
    phi2 = exp(-(DV.^2)./(2*S.^2)); 
    phi2 = cat(2,ones(nvl,1),phi2);
    [~,Ynet] = max(W'*phi2',[],1);
    Yp = Ynet'; 
end