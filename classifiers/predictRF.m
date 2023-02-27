% PREDICTRF Predict the class label using a RF model.
function Yp = predictRF(model,Xtt)
    B = model.Params.B;
    C = model.Params.C;
    N = size(Xtt,1);
    % Classify out-of-bag
    Ypb = zeros(N,B);
    for i = 1:B
        Ypb(:,i) = predict(model.Forest{i},Xtt);
    end
    % Majority vote
    V = zeros(N,C);
    for i = 1:N
        Ypi = Ypb(i,:);
        votes = accumarray(Ypi',ones(numel(Ypi),1),[C 1],@sum,0);
        V(i,:) = votes;
    end
    [~,Yp] = max(V,[],2);
end 