
% Evaluate individual in SOEIS-SBE

function fitpop = evalind_sga_sbe(X,Y,xbin,params)
    % Get subset 
    XS = X(logical(xbin),:);
    YS = Y(logical(xbin));
    % Penalize individual if a class is deleted
    NC = numel(unique(YS));
    if NC ~= params.NC
        fitpop = 1;  
        return; 
    end  
    % Get PDFs matrix of subset selected 
    PS = getPDF(XS,YS,params);
    hd = zeros(params.NC,params.D);
    for c = 1:params.NC
        for d = 1:params.D
            hd(c,d) = sqrt(0.5*trapz(params.xh, (sqrt(params.P{c,d})-sqrt(PS{c,d})).^2));
        end
    end
    % Compute objective function 
    w = params.weight; 
    f1 = mean2(hd);
    f2 = size(XS,1)/size(X,1); 
    fitpop = w*f1 + (1-w)*f2;
    clear XS YS PS hd; 
end
