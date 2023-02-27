
function [f1,f2] = evalind_sbe_capsule(X,Y,xbin,params)
    
    % Get subset 
    XS = X(logical(xbin),:);
    YS = Y(logical(xbin));
    
    % Penalize individual if a class is deleted
    NCS = numel(unique(YS));
    if NCS ~= params.NC
        f1 = 1; 
        f2 = 1;  
        return; 
    end  
    
    PS = getPDF(XS,YS,params);
    hd = zeros(params.NC,params.D);
    for c = 1:params.NC
        for d = 1:params.D
            hd(c,d) = sqrt(0.5*trapz(params.xh, (sqrt(params.P{c,d})-sqrt(PS{c,d})).^2));
        end
    end
    
    f1 = mean2(hd);
    f2 = size(XS,1)/size(X,1); 
    clear XS YS PS hd; 
    
end
