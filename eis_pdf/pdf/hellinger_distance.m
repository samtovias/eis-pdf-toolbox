% HELLINGER_DISTANCE Hellinger distance. 
%   HD = HELLINGER_DISTANCE(XH,PDF1,PDF2) computes the Hellinger distance between 
%   two univariate probability density functions. XH is a linearly spaced vector 
%   in range [-1.5,1.5] of 1-by-100 size and stands for the sample space. PDF1
%   and PDF2 are two estimations of probability density functions. HD is a similarity 
%   measure in the range [0,1], where 0 indicates that both PDFs are the same 
%   and 1 shows that both PDFs do not overlap in any point of XH. 
%   
%   Example:
%   -------
%   load concentric3.mat                   % Load a dataset 
%   X = minmaxnorm(X);                     % Normalize the dataset
%   xh = linspace(-1.5,1.5,100);           % Linearly spaced vector
%   h = silverman(X(:,1));                 % Compute bandwidth for variable 1  
%   pdf1 = kde(X(:,1),xh,h);               % PDF estimation of variable 1   
%   h = silverman(X(:,2));                 % Compute bandwidth for variable 2 
%   pdf2 = kde(X(:,2),xh,h);               % PDF estimation of variable 2
%   hd = hellinger_distance(xh,pdf1,pdf2); % Hellinger distance 
%   
%   See also KDE SILVERMAN
%   
%   
%   Reference:
%   ---------
%   Adele Cutler & Olga I. Cordero-Braña (1996) Minimum Hellinger Distance 
%   Estimation for Finite Mixture Models, Journal of the American Statistical 
%   Association, 91:436, 1716-1723, DOI: 10.1080/01621459.1996.10476743
    
% ------------------------------------------------------------------------
%   Cinvestav-IPN (Mexico)
%   HELLINGER_DISTANCE Version 1.0 (Matlab R2022a)
%   June 2022
%   Copyright (c) 2022, Samuel Omar Tovias Alanis 
% ------------------------------------------------------------------------
    
function hd = hellinger_distance(xh,pdf1,pdf2)
    hd = sqrt(0.5*trapz(xh,(sqrt(pdf1)-sqrt(pdf2)).^2));
end 