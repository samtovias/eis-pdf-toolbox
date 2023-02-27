% MINMAXUNNORM Min-max unnormalization. 
%   X = MINMAXUNNORM(XN,STATS) un-normalizes the dataset XN of N-by-D size 
%   (N instances and D features) previously normalized in the range [-1,1]. 
%   STATS is a matrix of 2-by-D size that contains D column vectors of 2-by-1 
%   size with the information of the minimum and maximum values of each dimension
%   of the original dataset. X is the un-normalized dataset of N-by-D size.     
%   
%   Example:
%   -------
%   load concentric3.mat          % Load a dataset 
%   [XN,dmn,dmx] = minmaxnorm(X); % Normalize the dataset
%   stats = [dmn;dmx];            % Minimum and maximum values of the dataset
%   X2 = minmaxunnorm(XN,stats);  % Un-normalize the dataset 
%   
%   See also SETUP_MOISLT MINMAXNORM
%   
%   
%   Reference:
%   ---------
%   K. L. Priddy, P. E. Keller, Artificial Neural Networks: An Introduction.
%   Bellingham, WA: SPIE-The Int. Soc. Optical Eng., 2005.
    
% ------------------------------------------------------------------------
%   Cinvestav-IPN (Mexico)
%   MINMAXUNNORM Version 1.0 (Matlab R2022a)
%   June 2022
%   Copyright (c) 2022, Samuel Omar Tovias Alanis 
% ------------------------------------------------------------------------
    
function X = minmaxunnorm(X,stats)
dmn = stats(1,:);
dmx = stats(2,:);
N = size(X,1);
mx = repmat(dmx,N,1);
mn = repmat(dmn,N,1);
X = 0.5*(X+1).*(mx-mn) + mn;
end 