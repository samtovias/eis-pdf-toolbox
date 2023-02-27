% MINMAXNORM Min-max normalization. 
%   XN = MINMAXNORM(X) normalizes the dataset X of N-by-D size (N instances 
%   and D features) using the minimun and maximum values in each dimension. 
%   XN is the normalized dataset with scaled features in the range [-1,1].
%   
%   [XN,DMN,DMX] = MINMAXNORM(X) returns the minimum and maximum values for 
%   each dimension in the DMN and DMX vectors, each of 1-by-D size.
%   
%   XN = MINMAXNORM(X,[DMN;DMN]) normalizes the dataset X using the minimum 
%   and maximum values given in DMN and DMX.     
%   
%   Example:
%   -------
%   load concentric3.mat          % Load a dataset  
%   [Xn,dmn,dmx] = minmaxnorm(X); % Normalize the dataset X 
%   
%   See also SETUP_MOISLT MINMAXUNNORM
%   
%   
%   Reference:
%   ---------
%   K. L. Priddy, P. E. Keller, Artificial Neural Networks: An Introduction.
%   Bellingham, WA: SPIE-The Int. Soc. Optical Eng., 2005.
    
% ------------------------------------------------------------------------
%   Cinvestav-IPN (Mexico)
%   MINMAXNORM Version 1.0 (Matlab R2022a)
%   June 2022
%   Copyright (c) 2022, Samuel Omar Tovias Alanis 
% ------------------------------------------------------------------------
    
function [X,dmn,dmx] = minmaxnorm(X,stats)
if nargin == 1
    dmx = max(X,[],1);
    dmn = min(X,[],1);
elseif nargin == 2
    dmn = stats(1,:);
    dmx = stats(2,:);
end
N = size(X,1);
ind = dmx~=dmn;
mx = repmat(dmx,N,1);
mn = repmat(dmn,N,1);
X_aux = (2.*(X(:,ind) - mn(:,ind))./(mx(:,ind)-mn(:,ind)))-1;
X(:,ind) = X_aux;
% Avoid NaN values when variables have no range
zeros_vectors = zeros(N,sum(~ind));  
X(:,~ind) = zeros_vectors;
end 