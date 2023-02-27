% NICLASS Gets the number of instances in each class.
%   NI = NICLASS(Y) computes the number of instances for each class. Y is the 
%   class labels vector of N-by-1 size, where N is the number of instances of 
%   the dataset. NI is an integer vector of 1-by-C size where C is the number 
%   of classes. 
%   
%   Example: 
%   --------
%   load concentric3.mat; % Load a dataset 
%   ni = niclass(Y);      % Gets the number of instances in each class 
%   
%   See also SETUP_MOISLT GETNBITS
    
% ------------------------------------------------------------------------
%   Cinvestav-IPN (Mexico)
%   NICLASS Version 1.0 (Matlab R2022a)
%   June 2022
%   Copyright (c) 2022, Samuel Omar Tovias Alanis 
% ------------------------------------------------------------------------
    
function ni = niclass(Y)
    ni = accumarray(Y,ones(numel(Y),1),[numel(unique(Y)) 1],@sum,0)'; 
end 