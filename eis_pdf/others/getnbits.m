% GETNBITS Get the number of bits to encode each variable in a genetic algorithm.  
%   NBITS = GETNBITS(PARAMS) computes the number of bits required to encodes
%   each variable of an optimization problem. PARAMS is a structure that contains 
%   some required parameters. NBITS is an integer vector of 1-by-NVAR size, 
%   where NVAR is the number of variables of the problem.
%   
%   Example:
%   -------- 
%   load concentric3.mat;           % Load a dataset 
%   params = struct();              % Initialize the structure
%   params.c = numel(unique(Y));    % Number of classes
%   params.ni = niclass(Y);         % Number of instances per class
%   params.nvar = params.c;         % Number of optimization variables
%   params.lmin = ones(1,params.c); % Lower limit of each variable
%   params.lmax = params.ni-1;      % Upper limit of each variable 
%   params.precision = 0;           % Precision of the variables 
%   nbits = getnbits(params);       % Number of bits per variable  
%   
%   See also SETUP_MOISLT NICLASS
    
% ------------------------------------------------------------------------
%   Cinvestav-IPN (Mexico)
%   GETNBITS Version 1.0 (Matlab R2022a)
%   June 2022
%   Copyright (c) 2022, Samuel Omar Tovias Alanis 
% ------------------------------------------------------------------------
    
function nbits = getnbits(params)
nbits = zeros(1,params.nvar);
for i = 1:params.nvar 
    nbits(i) = floor(log2((params.lmax(i)-params.lmin(i))*(10^params.precision)))+1;
end 
end 