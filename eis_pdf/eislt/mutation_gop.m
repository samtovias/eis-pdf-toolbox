% MUTATION Bit-flip mutation for global optimization. 
%   [BPOP,NMUT] = MUTATION(BPOP,PARAMS,NMUT) applies a bit-flip mutation to
%   the offspring population. BPOP is a matrix of NP-by-CHRLEN size that contains 
%   the binary strings (chromosomes) of the offspring, where NP is the number 
%   of individuals in the population and CHRLEN is the length of the chromosomes. 
%   PARAMS is a structure that contains the parameters of the optimization problem. 
%   NMUT is a counter to save the number of mutations performed (optional). 
%   PARAMS.PM is the mutation probability. 
%     
%   Example:
%   -------
%   load concentric3.mat                    % Load a dataset 
%   [data,params] = setup_moislt(X,Y);      % Setup with the default values
%   X = data.X; mn = data.mn; mx = data.mx; % Normalized dataset
%   np = params.np; chrlen = params.chrlen; % Number of individuals and size of the chromosome 
%   nobj = params.nobj; nvar = params.nvar; % Number of objectives and problem variables
%   bpop = logical(randi([0,1],np,chrlen)); % Randomly initialize the binary values of the individuals 
%   pop = decode(bpop, params);             % Decodes the individuals of the population 
%   xsol = pop(:,nobj+1:nobj+nvar);         % Obtain the decodified cut-off levels
%   for i=1:np                              % Evaluate all the individuals of the population
%       [pop(i,1),pop(i,2)] = ...           
%       evalind(X,Y,xsol(i,:),params);      
%   end                                                        
%   pop = rank_crowddist(pop,nobj,nvar);    % Assign rank and crowding distance
%   [b1,id] = selection(pop,bpop,params);   % Binary tournament selection
%   b2 = crossover(b1,id,params);           % Two-point crossover operator                     
%   b2 = mutation(b2,params);               % Bit-flip mutation
%    
%   See also SELECTION CROSSOVER
    
% ------------------------------------------------------------------------
%   Cinvestav-IPN (Mexico)
%   MUTATION Version 1.0 (Matlab R2022a)
%   June 2022
%   Copyright (c) 2022, Samuel Omar Tovias Alanis 
% ------------------------------------------------------------------------
    
function [bpop,nmut] = mutation_gop(bpop,params,nmut)
np = params.np; 
nvar = params.nvar; 
nbits = params.nbits; 
pm = params.pm; 
if nargin < 3; nmut = 0; end
for i = 1:np
    ind = bpop(i,:); 
    count = 0;
    for j = 1:nvar
        for k = 1:nbits(j)
            if (rand <= pm)
                nmut = nmut + 1;
                if (ind(1,count+k) == 0)
                    ind(1,count+k) = 1;
                else
                    ind(1,count+k) = 0;
                end
            end
        end
        count = count + nbits(j);
    end
    bpop(i,:) = ind;
end 
end 