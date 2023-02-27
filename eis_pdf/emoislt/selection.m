% SELECTION Binary tournament selection. 
%   [BPARENTS,ID] = SELECTION(POP,BPOP,PARAMS) applies the binary tournament 
%   selection operator to select the parents that will be recombinated. POP 
%   is a matrix of NP-by-(NOBJ+NVAR+2) size, where NP is the number of individuals, 
%   NOBJ is the number of optimization objectives and NVAR is the number of 
%   variables of the optimization problem. For each row of POP, the first NOBJ 
%   elements are the fitness values regarding each objective function, the next 
%   NVAR numbers corresponds to the decode values of the solution (cut-off levels), 
%   and the last two values are the rank and the crowding distance, respectively. 
%   BPOP is a matrix of NP-by-CHRLEN size, where CHRLEN is the length of the 
%   binary individual. PARAMS is a structure that contains the parameters of 
%   the optimization problem. BPARENTS is the binary population of selected 
%   parents of NP-by-CHRLEN size. ID is an array of NP-by-1 size that contains 
%   the index of the selected parents in the current population POP.                                                
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
%                                           
%   See also CROSSOVER MUTATION             
    
% ------------------------------------------------------------------------
%   Cinvestav-IPN (Mexico)
%   SELECTION Version 1.0 (Matlab R2022a)
%   June 2022
%   Copyright (c) 2022, Samuel Omar Tovias Alanis 
% ------------------------------------------------------------------------
    
function [bparents,id] = selection(pop,bpop,params)
np = params.np;
nobj = params.nobj; 
bparents = false(size(bpop));
id = zeros(np,1);
% Generate random indices 
index = randperm(np);
% Add the first index to the final
index = [index index(1)];
% Example for 4 individuals:
% [3 4 1 2 3]: 3 vs 4, 4 vs 1, 1 vs 2, 2 vs 3 
for i=1:np
    ind1 = pop(index(i),:);
    ind2 = pop(index(i+1),:);
    % Binary tournament
    idx = tournament(ind1,ind2,nobj);
    if idx == 1 
        bparents(i,:) = bpop(index(i),:);
        id(i) = index(i);
    else
        bparents(i,:) = bpop(index(i+1),:);
        id(i) = index(i+1); 
    end
end 
end 