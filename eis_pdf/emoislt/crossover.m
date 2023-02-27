% CROSSOVER Two points crossover. 
%   [BOFFSPRING,NCROSS] = CROSSOVER(BPARENTS,ID,PARAMS,NCROSS) applies a two 
%   point crossover operator to generate a binary offspring population. BPARENTS 
%   is a matrix of NP-by-CHRLEN size of the binary strings corresponding to the 
%   parents selected to be recombinated, where NP is the number of individuals 
%   in the population and CHRLEN is the length of the binary individuals (chromosomes). 
%   ID is an array of NP-by-1 size that contains the location index of the selected 
%   parents in the original population. PARAMS is a structure that contains the 
%   parameters of the optimization problem. NCROSS is a counter to save the number 
%   of crossover operations performed (optional). BOFFSPRING is a binary matrix 
%   of NP-by-CHRLEN size that contains the chromosomes of the offspring population. 
%   PARAMS.PC is the crossover probability. 
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
%                        
%   See also SELECTION MUTATION
    
% ------------------------------------------------------------------------
%   Cinvestav-IPN (Mexico)
%   CROSSOVER Version 1.0 (Matlab R2022a)
%   June 2022
%   Copyright (c) 2022, Samuel Omar Tovias Alanis 
% ------------------------------------------------------------------------
    
function [boffspring,ncross] = crossover(bparents,id,params,ncross)
np = params.np;
pc = params.pc;
nvar = params.nvar; 
nbits = params.nbits; 
boffspring = false(size(bparents));
if nargin < 4; ncross = 0; end 
for i = 1:2:np
    ind1 = id(i);
    ind2 = id(i+1); 
    while ind1 == ind2 
        ind2 = id(randperm(numel(id),1));
    end
    pos = 0;
    for j = 1:nvar
        boffspring(i,pos+1:pos+nbits(j)) = bparents(ind1,pos+1:pos+nbits(j));
        boffspring(i+1,pos+1:pos+nbits(j)) = bparents(ind2,pos+1:pos+nbits(j));
        if (rand <= pc)
            ncross = ncross + 1;
            site1 = round(rand*nbits(j));  
            site2 = round(rand*nbits(j));
            if site1 == 0; site1 = 1; end
            while (site2 == 0 || site1 == site2)
                site2 = round(rand*nbits(j));
            end
            if (site1 > site2)
                temp = site1;
                site1 = site2;
                site2 = temp;
            end
            boffspring(i,pos+site1:pos+site2) = bparents(ind2,pos+site1:pos+site2);
            boffspring(i+1,pos+site1:pos+site2) = bparents(ind1,pos+site1:pos+site2);
        end
        pos = pos + nbits(j);
    end
end 
end 