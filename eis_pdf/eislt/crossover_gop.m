% CROSSOVER_GOP Two points crossover for global optimization. 
%   [BOFFSPRING,NCROSS] = CROSSOVER_GOP(BPARENTS,ID,PARAMS,NCROSS) applies a two 
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
    
%                        
%   See also MICRO_EISLT SELECTION_GOP
    
% ------------------------------------------------------------------------
%   Cinvestav-IPN (Mexico)
%   CROSSOVER_GOP Version 1.0 (Matlab R2022a)
%   June 2022
%   Copyright (c) 2022, Samuel Omar Tovias Alanis 
% ------------------------------------------------------------------------
    
function [boffspring,ncross] = crossover_gop(bparents,id,params,ncross)
np = params.np;
pc = params.pc;
chrlen = params.chrlen; 
boffspring = false(size(bparents));
if nargin < 4; ncross = 0; end 
for i = 1:2:np
    ind1 = id(i);
    ind2 = id(i+1); 
    while ind1 == ind2 
        ind2 = id(randperm(numel(id),1));
    end 
    boffspring(i,:) = bparents(ind1,:);
    boffspring(i+1,:) = bparents(ind2,:);
    if (rand <= pc)
        ncross = ncross + 1;
        site1 = randi([2 chrlen-2],1);
        site2 = randi([site1+1 chrlen-1],1);
        boffspring(i,site1:site2) = bparents(ind2,site1:site2);
        boffspring(i+1,site1:site2) = bparents(ind1,site1:site2);
    end
end 
end 