% FILL_NONDOMINATED_SORT Generate the next population.
%   [POP,BPOP] = FILL_NONDOMINATED_SORT(MIXED_POP,MIXED_BPOP,NOBJ,NVAR) generates 
%   the new population by selecting the best NP individuals in terms of rank 
%   and crowding distance. MIXED_POP is a matrix of (NPx2)-by-(NOBJ-NVAR+2) 
%   size that contains the parents and offspring population combined. MIXED_BPOP 
%   is a matrix of (NPx2)-by-CHRLEN size with the binary strings of the individuals 
%   in MIXED_POP. NP is the population size, NOBJ is the number of objectives, 
%   NVAR is the number of optimization variables and CHRLEN is the length of 
%   the binary individuals (chromosomes). POP is a matrix of NP-by-(NOBJ+NVAR+2) 
%   size that contains the selected individuals and BPOP is a matrix of NP-by-CHRLEN 
%   size with the chromosomes of the selected individuals. 
%   
%   Example:
%   -------
%   load concentric3.mat                      % Load a dataset 
%   [data,params] = setup_moislt(X,Y);        % Setup MOISLT with default values
%   X = data.X; mn = data.mn; mx = data.mx;   % Normalized data set in the range [-1,1]
%   np = params.np; chrlen = params.chrlen;   % Number of individuals and size of the chromosome 
%   nobj = params.nobj; nvar = params.nvar;   % Number of objectives and problem variables 
%   bpop = logical(randi([0,1],np,chrlen));   % Randomly initialize the binary values of the individuals 
%   pop = decode(bpop, params);               % Decodes the individuals of the population 
%   xsol = pop(:,nobj+1:nobj+nvar);           % Obtain the decodified cut-off levels
%   for i=1:np                                % Evaluate all the individuals of the population
%       [pop(i,1),pop(i,2)] = ...           
%       evalind(X,Y,xsol(i,:),params);      
%   end                                                        
%   pop = rank_crowddist(pop,nobj,nvar);      % Assign rank and crowding distance
%   [b1,id] = selection(pop,bpop,params);     % Binary tournament selection
%   b2 = crossover(b1,id,params);             % Two-point crossover operator                     
%   boffspring = mutation(b2,params);         % Bit-flip mutation
%   offspring = decode(boffspring,params);    % Decodes the offspring population
%   xsol = offspring(:,nobj+1:nobj+nvar);     % Obtain the decodified cut-off levels
%   for i=1:np                                % Evaluate all the individuals of the offspring population 
%       [offspring(i,1),offspring(i,2)] = ...           
%       evalind(X,Y,xsol(i,:),params);        
%   end                                                                               
%   mpop = [pop;offspring];                   % Combine parent and offspring population 
%   mbpop = [bpop;boffspring];                % Combine binary parent and offspring population
%   mpop = rank_crowddist(mpop,nobj,nvar);    % Assign rank and crowding distance 
%   [pop,bpop] = fill_nondominated_sort(...   % Generate the next population 
%                mpop,mbpop,nobj,nvar);
%   
%   See also MOISLT RANK_CROWDDIST 

% ------------------------------------------------------------------------
%   Cinvestav-IPN (Mexico)
%   FILL_NONDOMINATED_SORT Version 1.0 (Matlab R2022a)
%   June 2022
%   Copyright (c) 2022, Samuel Omar Tovias Alanis 
% ------------------------------------------------------------------------
    
function [pop, bpop] = fill_nondominated_sort(mixed_pop,mixed_bpop,nobj,nvar)
m = size(mixed_pop,1);
fitsize = m/2;
id_parent = []; 
front = 1; 
front_array = []; 
[sorted_ranks, id_ranks] = sort(mixed_pop(:,nobj+nvar+1));
while numel(id_parent) < fitsize
    aux = id_ranks(sorted_ranks==front); 
    front_array = cat(1,front_array,aux);
    size_check = numel(id_parent) + numel(front_array);
    if size_check == fitsize
        id_parent = cat(1,id_parent,front_array);
        break
    elseif size_check < fitsize
        id_parent = cat(1,id_parent,front_array);
        front = front + 1;
        front_array = [];
    else
        % Sort crowding distance and select the individuals to fill population  
        miss_size = fitsize - (size_check - numel(front_array));
        rank = nobj + nvar + 1;
        crowddist = mixed_pop(front_array,rank+1);
        [~,id] = sort(crowddist,'descend');
        completes_pop = front_array(id(1:miss_size)); 
        id_parent = cat(1,id_parent,completes_pop);
    end
end 
pop = mixed_pop(id_parent,:);
bpop = mixed_bpop(id_parent,:);
end 