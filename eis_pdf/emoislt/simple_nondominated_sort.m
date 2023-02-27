% SIMPLE_NONDOMINATED_SORT Simple non-dominated sorting (naive approach). 
%   P = SIMPLE_NONDOMINATED_SORT(NP,FPOP) applies a non-dominated sorting to 
%   the individuals of the population. NP is the number of individuals. FPOP 
%   is an array of fitness values of NP-by-NOBJ size, where NOBJ is the number 
%   of objective functions in the multioptimization problem. P is an array that 
%   contains the index of the non-dominated solutions in the current population. 
%   
%   Example: 
%   -------
%   load concentric3.mat                    % Load a dataset 
%   [data,params] = setup_moislt(X,Y);      % Setup with the default values
%   X = data.X; mn = data.mn; mx = data.mx; % Normalized dataset
%   np = params.np; chrlen = params.chrlen; % Number of individuals and size of the chromosome 
%   nobj = params.nobj; nvar = params.nvar; % Number of objectives and optimization variables 
%   bpop = logical(randi([0,1],np,chrlen)); % Randomly initialize the binary values of the individuals 
%   pop = decode(bpop, params);             % Decodes the individuals of the population 
%   xsol = pop(:,nobj+1:nobj+nvar);         % Obtain the decodified cut-off levels
%   for i=1:np                              % Evaluate all the individuals of the population
%       [pop(i,1),pop(i,2)] = ...           
%       evalind(X,Y,xsol(i,:),params);      
%   end   
%   fpop = pop(:,1:nobj);                  % Extract the fitness values 
%   p = simple_nondominated_sort(np,fpop); % Obtain the non-dominated set
%   plot_nondomset(pop,p);                 % Plot the non-dominated set
%   
%   See also RANK_CROWDDIST DOMINANCE
    
% ------------------------------------------------------------------------
%   Cinvestav-IPN (Mexico)
%   SIMPLE_NONDOMINATED_SORT Version 1.0 (Matlab R2022a)
%   June 2022
%   Copyright (c) 2022, Samuel Omar Tovias Alanis 
% ------------------------------------------------------------------------
    
function p = simple_nondominated_sort(np,fpop)
p = [];
dom = 0;
for i=1:np
    for j=1:np
        if j ~= i
            % Check the non-dominance status of the solutions j and i 
            dom = dominance(fpop(j,:),fpop(i,:));
            % If j dominates i 
            if dom == 1
                break
            end
        end
    end
    % If none j dominates the ith solution, then i is added to non-dominated set 
    if dom ~= 1
        p = cat(2,p,i);
    end
    dom = 0;
end 
end 