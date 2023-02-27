% RANK_CROWDDIST Computes the rank and crowding distance of the population. 
%   POP = RANK_CROWDDIST(POP,NOBJ,NVAR) computes and assigns the rank and crowding 
%   distance to all individuals of the population. POP is the decoded population 
%   matrix of NP-by-(NOBJ+NVAR+2) size, where NP is the population size, NOBJ 
%   and NVAR indicates the number of objective functions and variables of the 
%   optimization problem. For each row of POP, the first NOBJ elements are the 
%   fitness values regarding each objective function, the next NVAR numbers 
%   corresponds to the values of the solution (cut-off levels), and the last 
%   two values are the rank and the crowding distance, respectively. 
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
%   pop = rank_crowddist(pop,nobj,nvar);    % Compute rank and crowding distance   
%   plot_nondomset(pop);                    % Plot the non-dominated set  
%   
%   See also MOISLT SIMPLE_NONDOMINATED_SORT 
%   
%   
%   Reference:
%   ---------
%   K. Deb, A. Pratap, S. Agarwal and T. Meyarivan, "A fast and elitist 
%   multiobjective genetic algorithm: NSGA-II," in IEEE Transactions on 
%   Evolutionary Computation, vol. 6, no. 2, pp. 182-197, April 2002, 
%   doi: 10.1109/4235.996017.
    
% ------------------------------------------------------------------------
%   Cinvestav-IPN (Mexico)
%   RANK_CROWDDIST Version 1.0 (Matlab R2022a)
%   June 2022
%   Copyright (c) 2022, Samuel Omar Tovias Alanis 
% ------------------------------------------------------------------------
    
function pop = rank_crowddist(pop,nobj,nvar)
np = size(pop,1);      % Size of the population 
fpop = pop(:,1:nobj);  % Fitness values of the population 
fmax = max(fpop);      % Maximum values of each objective function  
fmin = min(fpop);      % Minimum values of each objective function 
index = (1:np)';       % Index of individuals
popaux = [fpop index]; % Auxiliar population 
rank = 1;              % Rank  
% Compute the rank and crowding distance 
while np ~= 0
    % Simple non-dominated sorting 
    nd = simple_nondominated_sort(np,fpop);
    id = popaux(nd,end);
    % Assign the rank to the current non-dominated set
    pop(id,nobj+nvar+1) = rank;
    % Increase the rank counter
    rank = rank + 1;
    % Compute the crowding distance of the current non-dominated set
    l = numel(nd);
    if l == 1
        pop(id,end) = inf;
    elseif l == 2
        pop(id(1),end) = inf;
        pop(id(2),end) = inf;
    else
        crowddist = zeros(l,1);
        for m=1:nobj
            [~, ind] = sort(pop(id,m));
            ind = id(ind);
            crowddist(1) = inf; 
            crowddist(end) = inf;
            for j=2:l-1
                crowddist(j,1) = crowddist(j,1) + ... 
                (pop(ind(j+1),m) - pop(ind(j-1),m)) / (fmax(m) - fmin(m));
            end
        end
        pop(ind,end) = crowddist;
    end
    % Remove the current non-dominated set from the auxiliar population 
    np = np - l;
    popaux(nd,:) = [];
    fpop = popaux(:,1:nobj); 
end 