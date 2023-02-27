% DOMINANCE Check the dominance status.
%   DOM = DOMINANCE(A,B) checks the dominance status of two solutions. A and
%   B are vectors of 1-by-NOBJ size with the fitness values of the solutions. 
%   NOBJ is the number of objectives of the optimization problem. DOM is an 
%   integer value that takes the value of 1 if the solution A dominates B, -1 
%   if B dominates A, and 0 if both solutions are non-dominated. 
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
%   fpop = pop(:,1:nobj);                   % Extract the fitness values 
%   dom = dominance(fpop(1,:),fpop(2,:));   % Checking the dominance status of the first two individuals
%   
%   See also SIMPLE_NONDOMINATED_SORT RANK_CROWDDIST 

% ------------------------------------------------------------------------
%   Cinvestav-IPN (Mexico)
%   DOMINANCE Version 1.0 (Matlab R2022a)
%   June 2022
%   Copyright (c) 2022, Samuel Omar Tovias Alanis 
% ------------------------------------------------------------------------
    
function dom = dominance(a, b)
% Flags initialized to zero
flag1 = 0;
flag2 = 0;
% If solution a has at least 1 objective value better than b
if sum(a < b)
    flag1 = 1; 
end
% If solution b has at least 1 objective value better than a
if sum(a > b)
    flag2 = 1; 
end
% Dominance evaluation  
if (flag1 == 1) && (flag2 == 0)
    % Solution a dominates b 
    dom = 1;
elseif (flag1 == 0) && (flag2 == 1)
    % Solution b dominates a
    dom = -1;
else 
    % Both solutions are non-dominated
    dom = 0;
end
end             