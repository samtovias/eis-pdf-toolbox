% TOURNAMENT Binary tournament.
%   ID = TOURNAMENT(IND1,IND2,NOBJ) performs a binary tournament between two
%   individuals. IND1 and IND2 are two vectors of 1-by-(NOBJ+NVAR+2) size that
%   contain NOBJ fitness values, NVAR cut-off levels, the rank and the crowding 
%   distance of the individuals, respectively. NOBJ is the number of objectives 
%   of the optimization problem. ID indicates the winner individual. 
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
%   ind1 = pop(1,:);                        % Selects the first individual 
%   ind2 = pop(2,:);                        % Selects the second individual
%   id = tournament(ind1,ind2,nobj);        % Applies the binary tournament
%                                           
%   See also DOMINANCE SELECTION
    
% ------------------------------------------------------------------------
%   Cinvestav-IPN (Mexico)
%   TOURNAMENT Version 1.0 (Matlab R2022a)
%   June 2022
%   Copyright (c) 2022, Samuel Omar Tovias Alanis 
% ------------------------------------------------------------------------
    
function id = tournament(ind1,ind2,nobj)
flag = dominance(ind1(1:nobj),ind2(1:nobj));
% Check dominance 
if (flag == 1)
    id = 1; 
elseif (flag == -1)
    id = 2; 
% Check crowding distance 
elseif (ind1(1,end) > ind2(1,end))
    id = 1; 
elseif (ind2(1,end) > ind1(1,end))
    id = 2; 
% Random selection    
elseif (rand <= 0.5)
    id = 1; 
else
    id = 2; 
end 
end 