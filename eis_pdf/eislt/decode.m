% DECODE Decode the binary strings of the population.
%   POP = DECODE(POP_STRINGS,PARAMS) decodes the binary individuals of the population.  
%   POP_STRINGS is a binary matrix of NP-by-CHRLEN size, where NP is the number of       
%   individuals in the population and CHRLEN is the length of the binary individuals           
%   (chromosomes). PARAMS is a structure that contains the input parameters of the 
%   optimization problem. POP is a matrix of NP-by-(NOBJ+NVAR+L) size that contains the 
%   decoded individuals, where NOBJ and NVAR are the number of objectives and variables 
%   of the optimization problem, respectively. For each row of POP, the first NOBJ 
%   elements are the fitness values regarding each objective function, and the next NVAR 
%   values corresponds to the decode values of the solution (cut-off levels). Finally, 
%   for multiobjetive optimization problems L=2, since the last two columns of POP 
%   contains the values of the rank and the crowding distance of the NSGAII algorithm, 
%   respectively; for global optimization problems, NOBJ=1 and L=0. 
%   
%   Example:
%   -------
%   load concentric3.mat                       % Load a dataset 
%   [data,params] = setup_moislt(X,Y);         % Setup with the default values
%   X = data.X; mn = data.mn; mx = data.mx;    % Normalized dataset
%   np = params.np; chrlen = params.chrlen;    % Number of individuals and size of the chromosome 
%   bpop = logical(randi([0,1],np,chrlen));    % Randomly initialize the binary values of the individuals 
%   pop = decode(bpop, params);                % Decodes the individuals of the population 
%   
%   See also SETUP_MOISLT MINMAXNORM
    
% ------------------------------------------------------------------------
%   Cinvestav-IPN (Mexico)
%   DECODE Version 1.0 (Matlab R2022a)
%   June 2022
%   Copyright (c) 2022, Samuel Omar Tovias Alanis 
% ------------------------------------------------------------------------
    
function pop = decode(bpop,params)
l=0; if params.nobj == 2; l = 2; end 
ncol = params.nobj + params.nvar + l;
pop = zeros(params.np,ncol);
for i=1:params.np
    stop = 0;
    for j=1:params.nvar
        start = stop + 1; 
        stop = start + params.nbits(j) - 1;
        integer = sum(bpop(i,start:stop).*2.^(params.nbits(j)-1:-1:0));
        jj = params.nobj + j; 
        pop(i,jj) = params.lmin(j) + ... 
                   (params.lmax(j) - params.lmin(j)) / (2^params.nbits(j)-1) * integer;
        pop(i,jj) = round(pop(i,jj));
    end
end 
end 