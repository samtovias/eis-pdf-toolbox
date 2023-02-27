% EVALIND Evaluate the fitness of an individual in the objective functions. 
%   [F1,F2] = EVALIND(X,Y,XSOL,PARAMS) evaluates the fitness of an individual 
%   in the multiobjective optimization problem. X is a dataset of N-by-D size 
%   (N instances and D features). Y is the class labels vector of N-by-1 size.  
%   XSOL is a decoded solution (cut-off levels) of 1-by-C size, where C is the 
%   number of classes of the problem. PARAMS is a structure that contains the 
%   input parameters of the optimization problem. F1 and F2 are the fitness 
%   values of XSOL regarding each objective function.
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
%   xsol = pop(:,nobj+1:nobj+nvar);         % Obtain the cut-off levels values
%   for i=1:np                              % Evaluate the individuals of the population
%       [pop(i,1),pop(i,2)] = ...           
%       evalind(X,Y,xsol(i,:),params);      
%   end                                     
%                                            
%   See also LTIS HELLINGER_DISTANCE         
    
% ------------------------------------------------------------------------
%   Cinvestav-IPN (Mexico)
%   EVALIND Version 1.0 (Matlab R2022a)
%   June 2022
%   Copyright (c) 2022, Samuel Omar Tovias Alanis 
% ------------------------------------------------------------------------
    
function [f1,f2] = evalind(X,Y,xsol,params)
% Instance selection based on linkage trees 
[XS,YS] = ltis(X,Y,xsol,params);    
% Probability density functions of the selected subset
PDF = pdfs(XS,YS,params);
% Hellinger distance 
hd = zeros(params.c,params.d);
for i = 1:params.c
    for j = 1:params.d
         hd(i,j) = hellinger_distance(params.xh,params.PDF{i,j},PDF{i,j});
    end
end 
% Normalized cut-off levels 
rnorm = 1-((xsol-params.lmin)./(params.lmax-params.lmin))';
% Alpha factor for penalization 
ns = accumarray(YS,ones(numel(YS),1),[params.c 1],@sum,0);
ni = params.ni'; 
alpha = abs((ns/sum(ns))-(ni/sum(ni)));
% Compute the fitness values
f1 = mean(mean(hd,2).^(1-alpha));
f2 = mean(rnorm.^(1-alpha));
end 