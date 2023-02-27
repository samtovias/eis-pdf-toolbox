% EVALIND_GOP Evaluate the fitness of an individual. 
%   FIT = EVALIND_GOP(X,Y,XSOL,PARAMS) evaluates the fitness of an individual 
%   in a global optimization problem. X is a dataset of N-by-D size (N instances 
%   and D features). Y is the class labels vector of N-by-1 size. XSOL is a decoded
%   solution (cut-off levels) of 1-by-C size, where C is the number of classes
%   of the problem. PARAMS is a structure that contains the input parameters
%   of the optimization problem. FIT is the fitness value of the solution.
%   
%   Example:
%   -------
%   load concentric3.mat                    % Load a dataset 
%   [data,params] = setup_micro_eislt(X,Y); % Setup with the default values
%   X = data.X; mn = data.mn; mx = data.mx; % Normalized dataset
%   np = params.np; chrlen = params.chrlen; % Number of individuals and size of the chromosome 
%   nobj = params.nobj; nvar = params.nvar; % Number of objectives and optimization variables
%   bpop = logical(randi([0,1],np,chrlen)); % Randomly initialize the binary values of the individuals 
%   pop = decode(bpop, params);             % Decodes the individuals of the population
%   xsol = pop(:,nobj+1:nobj+nvar);         % Obtain the cut-off levels values
%   for i=1:np                              % Evaluate the individuals of the population
%       pop(i,1) = ...           
%       evalind_gop(X,Y,xsol(i,:),params);      
%   end                                     
%                                            
%   See also LTIS HELLINGER_DISTANCE         
    
% ------------------------------------------------------------------------
%   Cinvestav-IPN (Mexico)
%   EVALIND_GOP Version 1.0 (Matlab R2022a)
%   June 2022
%   Copyright (c) 2022, Samuel Omar Tovias Alanis 
% ------------------------------------------------------------------------
    
function f = evalind_gop(X,Y,xsol,params)
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
znorm = 1-((xsol-params.lmin)./(params.lmax-params.lmin))';
% Alpha factor for penalization 
ns = accumarray(YS,ones(numel(YS),1),[params.c 1],@sum,0);
ni = params.ni'; 
alpha = 1-abs((ns/sum(ns))-(ni/sum(ni)));
% Compute the fitness value
w = params.w; 
f = mean(mean(bsxfun(@plus,w*hd,(1-w)*znorm),2).^(alpha));
end 