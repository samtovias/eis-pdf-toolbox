% DISPLAY_STATUS Display the optimization status. 
%   DISPLAY_STATUS(PARAMS,G,TEND,FITNESS) display the status of an evolutionary 
%   optimization algorithm. PARAMS is a structure that contains the input parameters 
%   of the optimization problem. G is the number of the current generation, 
%   TEND is the record of the execution time of the current generation and FITNESS 
%   contains the fitness value. In a single objective optimization problem, 
%   FITNESS is the objective function value evaluated with the best solution 
%   of the generation G. For a multiobjective optimization problem, FITNESS is 
%   a vector of 1-by-2 size that contains the fitness values of the NSOL solution,
%   which is the nearest solution to the origin in the objective space. 
%    
%   Example:
%   --------
%   load concentric3.mat;                      % Load a dataset   
%   [data,params] = setup_moislt(X,Y,{10,10}); % Set parameters of the problem
%   X = data.X; mn = data.mn; mx = data.mx;    % Normalized dataset  
%   np = params.np; chrlen = params.chrlen;    % Number of individuals and size of the chromosome 
%   nobj = params.nobj; nvar = params.nvar;    % Number of objectives and binary variables
%   tStart = tic;                              % Take the execution time of the first generation 
%   bpop = logical(randi([0,1],np,chrlen));    % Randomly initialize the binary values of the individuals 
%   pop = decode(bpop, params);                % Decodes the individuals of the population 
%   xsol = pop(:,nobj+1:nobj+nvar);            % Obtain the decodified cut-off levels
%   for i=1:np                                 % Evaluate all the individuals of the population
%       [pop(i,1),pop(i,2)] = ...              
%       evalind(X,Y,xsol(i,:),params);         
%   end                                                                                         
%   pop = rank_crowddist(pop,nobj,nvar);       % Assign rank and crowding distance
%   fpop = pop(:,1:nobj);                      % Fitness values of the population 
%   [mdistance,ind]= min(dist([0 0],fpop'));   % Distance of NSOL to origin in the objective space
%   nsol = pop(ind,:);                         % NSOL solution 
%   nbsol = bpop(ind,:);                       % Binary string of NSOL  
%   nfit = fpop(ind,:);                        % Fitness of NSOL  
%   tEnd = toc(tStart);                        % Save the execution time of the first generation
%   display_status(params,1,tEnd,nfit);        % Display the optimization status
%                        
%   See also MOISLT SETUP_MOISLT
    
% ------------------------------------------------------------------------
%   Cinvestav-IPN (Mexico)
%   DISPLAY_STATUS Version 1.0 (Matlab R2022a)
%   June 2022
%   Copyright (c) 2022, Samuel Omar Tovias Alanis 
% ------------------------------------------------------------------------
    
function display_status(params,g,tEnd,fitness)
if nargin == 3 
    str1 = sprintf('Gen: %d/%d',g,params.gen);
elseif numel(fitness) == 1  
    str1 = sprintf('Gen: %d/%d | Fit = %.4f',g,params.gen,fitness);
elseif numel(fitness) == 2 
    str1 = sprintf('Gen: %d/%d | Fit = [%.4f,%.4f]',g,params.gen,fitness(1),fitness(2)); 
end 
if ~isempty(params.dataset) && ~isempty(params.numexp) 
status = sprintf('%s | Dataset: %s | Exp: %s | %s | Dist: %s (p = %s) | Time %f s.\n',...
         params.algorithm,params.dataset,num2str(params.numexp),str1,params.distance,params.p,tEnd);
elseif ~isempty(params.dataset) && isempty(params.numexp)
status = sprintf('%s | Dataset: %s | %s | Dist: %s (p = %s) | Time %f s.\n',...
         params.algorithm,params.dataset,str1,params.distance,params.p,tEnd);
elseif isempty(params.dataset) && isempty(params.numexp)
status = sprintf('%s | %s | Dist: %s (p = %s) | Time %.4f s.\n',...
         params.algorithm,str1,params.distance,params.p,tEnd);
end
fprintf(status);
end 