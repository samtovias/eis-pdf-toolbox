% EISLT Evolutionary instance selection based on linkage trees.  
%   OUT = EISLT(X,Y,VARARGIN) genetic algorithm that optimizes cut-off 
%   levels of a linkage tree by each class to select a subset of instances using 
%   the LTIS algorithm. X is a dataset of N-by-D size (N instances and D dimensions). 
%   Y is a class label vector of N-by-1 size. VARARGIN is a variable-length input 
%   argument list that contains the input parameters in the following order: 
%      
%   NP:       Size of the population
%   GEN:      Number of generations
%   W:        Weight of the objective function
%   PC:       Crossover probability, real value in the range [0,1]
%   PM:       Mutation probability, real value in the range [0,1]
%   DISTANCE: Distance metric to build linkage trees, options: {'Minkowsky','Yang'}
%   P:        String with the value of the distance metric order, options: {'05', '2', 'INF'} (Minkowski), {'05', '1', '2', 'INF'} (Yang)
%   HOPT:     Method to compute the bandwidth used in the PDF estimation, options: {'Silverman','DPI'}
%   DATASET:  String with the dataset name (optional)
%   NUMEXP:   String or integer with the number of the current experiment (optional) 
%   
%   OUT is a structure that contains the results of the EISLT algorithm with the following fields:  
%   
%   POP:              Final population 
%   PARAMS:           Parameters of the optimization problem
%   CURVES:           Fitness values of the best solution 
%   NCROSS:           Counter of the number of crossover operations performed
%   NMUT:             Counter of the number of mutation operations performed
%   MINMAXNORM_STATS: Minimum and maximum values of each dimension of X
%   EXECTIME:         Execution time of the whole process
%   DATASET:          String with the dataset name
%   NUMEXP:           String or integer with the number of the current experiment 
%   
%   Example:
%   -------
%   load concentric3.mat                     % Load a dataset 
%   out = eislt(X,Y);                        % Find a set of cut-off levels using micro_eislt   
%   params = out.params;                     % Parameters of the optimization problem  
%   fpop = out.pop(:,1:params.nobj);         % Fitness values of the final population 
%   [~,ibest] = min(fpop);                   % Index of the best solution 
%   xsol = out.pop(ibest,2:params.nvar+1);   % Cut-off levels of the best solution
%   Xn = minmaxnorm(X);                      % Dataset normalization
%   [XSn,YS] = ltis(Xn,Y,xsol,params);       % Instance selection based on linkage trees (selected subset normalized) 
%   stats = out.minmaxnorm_stats;            % Minimum and maximum values of each dimension of the original dataset
%   XS = minmaxunnorm(XSn,stats);            % Unnormalization of the selected subset
%   plot_datasets(X,Y,XS,YS,'eislt');        % Plot the original dataset and the selected subset
%   
%   See also SETUP_EISLT LTIS
%   
%   
%   Reference:
%   ---------
    
% ------------------------------------------------------------------------
%   Cinvestav-IPN (Mexico)
%   EISLT Version 1.0 (Matlab R2022a)
%   June 2022
%   Copyright (c) 2022, Samuel Omar Tovias Alanis 
% ------------------------------------------------------------------------    
    
function out = eislt(X,Y,varargin) 
% Open a parallel pool 
v = ver; parfor_flag = 0;
if any(strcmp('Parallel Computing Toolbox', {v.Name}))
    if isempty(gcp('nocreate'))
        parpool;
        off_parpool = 1;
    else
        off_parpool = 0; 
    end
    poolobj = gcp('nocreate');
    parfor_flag = 1; 
end 
% Set parameters of the optimization problem 
[data,params] = setup_eislt(X,Y,varargin);
% Normalized data set in the range [-1,1]
X = data.X; mn = data.mn; mx = data.mx;
% Starts the micro evolutionary optimization algorithm 
% Take the execution time of the whole process 
exectime_start = tic; 
% Initialization of the population  
bpop = logical(randi([0,1],params.np,params.chrlen));
% Decodes the individuals of the population
pop = decode(bpop, params);
disp('Initialization of EISLT done...');
% Evaluate the initial population 
xsol = pop(:,params.nobj+1:params.nobj+params.nvar);
fpop = zeros(params.np,1);
if parfor_flag
    parfor i=1:params.np
        fpop(i) = evalind_gop(X,Y,xsol(i,:),params);
    end
else
    for i=1:params.np
        fpop(i) = evalind_gop(X,Y,xsol(i,:),params);
    end
end
pop(:,1) = fpop;     
% Current best solution  
[fbest,ibest] = min(fpop);
bbest = bpop(ibest,:); 
xbest = pop(ibest,2:params.nvar+1);
% Counts of crossovers
ncross = 0;
% Counts of mutations
nmut = 0;
% Curves
curves = zeros(params.gen+1,2);
curves(1,1) = fbest;
curves(1,2) = mean(fpop);
% Main loop  
for g = 1:params.gen   
    % Take the execution time of each generation 
    tStart = tic;  
    % Binary tournament selection  
    [bparents,id] = selection_gop(pop,bpop,params);
    % Crossover 
    [bpop,ncross] = crossover_gop(bparents,id,params,ncross);
    % Mutation 
    [bpop,nmut] = mutation_gop(bpop,params,nmut);
    % Decode 
    pop = decode(bpop,params);
    % Evaluate the current population 
    xsol = pop(:,params.nobj+1:params.nobj+params.nvar);
    fpop = zeros(params.np,1);
    if parfor_flag 
        parfor i=1:params.np
            fpop(i) = evalind_gop(X,Y,xsol(i,:),params);
        end
    else
        for i=1:params.np
            fpop(i) = evalind_gop(X,Y,xsol(i,:),params);
        end
    end
    pop(:,1) = fpop;    
    % Best solution of the offspring   
    [fit,ind] = min(fpop);
    % Update the best solution   
    if (fit < fbest)    
        fbest = fit; 
        bbest = bpop(ind,:);
        xbest = pop(ind,2:params.nvar+1);
    else
        % Elitism  
        ind = randi([1,params.np],1,1);
        fpop(ind) = fbest; 
        bpop(ind,:) = bbest;
        pop(ind,2:params.nvar+1) = xbest;
        pop(ind,1) = fbest; 
    end
    % Save the execution time of the gth generation 
    tEnd = toc(tStart);
    % Display optimization status
    display_status(params,g,tEnd,fbest);  
    % Save curves
    curves(g+1,1) = fbest;
    curves(g+1,2) = mean(fpop);
end
% Save the execution time of the whole process
exectime = toc(exectime_start);
% Save the results 
out.pop = pop;
out.params = params;
out.curves = curves;
out.ncross = ncross;
out.minmaxnorm_stats = [mn;mx];
out.exectime = exectime;
out.dataset = params.dataset; 
out.numexp = params.numexp; 
% Delete the parallel pool 
if parfor_flag && off_parpool
    delete(poolobj);
end 
end 