% EMOIS-LT Evolutionary multiobjective instance selection based on linkage trees.  
%   OUT = MOISLT(X,Y,VARARGIN) obtains a set of solutions on a Pareto front 
%   representing cut-off levels of a linkage tree for each class to select a 
%   subset of instances using the LTIS algorithm. X is a dataset of N-by-D 
%   size (N instances and D dimensions). Y is a class label vector of N-by-1 
%   size. VARARGIN is a variable-length input argument list that contains the 
%   input parameters in the following order: 
%          
%   NP:       Size of the population 
%   GEN:      Number of generations
%   PC:       Crossover probability, real value in the range [0,1]
%   PM:       Mutation probability, real value in the range [0,1]
%   DISTANCE: Distance metric to build linkage trees, options: {'Minkowsky','Yang'}
%   P:        String with the value of the distance metric order, options: {'05', '2', 'INF'} (Minkowski), {'05', '1', '2', 'INF'} (Yang)
%   HOPT:     Method to compute the bandwidth used in the PDF estimation, options: {'Silverman','DPI'}
%   DATASET:  String with the dataset name (optional)
%   NUMEXP:   String or integer with the number of the current experiment (optional) 
%   
%   OUT is a structure that contains the results of the MOISLT algorithm with the following fields:  
%   
%   POP:              Final population 
%   PARAMS:           Parameters of the optimization problem
%   CURVES:           Fitness values of the nearest solution to the origin in the objective space (NSOL solution)
%   NCROSS:           Counter of the number of crossover operations performed
%   NMUT:             Counter of the number of mutation operations performed
%   MINMAXNORM_STATS: Minimum and maximum values of each dimension of X
%   EXECTIME:         Execution time of the whole process
%   DATASET:          String with the dataset name
%   NUMEXP:           String or integer with the number of the current experiment 
%   
%   Example:
%   -------
%   load concentric3.mat                % Load a dataset 
%   out = moislt(X,Y,10,10);            % Find a set of cut-off levels using MOISLT   
%   params = out.params;                % Parameters of the optimization problem  
%   nvar = params.nvar;                 % Number of variables of the optimization problem  
%   nobj = params.nobj;                 % Number of objectives of the optimization problem
%   fpop = out.pop(:,1:nobj);           % Fitness values of the population 
%   [~,ind]= min(dist([0 0],fpop'));    % Index of the nearest solution to the origin in the objective space (NSOL solution)
%   xpop = out.pop(ind,:);              % NSOL solution 
%   xsol = xpop(nobj+1:nobj+nvar);      % Cut-off levels of the NSOL solution
%   Xn = minmaxnorm(X);                 % Dataset normalization
%   [XSn,YS] = ltis(Xn,Y,xsol,params);  % Instance selection based on linkage trees (selected subset normalized) 
%   stats = out.minmaxnorm_stats;       % Minimum and maximum values of each dimension of the original dataset
%   XS = minmaxunnorm(XSn,stats);       % Unnormalization of the selected subset
%   plot_datasets(X,Y,XS,YS,'MOISLT');  % Plot the original dataset and the selected subset
%   
%   See also SETUP_MOISLT LTIS
      
% ------------------------------------------------------------------------
%   Cinvestav-IPN (Mexico)
%   MOISLT Version 1.0 (Matlab R2022a)
%   June 2022
%   Copyright (c) 2022, Samuel Omar Tovias Alanis 
% ------------------------------------------------------------------------    
    
function out = moislt(X,Y,varargin) 
% Open a parallel pool 
v = ver; parfor_flag = 0;
if any(strcmp('Parallel Computing Toolbox', {v.Name}))
    if isempty(gcp('nocreate'))
        parpool; 
    end
    poolobj = gcp('nocreate');
    parfor_flag = 1; 
end 
% Set parameters of the optimization problem 
[data,params] = setup_moislt(X,Y,varargin);
% Normalized data set in the range [-1,1]
X = data.X; mn = data.mn; mx = data.mx;
% Starts the multiobjective optimization algorithm 
% Take the execution time of the whole process 
exectime_start = tic; 
% Take the execution time of the first generation 
tStart = tic; 
% Initialization of the population  
bpop = logical(randi([0,1],params.np,params.chrlen));
% Decodes the individuals of the population
pop = decode(bpop, params);
disp('Initialization done, now performing first generation...');
% Evaluate the initial population 
xsol = pop(:,params.nobj+1:params.nobj+params.nvar);
f1 = zeros(params.np,1); 
f2 = zeros(params.np,1);
if parfor_flag == 1 
    parfor i=1:params.np
        [f1(i),f2(i)] = evalind(X,Y,xsol(i,:),params);
    end
else
    for i=1:params.np
        [f1(i),f2(i)] = evalind(X,Y,xsol(i,:),params);
    end
end
pop(:,1) = f1; 
pop(:,2) = f2; 
% Assign rank and crowding distance                                  
pop = rank_crowddist(pop,params.nobj,params.nvar);
% Find the nearest solution to the origin in the objective space (nsol)
fpop = pop(:,1:params.nobj);
[mdistance,ind]= min(dist([0 0],fpop'));
nsol = pop(ind,:); 
nbsol = bpop(ind,:);
nfit = fpop(ind,:);
% Save the execution time of the first generation
tEnd = toc(tStart);
% Display optimization status
display_status(params,1,tEnd,nfit);
% Counts of crossovers and mutations
nmut = 0;
ncross = 0;
% Save curves of nsol 
curves = struct();
curves.nfit = zeros(params.gen,2);
curves.mdistance = zeros(params.gen,1); 
curves.nfit(1,:) = nfit; 
curves.mdistance(1) = mdistance; 
% Main loop 
for g = 2:params.gen
    % Take the execution time of the gth generation 
    tStart = tic; 
    % Binary tournament selection   
    [bparents,id] = selection(pop,bpop,params);
    % Two point crossover per binary variable   
    [boffspring,ncross] = crossover(bparents,id,params,ncross); 
    % Bit flip mutation 
    [boffspring,nmut] = mutation(boffspring,params,nmut);
    % Decode the offspring population
    offspring = decode(boffspring,params);
    % Evaluate the offspring population   
    xsol = offspring(:,params.nobj+1:params.nobj+params.nvar);
    f1 = zeros(params.np,1);
    f2 = zeros(params.np,1);
    if parfor_flag == 1
        parfor i = 1:params.np
           [f1(i),f2(i)] = evalind(X,Y,xsol(i,:),params);
        end
    else
        for i = 1:params.np
            [f1(i),f2(i)] = evalind(X,Y,xsol(i,:),params);
        end
    end
    offspring(:,1) = f1; 
    offspring(:,2) = f2; 
    % Combine parent and offspring population 
    mixed_pop = [pop;offspring];
    mixed_bpop = [bpop;boffspring];
    % Assign rank and crowding distance 
    mixed_pop = rank_crowddist(mixed_pop,params.nobj,params.nvar);
    % Selects the next population
    [pop, bpop] = fill_nondominated_sort(mixed_pop,mixed_bpop,params.nobj,params.nvar);
    % Find nsol of the next population  
    fpop = pop(:,1:params.nobj);
    [distance,ind]= min(dist([0 0],fpop'));
    % Check to replace nsol
    if distance < mdistance 
        % Obtain the new nsol
        mdistance = distance; 
        nsol = pop(ind,:); 
        nbsol = bpop(ind,:);
        nfit = fpop(ind,:);
    else 
        % Maintain nsol in the next generation 
        ranks = pop(:,end-1); 
        if sum(ranks == 1) == params.np 
            replace = find(pop(:,end)~=Inf);
        else 
            replace = find(pop(:,end-1)~=1);
        end 
        ind = replace(randperm(numel(replace),1)); 
        pop(ind,:) = nsol; 
        bpop(ind,:) = nbsol; 
    end
    % Save the execution time of the gth generation 
    tEnd = toc(tStart);
    % Display optimization status
    display_status(params,g,tEnd,nfit); 
    % Save curves of nsol 
    curves.nfit(g,:) = nfit; 
    curves.mdistance(g) = mdistance; 
end
% Save the execution time of the whole process
exectime = toc(exectime_start);
% Save the results 
out.pop = pop;
out.params = params;
out.curves = curves;
out.ncross = ncross;
out.nmut = nmut;
out.minmaxnorm_stats = [mn;mx];
out.exectime = exectime;
out.dataset = params.dataset; 
out.numexp = params.numexp; 
% Delete the parallel pool 
if parfor_flag==1
    delete(poolobj);
end 
end 