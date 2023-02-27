
% NSGA-II for Instance Selection with SBE encoding

function out = nsga2_sbe(X,Y,nsga2,dataset,numexp)

% Normalize in the range [-1,1] 
[X,mn,mx] = minmaxnorm(X);               

% SBE PARAMS                      
params.N  = size(X,1);               % Cardinality of dataset    
params.D  = size(X,2);	             % Dimensionality of dataset
params.NC = max(Y);                  % Number of classes
params.xh = linspace(-1.5,1.5,100);  % PDF grid              
params.h  = h_estimate(X,Y,params);  % Compute bandwidth for PDF estimate 
params.P  = getPDF(X,Y,params);      % Obtain original PDFs                           

% NSGA-II params
popsize = nsga2.popsize; 
ngen = nsga2.ngen; 
nobj = nsga2.nobj; 
ncon = nsga2.ncon; 
nreal = nsga2.nreal; 
nsga2.nbin = 1; 
nbits = zeros(1,numel(nsga2.nbin));
for i = 1:numel(nbits)
    nsga2.nbits(i) = params.N;
    nsga2.min_binvar(i) = 0;
    nsga2.max_binvar(i) = 1;
end
nsga2.bitlength = sum(nsga2.nbits);
nbin = nsga2.nbin; 
nbits = nsga2.nbits;
%min_binvar = nsga2.min_binvar; % (Not necessary in SBE)
%max_binvar = nsga2.max_binvar; 
pcross_bin = nsga2.pcross_bin;
pmut_bin = nsga2.pmut_bin;
bitlength = nsga2.bitlength;

% Global constants 
global INF EPS E PI
INF = 1.0e14;
EPS = 1.0e-14;
E = 2.71828182845905;
PI = 3.14159265358979;

% Counts of crossovers and mutations
nbinmut = 0;
nbincross = 0;

% Initialization 
%ncolumn = nobj+ncon+nreal+nbin+1+1+1; % cons_viol+rank+crowd_dist (original)
ncolumn = nobj+ncon+nreal+nbits+1+1+1; % cons_viol+rank+crowd_dist (SBE)
parent_pop = zeros(popsize,ncolumn);
% parent_strings = round(rand(popsize,bitlength)); (original)
typeinit = 'sparse';
parent_strings = initpop_nsga2_sbe(popsize,bitlength,typeinit); % (SBE)
parent_pop = decode_sbe(parent_pop, parent_strings, nobj, ncon, nreal);
disp(' Initialization done, now performing first generation');

% First generation
tStart = tic;

% Evaluation of the population 
% Instance selection with simple binary encoding 
%********************************* 
%parent_pop = evalind_sbe(parent_pop, popsize, nobj, ncon, nreal, nbin,...
%                                     X, Y, params); % (original)
%parent_pop = evalind_sbe(parent_pop, popsize, nobj, ncon, nreal, nbits,...
%                                     X, Y, params); % (SBE)                                
xbin = parent_pop(:, nobj+ncon+nreal+1:nobj+ncon+nreal+nbits);  
f1 = zeros(popsize,1);
f2 = zeros(popsize,1); 
parfor i=1:popsize
    [f1(i),f2(i)] = evalind_sbe_capsule(X,Y,xbin(i,:),params);
end
parent_pop(:,1) = f1; 
parent_pop(:,2) = f2; 
parent_pop(:,nobj+ncon+nreal+nbits+1) = 0.0;
clear f1 f2 
poolobj = gcp('nocreate');
%********************************* 

%parent_pop = assign_rank_and_crowding_distance(parent_pop, nobj, ncon,nreal, nbin); % (original) 
parent_pop = assign_rank_and_crowding_distance(parent_pop, nobj, ncon, nreal, nbits); % (SBE)
tEnd = toc(tStart);

% Display message
fprintf('NSGA2 | gen: %d - dataset: %s - fold: %d - exp: %s - time %f\n',...
        1,dataset,numexp,nsga2.tag,tEnd); 
    
for i=2:ngen
    tStart = tic; 
    [child_pop, child_strings, ~, nbincross] = selection(parent_pop, parent_strings, 0.0, pcross_bin, ...
        0.0, nbincross, 0, [], [], nbits, nobj, ncon, nreal, nbin);
    [child_pop, child_strings, nbinmut, ~] = mutation(child_pop, child_strings, popsize, nreal, 0.0, ...
        [], [], 0.0, nobj, ncon, nbin, nbits, pmut_bin, nbinmut, 0);
    if (nbin ~= 0)
        child_pop = decode_sbe(child_pop, child_strings, nobj, ncon, nreal);    
    end
    
    % Evaluation of the population 
    % Instance selection with simple binary encoding 
    %*********************************
    %child_pop = evalind_sbe(child_pop, popsize, nobj, ncon, nreal, nbits,...
    %                                 X, Y, params);
    xbin = child_pop(:, nobj+ncon+nreal+1:nobj+ncon+nreal+nbits);  
    f1 = zeros(popsize,1);
    f2 = zeros(popsize,1); 
    parfor j=1:popsize
        [f1(j),f2(j)] = evalind_sbe_capsule(X,Y,xbin(j,:),params);
    end
    child_pop(:,1) = f1; 
    child_pop(:,2) = f2; 
    child_pop(:,nobj+ncon+nreal+nbits+1) = 0.0;
    clear f1 f2 
    %*********************************                              
                                 
    mixed_pop = [parent_pop; child_pop];
    mixed_pop = assign_rank_and_crowding_distance(mixed_pop, nobj, ncon, nreal, nbits);
    %[parent_pop, parent_strings] = fill_nondominated_sort(mixed_pop, [parent_strings; child_strings], nobj, ncon, nreal, nbin);
    [parent_pop, parent_strings] = fill_nondominated_sort(mixed_pop, [parent_strings; child_strings], nobj, ncon, nreal, nbits);
    clear child_pop mixed_pop
    tEnd = toc(tStart);  
    fprintf('NSGA2 | gen: %d - dataset: %s - fold: %d - exp: %s - time %f\n',...
            i,dataset,numexp,nsga2.tag,tEnd);     
end
    
out.minmaxnorm_stats = [mn;mx];
out.nbinmut = nbinmut;
out.nbincross = nbincross;
out.parent_pop = parent_pop; 
out.dataset = dataset; 
out.numexp = numexp; 
out.params = params; 
out.nsga2 = nsga2; 

delete(poolobj);