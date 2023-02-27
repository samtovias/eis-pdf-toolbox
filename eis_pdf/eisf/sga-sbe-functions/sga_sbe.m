
% SGA for Instance Selection with Simple Binary Encoding (SBE)

function out = sga_sbe(X,Y,sga,dataset,numexp)

% Normalize in the range [-1,1] 
[X,mn,mx] = minmaxnorm(X);               

% SBE PARAMS   
params.weight = sga.weight;          % SGA weights for objective function
params.N  = size(X,1);               % Cardinality of dataset    
params.D  = size(X,2);	             % Dimensionality of dataset
params.NC = max(Y);                  % Number of classes
params.xh = linspace(-1.5,1.5,100);  % PDF grid              
params.h  = h_estimate(X,Y,params);  % Compute bandwidth for PDF estimate 
params.P  = getPDF(X,Y,params);      % Obtain original PDFs                           

% SGA params
popsize = sga.popsize; 
ngen = sga.ngen; 
sga.nbin = 1; 
nbits = zeros(1,numel(sga.nbin));
for i = 1:numel(nbits)
    sga.nbits(i) = params.N;
    sga.min_binvar(i) = 0;
    sga.max_binvar(i) = 1;
end
sga.bitlength = sum(sga.nbits);
bitlength = sga.bitlength;

% Initialization 
binpop = logical(round(rand(popsize,bitlength))); 

% Evaluate fitness 
fitpop = zeros(popsize,1,'single');
parfor i=1:popsize
    fitpop(i) = evalind_sga_sbe(X,Y,binpop(i,:),params);
end  
poolobj = gcp('nocreate');

% Current best solution  
[fitbest,indbest] = min(fitpop);
binbest = binpop(indbest,:); 

% Convergence curves
curves = zeros(2,ngen+1,'single');
curves(1,1) = fitbest;
curves(2,1) = mean(fitpop);
    
% Generations 
for gen = 1:ngen  
   
   tStart = tic;  
   
   % Selection  
   binpop = selection_sga(binpop,fitpop);
   
   % Crossover 
   [binpop,sga] = crossover_sbe_sga(binpop,sga);
   
   % Mutation  
   [binpop,sga] = mutation_sbe_sga(binpop,sga);
   
   % Evaluate   
   parfor i=1:popsize
        fitpop(i) = evalind_sga_sbe(X,Y,binpop(i,:),params);
   end  
   
   % Best solution after evaluation  
   past_fitbest = fitbest; 
   past_binbest = binbest; 
   [fitbest,indbest] = min(fitpop); 
   binbest = binpop(indbest,:); 
   
   % Elitism  
   index = 1:popsize; 
   index(indbest) = [];
   ind = single(datasample(index,1));
   binpop(ind,:) = past_binbest;
   fitpop(ind) = past_fitbest;
   
   % Get convergence values (best values of current generation)
   fitbest = min(past_fitbest,fitbest);
   curves(1,gen+1) = fitbest;
   curves(2,gen+1) = mean(fitpop);
   tEnd = toc(tStart);
   
    % Current generation      
    fprintf('SGA | gen: %d - dataset: %s - fold: %d - fit: %f - exp: %s - time %f\n',...
            gen,dataset,numexp,fitbest,sga.tag,tEnd);  
   
end

% Experiment final results
[~,i] = min(fitpop);
xbin = binpop(i,:);
XS = X(logical(xbin),:);
YS = Y(logical(xbin));

% Unnormalize selected subset 
XS = minmaxunnorm(XS,[mn;mx]);

% Save results
out.XS = XS; 
out.YS = YS; 
out.xbin = xbin; 
out.binpop = binpop; 
out.fitpop = fitpop; 
out.curves = curves;
out.params = params; 
out.sga = sga;
out.minmaxnorm_stats = [mn;mx];
out.dataset = dataset; 
out.numexp = numexp; 

delete(poolobj);