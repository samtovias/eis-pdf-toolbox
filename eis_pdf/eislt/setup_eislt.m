% SETUP_EISLT Initial configuration of EISLT.
%   [DATA,PARAMS] = SETUP_EISLT(X,Y,VARARGIN) sets the initial configuration 
%   of the EISLT algorithm. X is a dataset of N-by-D size (N instances and D 
%   features). Y is a class label vector of N-by-1 size. VARARGIN is a cell 
%   array of variable length that contains the input parameters in the following 
%   order: 
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
%   DATA and PARAMS are the normalized dataset in the range [-1,1] and the input 
%   parameters of the optimization problem, respectively.
%   
%   Example:
%   -------
%   load concentric3.mat                            % Load a dataset 
%   [data,params] = setup_eislt(X,Y,{10});          % Setup with: 
%                                                   % - 10 generations 
%                                                   % - All the other parameters are set to default values 
%   
%   See also EISLT MINMAXNORM
    
% ------------------------------------------------------------------------
%   Cinvestav-IPN (Mexico)
%   SETUP_EISLT Version 1.0 (Matlab R2022a)
%   June 2022
%   Copyright (c) 2022, Samuel Omar Tovias Alanis 
% ------------------------------------------------------------------------
    
function [data,params] = setup_eislt(X,Y,varargin)
% Default parameters  
np = 100;               % Size of the population  
gen = 100;              % Number of generations     
w = 0.5;                % Weight of the objective function
pc = 0.9;               % Crossover probability 
pm = '1/L';             % Mutation probability
distance = 'Minkowski'; % Distance metric to build linkage trees 
p = '2';                % Value of the distance metric order  
hopt = 'Silverman';     % Method to compute the bandwidth used in the PDF estimation       
dataset = '';           % Dataset name 
numexp = '';            % Number of the current experiment  
% Passing out the input parameters
if ~isempty(varargin) 
inputs = numel(varargin{1});
params = {'np','gen','w','pc','pm','distance','p','hopt','dataset','numexp'}; 
    for i = 1:inputs
        str = strcat(params{i},' = varargin{1}{1,',num2str(i),'};');
        eval(str);
    end
end
% Normalize the dataset in the range [-1,1] 
[X,mn,mx] = minmaxnorm(X);
data = struct(); data.X = X; data.mn = mn; data.mx = mx;  
% Setup parameters 
params = setparams_eislt(X,Y,np,gen,w,pc,pm,distance,p,hopt,dataset,numexp);
% Check mutation probability value 
if strcmp(params.pm,'1/L'); params.pm = 1/params.chrlen; end 
end 
%***************************
% Setup parameters of EISLT  
function params = setparams_eislt(X,Y,np,gen,w,pc,pm,distance,p,hopt,dataset,numexp)
% General information -------------------                          
params = struct();                      % Initialize the structure of the input parameters
params.algorithm = 'EISLT';             % Name of the proposed algorithm 
params.metaheuristic = 'GA-GOP';        % Global optimization algorithm  
params.numexp = numexp;                 % Number of the current experiment 
% General optimization parameters ------- 
params.np = np;                         % Size of the population 
params.gen = gen;                       % Number of generations 
params.w = w;                           % Weight of the objective function 
params.pc = pc;                         % Crossover probability 
params.pm = pm;                         % Mutation probability 
% Dataset information ------------------- 
params.dataset = dataset;               % Dataset name
params.n = size(X,1);                   % Number of instances 
params.d  = size(X,2);	                % Number of dimensions
params.c = numel(unique(Y));            % Number of classes
params.ni = niclass(Y);                 % Number of instances per class
% Specific optimization parameters -----  
params.nobj = 1;                        % Number of objectives
params.nvar = params.c;                 % Number of optimization variables
params.lmin = ones(1,params.c);         % Lower limit of each variable
params.lmax = params.ni-1;              % Upper limit of each variable 
params.precision = 0;                   % Decimal precision of the variables
params.nbits = getnbits(params);        % Number of bits per variable 
params.chrlen = sum(params.nbits);      % Length of the chromosome 
% PDF estimation parameters ------------- 
params.xh = linspace(-1.5,1.5,100);     % PDF grid                 
params.hopt = hopt;                     % Method to compute the bandwidth used in the PDF estimation
params.H = bandwidths(X,Y,params);      % Compute the bandwidths for the PDF estimation  
params.PDF = pdfs(X,Y,params);          % Obtain the original PDFs using the KDE algorithm 
% Linkage trees construction setup ------   
params.distance = distance;             % Distance metric to build linkage trees 
params.p = p;                           % Value of the distance metric order 
params.cofmin = params.lmin;            % Minimum cutoff values
params.cofmax = params.lmax;            % Maximum cutoff values
params.LT = buildlt(X,Y,params);        % Build linkage trees     
end                                       