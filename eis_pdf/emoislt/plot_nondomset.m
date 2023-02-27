% PLOT_NONDOMSET Plot the non-dominated solutions. 
%   PLOT_NONDOMSET(POP,P) plots the non-dominated solutions. POP is the decoded 
%   population matrix of NP-by-(NOBJ+NVAR+2) size, where NP is the population 
%   size, NOBJ and NVAR indicates the number of objective functions and variables 
%   of the optimization problem. For each row of POP, the first NOBJ elements 
%   are the fitness values regarding each objective function, the next NVAR 
%   numbers corresponds to the values of the solution (cut-off levels), and 
%   the last two values are the rank and the crowding distance, respectively.  
%   P is an array that contains the index of the non-dominated solutions in 
%   the current population.
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
%   xsol = pop(:,nobj+1:nobj+nvar);         % Obtain the decodified cut-off levels
%   for i=1:np                              % Evaluate all the individuals of the population
%       [pop(i,1),pop(i,2)] = ...           
%       evalind(X,Y,xsol(i,:),params);      
%   end                                                        
%   pop = rank_crowddist(pop,nobj,nvar);    % Compute rank and crowding distance   
%   fpop = pop(:,1:nobj);                   % Extract the fitness values 
%   p = simple_nondominated_sort(np,fpop);  % Obtain the non-dominated set
%   plot_nondomset(pop,p);                  % Plot the non-dominated set
%   
%   See also RANK_CROWDDIST SIMPLE_NONDOMINATED_SORT
    
% ------------------------------------------------------------------------
%   Cinvestav-IPN (Mexico)
%   PLOT_NONDOMSET Version 1.0 (Matlab R2022a)
%   June 2022
%   Copyright (c) 2022, Samuel Omar Tovias Alanis 
% ------------------------------------------------------------------------

function plot_nondomset(pop,p)
f = figure; 
f.Position = [40 350 1150 450];
if nargin == 2
    rank = zeros(size(pop,1),1); 
    rank(p) = 1;
    legend1 = 'Non-dominated solutions';
    legend2 = 'Dominated solutions';
else
    rank = pop(:,end-1);
    legend1 = 'Solutions with rank = 1';
    legend2 = 'Solutions with rank > 1';
end
% Plot all the individuals of pop 
subplot(1,2,1);
scatter(pop(:,1),pop(:,2),...
     'MarkerEdgeColor','k');
legend('All solutions');
xlabel('$f_1$','Interpreter','latex','FontSize',14);
ylabel('$f_2$','Interpreter','latex','FontSize',14);
box on; 
% Plot the non-dominated set of pop 
subplot(1,2,2);
scatter(pop(rank==1,1),pop(rank==1,2),'MarkerEdgeColor','b');
hold on;
scatter(pop(rank~=1,1),pop(rank~=1,2),'MarkerEdgeColor','r');
hold on;
legend(legend1,legend2);
xlabel('$f_1$','Interpreter','latex','FontSize',14);
ylabel('$f_2$','Interpreter','latex','FontSize',14);
box on; 
end