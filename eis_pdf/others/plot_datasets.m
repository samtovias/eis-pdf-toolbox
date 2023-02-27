% PLOT_DATASETS Plot the original dataset and the selected subset. 
%   PLOT_DATASETS(X,Y,XS,YS,METHOD) plots the instances of the original dataset and 
%   the selected subset. X is the original dataset of N-by-D size (N instances 
%   and D dimensions) and Y is the class labels vector of N-by-1 size. XS is 
%   the selected subset of Ns-by-D size (Ns selected instances) and YS is the
%   class labels vector of Ns-by-1 size. METHOD is a string with the name of
%   the instance selection algorithm used (optional). 
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
%   See also PLOT_PDFS PLOT_NONDOMSET
    
% ------------------------------------------------------------------------
%   Cinvestav-IPN (Mexico)
%   PLOT_DATASETS Version 1.0 (Matlab R2022a)
%   June 2022
%   Copyright (c) 2022, Samuel Omar Tovias Alanis 
% ------------------------------------------------------------------------ 
    
function plot_datasets(X,Y,XS,YS,method)
if nargin < 5 
    method = 'Instance selection'; 
else
    method = upper(method); 
end
% Classes 
c = max(Y); 
% Define colors 
colors = cell(c,1); 
colors{1} = [0 0.4470 0.7410];      % Blue
colors{2} = [0.8500 0.3250 0.0980]; % Orange 
colors{3} = [0.9290 0.6940 0.1250]; % Yellow 
if c > 3
    for i = 4:c 
        colors{i} = [rand rand rand]; 
    end
end 
% Original dataset 
f = figure; 
f.Position = [167 389 833 409]; 
subplot(1,2,1);
for i = 1:c 
    scatter(X(Y==i,1),X(Y==i,2),'MarkerFaceColor',colors{i},'MarkerEdgeColor',colors{i});
    hold on; 
end  
title('Original dataset (X,Y)'); 
legend('Class 1','Class 2','Class 3','location','southwest');
axis square;
box on; 
% Selected subset  
subplot(1,2,2);
for i = 1:c 
    scatter(XS(YS==i,1),XS(YS==i,2),'MarkerFaceColor',colors{i},'MarkerEdgeColor',colors{i});
    hold on; 
end 
title('Selected subset (XS,YS)');
legend('Class 1','Class 2','Class 3','location','southwest');
axis square; 
box on; 
gtitle = strcat(method,{' '},'example');
sgtitle(gtitle);
end 