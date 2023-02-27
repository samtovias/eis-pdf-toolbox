% PLOT_PDFS Plot the univariate PDFs of a dataset.
%   PLOT_PDFS(PDF,PARAMS) Plots the univariate PDFs of a dataset with 
%   C classes and D features. PDF is a cell array of C-by-D size that 
%   contains all the PDFs per class and dimension. PARAMS is a structure 
%   that contains the required parameters to plot.
%   
%   Example:
%   --------
%   load concentric3.mat;              % Load dataset 
%   [data,params] = setup_moislt(X,Y); % Setup parameters  
%   PDF = params.PDF;                  % PDFs of the dataset 
%   plot_pdfs(PDF,params);             % Plot PDFs  
%   
%   See also SETUP_MOISLT PLOT_DATASETS  
    
% ------------------------------------------------------------------------
%   Cinvestav-IPN (Mexico)
%   PLOT_PDFS Version 1.0 (Matlab R2022a)
%   June 2022
%   Copyright (c) 2022, Samuel Omar Tovias Alanis 
% ------------------------------------------------------------------------
    
function plot_pdfs(PDF,params) 
figure; 
cont = 0; 
rows = min(params.c,8); 
cols = min(params.d,8);
for i = 1:rows
    for j = 1:cols 
        cont = cont + 1; 
        subplot(rows,cols,cont); 
        plot(params.xh,PDF{i,j}); 
        if i == 1 
            title(strcat('Feature',{' '},num2str(j)))
        end
        if j == 1 
            ylabel(strcat('Class',{' '},num2str(i))); 
        end    
    end 
end
sgtitle('Probability density functions'); 
end 