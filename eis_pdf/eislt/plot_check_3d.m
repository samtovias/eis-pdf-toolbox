

[~,X2] = pca(X);
[~,XS2] = pca(XS);
method = "MICRO_EISLT";
%plot_datasets(X2,Y,XS2,YS,method);

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
    scatter3(X2(Y==i,1),X2(Y==i,2),X2(Y==i,3),'MarkerFaceColor',colors{i},'MarkerEdgeColor',colors{i});
    hold on; 
end  
title('Original dataset (X,Y)'); 
legend('Class 1','Class 2','Class 3','location','southwest');
axis square;
box on; 
% Selected subset  
subplot(1,2,2);
for i = 1:c 
    scatter3(XS2(YS==i,1),XS2(YS==i,2),XS2(YS==i,3),'MarkerFaceColor',colors{i},'MarkerEdgeColor',colors{i});
    hold on; 
end 
title('Selected subset (XS,YS)');
legend('Class 1','Class 2','Class 3','location','southwest');
axis square; 
box on; 
gtitle = strcat(method,{' '},'example');
sgtitle(gtitle);
