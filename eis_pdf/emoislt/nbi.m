% NBI Normal boundary intersection
%   Select a solution from the pareto front with the maximal distance of 
%   the orthogonal projection.
function [ind,coord,distance] = nbi(pareto_front)
[~,mx_x] = max(pareto_front(:,1));
[~,mx_y] = max(pareto_front(:,2));
po1 = pareto_front(mx_x,:);
po2 = pareto_front(mx_y,:);
x = [po1(1) po2(1)];
y = [po1(2) po2(2)];
m = (po2(2) - po1(2))/(po2(1) - po1(1));
b = y(1) - m*x(1);
m2 = -1/m;
x2 = zeros(numel(pareto_front(:,1)),2);
y2 = zeros(numel(pareto_front(:,1)),2);
distance = zeros(numel(pareto_front(:,1)),1);
for i = 1:numel(pareto_front(:,1))
    b2 = pareto_front(i,2) - m2*pareto_front(i,1); 
    A = [-m 1;-m2 1]; 
    vb = [b;b2];
    po3 = A\vb; 
    x2(i,:) = [pareto_front(i,1) po3(1)];
    y2(i,:) = [pareto_front(i,2) po3(2)];
    distance(i) = dist(pareto_front(i,:),po3); 
end
[~,ind] = max(distance);
coord.x = x2; 
coord.y = y2;
coord.po1 = po1; 
coord.po2 = po2; 