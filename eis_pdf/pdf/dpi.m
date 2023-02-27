% DPI Direct plug-in rule. 
%   H = DPI(XI) computes the bandwidth with the direct plug-in rule for PDF
%   estimation using KDE. XI is an univariate continuous variable of N-by-1
%   or 1-by-N size. H is the bandwidth.
%   
%   Example:
%   ------- 
%   load concentric3.mat            % Load a dataset   
%   xi = X(:,1);                    % Extract values of one feature
%   xi = minmaxnorm(xi);            % Normalize the univariate variable
%   h = dpi(xi);                    % Compute the bandwidth with the DPI rule 
%   xh = linspace(-1.5,1.5,100);    % Generate a linearly space vector 
%   ph = kde(xi,xh,h);              % Compute the density estimation 
%   plot(xh,ph);                    % Plot the estimation 
%   t1 = strcat('KDE with',{' '}); 
%   t2 = upper('DPI');
%   t3 = ' rule, h =';
%   t = strcat(t1,t2,t3,{' '},num2str(h)); 
%   title(t);
%   xlabel('xi'); 
%   ylabel('PDF'); 
%   legend('PDF of xi');
%   
%   See also KDE SILVERMAN
%   
%   
%   References:
%   ---------
%   Wand, M. P., & Jones, M. C. (1994). Kernel Smoothing. Chapman & Hall/CRC. 
%   pp. 72
%   
%   Sheather, S.J. and Jones, Chris (1991). A reliable data-based bandwidth 
%   selection method for kernel density estimation. Journal of the Royal Statistical 
%   Society: Series B (Statistical Methodology), 53(3) pp. 683–690.
    
% ------------------------------------------------------------------------
%   Cinvestav-IPN (Mexico)
%   DPI Version 1.0 (Matlab R2022a)
%   June 2022
%   Copyright (c) 2022, Samuel Omar Tovias Alanis 
% ------------------------------------------------------------------------
    
function h = dpi(xi)
% Gaussian kernel derivatives 
dK4 = @(x)exp(x.^2.*(-1.0./2.0)).*1.196826841204298-x.^2.*exp(x.^2.*(-1.0./2.0)).*2.393653682408596+x.^4.*exp(x.^2.*(-1.0./2.0)).*3.989422804014327e-1;
dK6 = @(x)exp(x.^2.*(-1.0./2.0)).*(-5.984134206021491)+x.^2.*exp(x.^2.*(-1.0./2.0)).*1.795240261806447e1-x.^4.*exp(x.^2.*(-1.0./2.0)).*5.984134206021491+x.^6.*exp(x.^2.*(-1.0./2.0)).*3.989422804014327e-1; 
% Number of points 
n = length(xi);
% Difference between each pair of points 
xij = bsxfun(@minus,xi,xi'); 
% Step 1: Estimate Psi_8 with median absoulte deviation   
Psi8_NS = 105/(32*pi^0.5*mad(xi,1)^9);
% Step 2: Estimate Psi_6
g1 = (11.9683/(Psi8_NS*n))^(1/9);
Psi_6 = psi_r(dK6,6,g1,xij,n);
% Step 3: Estimate Psi_4 
g2 = (-2.3937/(Psi_6*n))^(1/7);
Psi_4 = psi_r(dK4,4,g2,xij,n);
% Step 4: Compute bandwidth h
h = (0.2821/(Psi_4*n))^(1/5);
end
%***************************
% PSI estimator for DPI rule
function Psi = psi_r(dKr,r,g,xij,n)
xij = xij./g;
S = sum(dKr(xij(:)));
Psi = (n*(n-1))^(-1)*g^(-r-1)*S; 
end   