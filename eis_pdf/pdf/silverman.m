% SILVERMAN Silverman's rule of thumb
%   H = SILVERMAN(XI) computes the bandwidth with the Silverman's rule of thumb 
%   for PDF estimation using KDE. XI is an univariate continuous variable of 
%   N-by-1 or 1-by-N size. H is the bandwidth.  
%   
%   Example:
%   ------- 
%   load concentric3.mat                % Load a dataset   
%   xi = X(:,1);                        % Extract values of one feature
%   xi = minmaxnorm(xi);                % Normalize the univariate variable
%   h = silverman(xi);                  % Compute the bandwidth with the Silverman's rule
%   xh = linspace(-1.5,1.5,100);        % Generate a linearly space vector 
%   ph = kde(xi,xh,h);                  % Compute the density estimation 
%   plot(xh,ph);                        % Plot the estimation 
%   t = 'KDE with Silverman rule, h =';
%   t = strcat(t,{' '},num2str(h)); 
%   title(t);
%   xlabel('xi'); 
%   ylabel('PDF'); 
%   legend('PDF of xi');
%   
%   See also KDE DPI
%   
%   
%   References:
%   ---------
%   Silverman, B. W. (1986). Density estimation for statistics and data 
%   analysis. Springer. 
     
% ------------------------------------------------------------------------
%   Cinvestav-IPN (Mexico)
%   SILVERMAN Version 1.0 (Matlab R2022a)
%   June 2022
%   Copyright (c) 2022, Samuel Omar Tovias Alanis 
% ------------------------------------------------------------------------
    
function h = silverman(xi)
s = std(xi);
n = numel(xi);
h = ((4*s^5)/(3*n))^(1/5);
end 