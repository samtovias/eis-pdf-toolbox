% KDE Kernel density estimation 
%   PH = KDE(XI,XH,H) computes the KDE algorithm to obtain a probability density 
%   function. XI is a random variable of N-by-1 size. XH is the linearly spaced 
%   vector that contains the center of the regions for the KDE method. H is 
%   the bandwidth. PH is a vector of 1-by-100 size that contains the values 
%   of the estimated probability density function.  
%   
%   Example: 
%   --------
%   load concentric3.mat                % Load a dataset
%   xi = X(:,1);                        % Extract the values of the first variable 
%   xi = minmaxnorm(xi);                % Normalize the random variable
%   xh = linspace(-1.5,1.5,100);        % Generate a linearly space vector  
%   h = dpi(xi);                        % Computes the bandwidth with the DPI rule 
%   ph = kde(xi,xh,h);                  % Estimate the PDF with the KDE method  
%   plot(xh,ph);                        % Plot the PDF  
%   title('Kernel density estimation');                       
%   xlabel('xi');                   
%   ylabel('PDF');                  
%   legend('PDF of xi');            
%   
%   See also DPI SILVERMAN 
%   
%   
%   Reference:
%   ---------
%   Duda, R. O., Hart, P. E., & Stork, D. G. (2000). 
%   Pattern Classification (2a ed.). John Wiley & Sons. pp. 164
    
% ------------------------------------------------------------------------
%   Cinvestav-IPN (Mexico)
%   KDE Version 1.0 (Matlab R2022a)
%   June 2022
%   Copyright (c) 2022, Samuel Omar Tovias Alanis 
% ------------------------------------------------------------------------
    
function ph = kde(x,xh,h)
k = kpdf(x,h);
ph = k(xh);
end 
%****************************
% KDE elementwise application  
function k = kpdf(data,h)
% Gaussian kernel
phi = @(x) exp(-.5*x.^2)/sqrt(2*pi);    
% Kernel density estimation
kernel = @(x) mean(phi((x-data)/h)/h);  
% Elementwise application
k = @(x) arrayfun(kernel,x);            
end 