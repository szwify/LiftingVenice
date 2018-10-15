function [res]=pp_terzaghi(xi,tau,varargin)

% SOLUTION OF THE 1D DIFFUSION EQUATION 
% with ct initial pore -pressure - drained at xi =0 - symmetry at xi=1 and
% drained at xi=2 
% note timecale is t/4  because the domain is from 0 to 2.

    if nargin < 3
        n=200;
    else
        n=varargin{1};
    end
    
    res=0.*(xi'*tau);
   
    for k=1:2:n
        res_1= sin(pi * k * xi /2.) ;
        res_2 = exp(-tau*((pi*k)^2.));
        res=res + (4./(pi*k))*(res_1'*res_2);
    end
       
end
