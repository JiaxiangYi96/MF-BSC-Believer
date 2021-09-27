function obj = Branin_function(x)
% ----------------------------------------------------------------------------
% constrained Branin Test Function for constrained Nonlinear Optimization
% A. Forrester, A. Sobester, A. Keane, Engineering design via surrogate modelling: 
% a practical guide, John Wiley & Sons, 2008.
% 0 <= x1 <= 1
%  0 <= x2 <= 1
% fmin = 5.5757 where xmin =   [0.96773, 0.20667]
% ----------------------------------------------------------------------------
x1 = 15*x(:,1) - 5;
x2 = 15*x(:,2);
obj = (x2 - (5.1/(4*pi^2))*x1.^2 + (5/pi)*x1 -6).^2 + 10*(1-1/(8*pi))*cos(x1) +10 + 5*x(:,1)-50;


end