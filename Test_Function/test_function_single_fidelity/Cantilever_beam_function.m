function obj=Cantilever_beam_function(x)
% x1 and x2 is the input variables 
x1=x(:,1);
x2=x(:,2);
% g(x)=l/325-w*b*l4/(E*I) I=bh^3/12
% where I=bh3/12, where w, b, l, E, I and h are load per unit, width, span, modulus of elasticity, moment
% of inertia of the cross section and depth, respectively. Assuming E and l are 2.6¡Á104 MPa and 6m respectively, 
% the LSF is reduced to
obj=18.46154-74769.23.*x1./(x2.^3);

end