function obj=Cantilever_tube_function(x)
% x is a 9 dimensions input vairables 
t=x(:,1);
d=x(:,2);
F1=x(:,3);
F2=x(:,4);
T=x(:,5);
Sy=x(:,6);
P=x(:,7);
L1=x(:,8);
L2=x(:,9);

% the constants
theta1=5*pi/180;
theta2=10*pi/180;

% calculate the Intermediate variables 
I=pi/64*(d.^4-(d-2*t).^4);
J=2*I;
taozx=(T.*d)./(2.*J);
M=F1.*L1.*cos(theta1)+F2.*L2.*cos(theta2);
A=pi/4.*(d.^2-(d-2*t).^2);
h=d/2;
sigmax=(P+F1.*sin(theta1)+F2.*sin(theta2))./A+(M.*h)./I;
% calculate the objective function values 
obj=Sy-sqrt(sigmax.^2+3*taozx.^2);
end 