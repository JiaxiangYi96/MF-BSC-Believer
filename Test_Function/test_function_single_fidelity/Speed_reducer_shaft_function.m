function obj=Speed_reducer_shaft_function(x)
% x1 to x6 is the input variables 
D=x(:,1); % Diameter D(mm)
L=x(:,2); % Span L(mm)
F=x(:,3); % External Force F(N)
T=x(:,4); % Torques T(Nmm)
S=x(:,5); % Strength S(MPa)
obj=S-32/(pi.*D.^3).*sqrt(F.^2.*L.^2/16+T.^2);

end