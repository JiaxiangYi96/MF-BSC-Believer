function obj=Modified_Rastrigin_function(x)

x1=x(:,1);
x2=x(:,2);

obj1=x1.^2-5*cos(2*pi*x1);

obj2=x2.^2-5*cos(2*pi*x2);

obj=10-obj1-obj2;

end