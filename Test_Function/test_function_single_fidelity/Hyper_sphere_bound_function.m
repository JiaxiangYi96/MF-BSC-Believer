function obj=Hyper_sphere_bound_function(x)
x1=x(:,1);
x2=x(:,2);
obj=1-x1.^3-x2.^3;
end