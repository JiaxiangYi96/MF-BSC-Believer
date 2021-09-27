function obj=Circular_Pipe_function(x)
% define the consistent number
R=0.3377;  % radius of the pipe
t=0.03377; % the thickness of the pipe
M=3;     % applied bending moment
% define the random number
sigmaf=x(:,1); % flow stress
theta=x(:,2);  % half-crack angle
  % defining the objective funtion 
  obj=4.*t.*sigmaf.*R.^2.*(cos(theta./2)-0.5.*sin(theta))-M;
end