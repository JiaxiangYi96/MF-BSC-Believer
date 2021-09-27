clc;
clear all;
result=zeros(11,1);
x=linspace(0,1,11);
test_hf_function='roof_truss';
test_lf_function='Cheap_roof_truss';
[num_vari,z,mu,sigma,design_space]=test_function_multi_fidelity(test_hf_function);
for i=1:1:11
    A=x(i);
num_vari=6;
% number of initial sampling 
num_initial_sample=500;

% sampling and calculate the response 
sample_x = repmat(design_space(1,:),num_initial_sample,1)+repmat(design_space(2,:)-design_space(1,:),num_initial_sample,1).*...
    lhsdesign(num_initial_sample,num_vari,'criterion','maximin','iterations',1000);
y1= feval(test_hf_function,sample_x);
y2=feval(test_lf_function,sample_x,A);

R=r(y1,y2,num_initial_sample);

result(i,:)=R;

end
figure(1)
plot(x,result,'k','LineWidth',1);
set(gca,'FontSize',14,'Fontname', 'Times New Roman');
hold on

