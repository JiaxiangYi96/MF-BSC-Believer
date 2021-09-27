function obj=Infill_MFRU(x,Kriging_model_lf,kriging_model_discrepancy,sample_xh,sample_xl)
[y, mse] = estimated_value(x,Kriging_model_lf,kriging_model_discrepancy);
s=sqrt(max(0,mse));
y=abs(y);
b=[1 2];
  X_temp=y./s;
    p_fail=Gaussian_PDF(X_temp);
   ru=y./(p_fail.*s);
% distance function 
sample_x=[sample_xh;sample_xl];
d0=zeros(size(sample_x,1),size(sample_x,1));

for ii=1:size(sample_x,1)
x_temp=sample_x(ii,:);
x_temp=repmat(x_temp,size(sample_x,1),1);
   d0(ii,:)=sqrt(sum((x_temp-sample_x).^2,2))';
end
d0(d0==0)=Inf;
   
d_thete=min(min(d0));
d1=zeros(size(x,1),1);
for i = 1:size(x,1)
  
     d2=sqrt(sum((repmat(x(i,:),size(sample_x,1),1)-sample_x).^2,2));
     d2= d_thete-min(d2);
       d1(i,:)=max(d2,0);
    
end

obj=ru+d1*10^5;


obj=ru;


end