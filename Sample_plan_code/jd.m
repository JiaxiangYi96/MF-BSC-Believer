function [J,distinct_d]=jd(X,p)
% Computes the distances between all pairs of points in a sampling
% plan X using the p�Cnorm, sorts them in ascending order and
% removes multiple occurrences.
%%
%Inputs:
% X �C sampling plan being evaluated
% p �C distance norm (p=1 rectangular �C default, p=2 Euclidean)
%%
%Outputs:
% J �C multiplicity array (that is, the number of pairs
% separated by each distance value).
% distinct_d �C list of distinct distance values
if ~exist('p','var')
p=1;
end
% Number of points in the sampling plan
n=size(X,1);
% Compute the distances between all pairs of points
d=zeros(1,n*(n-1)/2);
for i=1:n-1
for j=i+1:n
% Distance metric: p�Cnorm
d((i-1)*n-(i-1)*i/2+j-i)=norm(X(i,:)-X(j,:),p);
end
end
% Remove multiple occurrences
distinct_d=unique(d);
% Pre-allocate memory for J
J =zeros(size(distinct_d));
% Generate multiplicity array
for i=1:length(distinct_d)
% J(i) will contain the number of pairs separated
% by the distance distinct_d(i)
J(i)=sum(ismember(d,distinct_d(i)));
end