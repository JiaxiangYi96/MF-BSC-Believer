function Phiq=mmphi(X, q, p)
% Calculates the sampling plan quality criterion of Morris and
% Mitchell.
%%
%Inputs:
% X �C sampling plan
% q �C exponent used in the calculation of the metric
% p �C the distance metric to be used (p=1 rectangular �C
% (default, p=2 Euclidean)
%%
%Output:
% Phiq �C sampling plan ��space�Cfillingness�� metric
% Assume defaults if arguments list incomplete
if ~exist('p','var')
p=1;
end
if ~exist('q','var')
q=2;
end
% Calculate the distances between all pairs of
% points (using the p�Cnorm) and build multiplicity array J
%
[J,d]=jd(X,p);
% The sampling plan quality criterion
Phiq=sum(J.*(d.^(-q)))^(1/q);