function [Hstar,Hupper]= Rinott_Number(k, n0, pstar, conf, rep)
%  function to return estimate of Rinott h and upper confidence bound on it
%   k = number of systems
%   n0 = first-stage sample size
%   pstar = 1-alpha value (PCS) (make this an array for fixed k)
%   rep = number of replications to use for estimate
%   conf = confidence level on upper bound
Z = randn(rep,k-1);
Y = chi2rnd(n0-1,rep,k-1);
C = chi2rnd(n0-1,rep,1);
Cmat=repmat(C,1,k-1);
denom = sqrt((n0 - 1.) * (1./Y + 1./Cmat));
H=sort(max(Z.*denom,[],2));
Hstar=quantile(H,pstar);
upper=ceil(pstar * rep +  norminv(conf) * sqrt(pstar * (1. - pstar) *rep) + 0.5)-1;
Hupper= H(upper);


% generating for a range of k:

% fix k_max 

% slice Z,Y matrices for values of k < k_max
end

