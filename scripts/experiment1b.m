% in this experiment we compare how STTB compares to BiPASS with a fixed
% total budget.

% the parameters of the experiment are

% means: the means of k systems to be tested
% variance: the variance of these systems
% Nmax: total budget, i.e. the total number of simulation outputs to be sampled 
% M: number of macroreplications
%
% for STTB:
%     
% alpha: significance level
% max_budget: samples per system
% output: precomputed simulation output 
%
% for BiPASS
%     
% c: parameter for controlling EFER guarantee
% n0: initial budget per system
% dn: budget per iteration
% simulation: encapsulated function for simulation output

% General Parameters
clear
clc

add_rm_paths('add');

% beta distribution parameters
a = 4;
b = 2;

k = 1000;
means = betainv(linspace(0,1,k),a,b);

variance = 0.5;
Nmax = 1000000;
M = 40;

% marginal and uniform STTB variables
alpha = 0.05;
max_budget = ceil(Nmax/k);
output = sim_output(means, variance, max_budget, M);

t_valuesSTTBm = arrayfun(@(n) tinv((1-alpha)^(1/(k-1)),n-1), 1:max_budget);
t_valuesSTTBu = arrayfun(@(n) tinv((1-alpha)^(1/((k-1)*k)),n-1), 1:max_budget);
t_valuesDSTTBm = arrayfun(@(n) tinv((1-alpha)^(1/k),n-1), 1:max_budget);
t_valuesDSTTBu = arrayfun(@(n) tinv( (1 + (1-alpha)^(1/k))/2 ,n-1), 1:max_budget);


% BiPASS variables
% these values of c and n0 gives an EFER of 0.05 (Nelson masterclass)
c = 5;
n0 = 10;
dn = 1;

% Screening for STTB
subsetSTTBm = zeros(max_budget-1,k,M);
subsetSTTBu = zeros(max_budget-1,k,M);
subsetDSTTBm = zeros(max_budget-1,k,M);
subsetDSTTBu = zeros(max_budget-1,k,M);
for m = 1:M 
    for j = 2:max_budget
        mu = mean(output(1:j,:,m));
        sigma = mean(var(output(1:j,:,m)));
        subsetSTTBm(j-1,:,m) = STTB(mu,sigma,t_valuesSTTBm(j),j);
        subsetSTTBu(j-1,:,m) = STTB(mu,sigma,t_valuesSTTBu(j),j);
        subsetDSTTBm(j-1,:,m) = DSTTB(mu,sigma,t_valuesDSTTBm(j),j);
        subsetDSTTBu(j-1,:,m) = DSTTB(mu,sigma,t_valuesDSTTBu(j),j);
    end
    m
end

% Processing STTB

x = (1:max_budget) * k;

card_sttbm = zeros(M,max_budget);
card_sttbm(:,1) = k;
card_sttbm(:,2:end) = reshape(sum(subsetSTTBm,2),[max_budget-1 M])';

card_sttbu = zeros(M,max_budget);
card_sttbu(:,1) = k;
card_sttbu(:,2:end) = reshape(sum(subsetSTTBu,2),[max_budget-1 M])';

card_dsttbm = zeros(M,max_budget);
card_dsttbm(:,1) = k;
card_dsttbm(:,2:end) = reshape(sum(subsetDSTTBm,2),[max_budget-1 M])';

card_dsttbu = zeros(M,max_budget);
card_dsttbu(:,1) = k;
card_dsttbu(:,2:end) = reshape(sum(subsetDSTTBu,2),[max_budget-1 M])';


add_rm_paths('remove')

cd ..
save("data/exp1b_data.mat")



