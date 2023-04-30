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

% BiPASS variables
% these values of c and n0 gives an EFER of 0.05 (Nelson masterclass)
c = 5;
n0 = 10;
dn = 1;

% Screening for STTB
subsetSTTBm = zeros(max_budget-1,k,M);
for m = 1:M 
    for j = 2:max_budget
        mu = mean(output(1:j,:,m));
        sigma = mean(var(output(1:j,:,m)));
        subsetSTTBm(j-1,:,m) = STTB(mu,sigma,t_valuesSTTBm(j),j);
    end
    m
end

% Screening for BiPASS
subsetBiPASS = zeros(M,k);
for m = 1:M
    simulation = @(dn, active) sim_output(means(active), variance, dn, 1);
    subsetBiPASS(m,:) = BiPASS(simulation,c,n0,dn,Nmax,k);
    m
end

% Processing STTB

x_sttbm = (1:max_budget) * k;

card_sttbm = zeros(M,max_budget);
card_sttbm(:,1) = k;
card_sttbm(:,2:end) = reshape(sum(subsetSTTBm,2),[max_budget-1 M])';


avg_opt_gap_sttbm = zeros(M,max_budget);

avg_opt_gap_sttbm(:,1) = means(end) - mean(means);
for m = 1:M
    avg_opt_gap_sttbm(m,2:max_budget) = means(end) - (subsetSTTBm(:,:,m) * means' ./ sum(subsetSTTBm(:,:,m),2))';
end


% Processing BiPASS

x_bipass = unique(subsetBiPASS)';
n_transitions = size(x_bipass,2);


card_bipass = zeros(M,n_transitions);
subsets = subsetBiPASS;
subsets(subsets == 0) = Inf;
for n = 1:n_transitions
    card_bipass(:,n) = sum(subsets >= x_bipass(n),2);
end


avg_opt_gap_bipass = zeros(M,n_transitions);
    
for n = 1:n_transitions
    active_systems = subsets > x_bipass(n);
    avg_opt_gap_bipass(:,n) = means(end) - active_systems * means' ./ sum(active_systems,2);
end

add_rm_paths('remove')

cd ..
save("data/exp1_data.mat")



