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

%% General Parameters
k = 1000;
means = arrayfun(@sqrt, 1:k);
variance = 10;
Nmax = 100000;
M = 4;

% marginal and uniform STTB variables
alpha = 0.05;
max_budget = ceil(Nmax/k);
output = sim_output(means, variance, max_budget, M);

t_valuesSTTBm = arrayfun(@(n) tinv((1-alpha)^(1/(k-1)),n-1), 1:max_budget);
t_valuesSTTBu = arrayfun(@(n) tinv((1-alpha)^(1/((k-1)*k)),n-1), 1:max_budget);

% BiPASS variables
% these values of c and n0 gives an EFER of 0.05 (Nelson masterclass)
c = 5;
n0 = 10;
dn = 1;

% marginal and uniform ODP variables
t_valuesODPm = arrayfun(@(n) tinv((1-alpha)^(1/k),n-1), 1:max_budget);
t_valuesODPu = arrayfun(@(n) tinv(0.5*(1+(1-alpha)^(1/k)),n-1), 1:max_budget);

%% Screening for STTB
subsetSTTBm = zeros(max_budget-1,k,M);
subsetSTTBu = zeros(max_budget-1,k,M);
subsetODPm = zeros(max_budget-1,k,M);
subsetODPu = zeros(max_budget-1,k,M);
for m = 1:M 
    for j = 2:max_budget
        mu = mean(output(1:j,:,m));
        sigma = mean(var(output(1:j,:,m)));
        subsetSTTBm(j-1,:,m) = STTB(mu,sigma,t_valuesSTTBm(j),j);
        subsetSTTBu(j-1,:,m) = STTB(mu,sigma,t_valuesSTTBu(j),j);
        subsetODPm(j-1,:,m) = ODP(mu,sigma,t_valuesODPm(j),j);
        subsetODPu(j-1,:,m) = ODP(mu,sigma,t_valuesODPu(j),j);
    end
end

%% Screening for BiPASS
subsetBiPASS = zeros(M,k);
for m = 1:M
    simulation = @(dn, active) sim_output(means(active), variance, dn, 1);
    subsetBiPASS(m,:) = BiPASS(simulation,c,n0,dn,Nmax,k);
end

%% Plotting

figure
t = tiledlayout(1,2);


[card_sttbm, avg_opt_gap_sttbm] = procedure_curves(subsetSTTBm, max_budget, k, means, M);
[card_sttbu, avg_opt_gap_sttbu] = procedure_curves(subsetSTTBu, max_budget, k, means, M);
[card_odpm, avg_opt_gap_odpm] = procedure_curves(subsetODPm, max_budget, k, means, M);
[card_odpu, avg_opt_gap_odpu] = procedure_curves(subsetODPu, max_budget, k, means, M);

%% CARDINALITY PLOT BIPASS

transitions = unique(subsetBiPASS)';
n_transitions = size(transitions,2);
card_bipass = [transitions ; zeros(1,n_transitions)];
for n = 1:n_transitions
    card_bipass(2,n) = mean(sum(subsetBiPASS >= transitions(n),2));
end

% AVERAGE OPTIMALITY GAP BIPASS

avg_opt_gap_bipass = [transitions ; zeros(1,n_transitions)];

for n = 1:n_transitions
    macrorep = subsetBiPASS;
    macrorep(macrorep == 0) = Inf;
    macrorep = macrorep > transitions(n);
    avg_opt_gap_bipass(2,n) = means(end) - mean(macrorep * means' ./ sum(macrorep,2));
end

%%

p1 = nexttile;
hold on

lgd = legend('Location', 'NorthOutside', 'Orientation', 'Horizontal');
lgd.Layout.Tile = 'North';

xlabel(t, 'Budget Used ($n_0$)', 'interpreter', 'latex')

stairs(card_sttbm(1,:), card_sttbm(2,:))
stairs(card_sttbu(1,:), card_sttbu(2,:))
stairs(card_odpm(1,:), card_odpm(2,:))
stairs(card_odpu(1,:), card_odpu(2,:))
stairs(card_bipass(1,:), card_bipass(2,:))
plot(n0*k, k, 'r*');

legend(p1, 'marginal STTB', 'uniform STTB', 'marginal ODP', 'uniform ODP', 'BiPASS', 'End of BiPASS Initial Stage');

ylabel('$|S|$', 'interpreter', 'latex')
hold off



nexttile
hold on

stairs(avg_opt_gap_sttbm(1,:),avg_opt_gap_sttbm(2,:))
stairs(avg_opt_gap_sttbu(1,:),avg_opt_gap_sttbu(2,:))
stairs(avg_opt_gap_odpm(1,:),avg_opt_gap_odpm(2,:))
stairs(avg_opt_gap_odpu(1,:),avg_opt_gap_odpu(2,:))
stairs(avg_opt_gap_bipass(1,:), avg_opt_gap_bipass(2,:))
plot(n0*k, means(end) - mean(means), 'r*')
% lgd = legend('STTB', 'BiPASS', 'End of BiPASS Sampling');
% lgd.Interpreter = 'latex';
ylabel('$ \frac{1}{|S|} \sum_{i \in S}( \mu_k - \mu_i )$', 'interpreter', 'latex')
hold off