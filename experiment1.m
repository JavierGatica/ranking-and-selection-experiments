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
k = 10000;
means = arrayfun(@sqrt, 1:k);
variance = 10;
Nmax = 1000000;
M = 20;

% STTB variables
alpha = 0.05;
max_budget = ceil(Nmax/k);
outputSTTB = sim_output(means, variance, max_budget, M);

t_values = arrayfun(@(n) tinv((1-alpha)^(1/(k-1)),n-1), 1:max_budget);

% BiPASS variables
% these values of c and n0 gives an EFER of 0.05 (Nelson masterclass)
c = 5;
n0 = 10;
dn = 1;

%% Screening for STTB
subsetSTTB = zeros(max_budget-1,k,M);
for m = 1:M 
    for j = 2:max_budget
        mu = mean(outputSTTB(1:j,:,m));
        sigma = mean(var(outputSTTB(1:j,:,m)));
        subsetSTTB(j-1,:,m) = STTB(mu,sigma,t_values(j),j);
    end
end

%% Screening for BiPASS
subsetBiPASS = zeros(M,k);
for m = 1:M
    simulation = @(dn, active) sim_output(means(active), variance, dn, 1);
    subsetBiPASS(m,:) = BiPASS(simulation,c,n0,dn,Nmax,k);
end

%% Plotting

% Cardinality Plot
figure
t = tiledlayout(1,2);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CARDINALITY PLOT STTB

card_sttb = [(0:max_budget) * k ;[k, k, mean(sum(subsetSTTB,2),3)']];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CARDINALITY PLOT BIPASS
transitions = unique(subsetBiPASS)';
n_transitions = size(transitions,2);
card_bipass = [transitions ; zeros(1,n_transitions)];
for n = 1:n_transitions
    card_bipass(2,n) = mean(sum(subsetBiPASS >= transitions(n),2));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOTTING
p1 = nexttile;
hold on

lgd = legend('Location', 'NorthOutside', 'Orientation', 'Horizontal');
lgd.Layout.Tile = 'North';

xlabel(t, 'Budget Used ($N$)', 'interpreter', 'latex')
%ylabel(t, 'Procedure Performance Measures', 'interpreter', 'latex')

stairs(card_sttb(1,:), card_sttb(2,:))
stairs(card_bipass(1,:), card_bipass(2,:))
plot(n0*k, k, 'r*');
% lgd = legend('STTB', 'BiPASS', 'End of BiPASS Sampling');
% lgd.Interpreter = 'latex';
legend(p1, 'STTB', 'BiPASS', 'End of BiPASS Initial Stage');

ylabel('$|S|$', 'interpreter', 'latex')
hold off

% Optimality Gap Plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OPTIMALITY GAP PLOT STTB
% 
% opt_gap_sttb = [(0:max_budget) * k ; zeros(1,max_budget+1)]; % optimality gaps
% systems_in_set = zeros(1,M);
% opt_gap_sttb(2,1:2) = ones(1,2) * means(1);
% for n = 2:max_budget
%     for m = 1:M
%         systems_in_set(m) = means(find(subsetSTTB(n-1,:,m),1,'first'));
%     end
%     opt_gap_sttb(2,n+1) = mean(systems_in_set);
% end
% opt_gap_sttb(2,:) = means(end) - opt_gap_sttb(2,:); % average optimality gap for any given budget

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OPTIMALITY GAP PLOT BIPASS

% dado un budget transitions(n), quiero retornar los M indices
% con el sistema con media mas pequeña que todavía esté en el conjunto
% 
% opt_gap_bipass = [transitions ; zeros(1,n_transitions)];
% opt_gap_bipass(2,1) = means(1);
% for n = 2:n_transitions
%     for m = 1:M
%         nonzero_value = find(subsetBiPASS(m,:),1,'last');
%         systems_in_set(m) = means(find(subsetBiPASS(m,1:nonzero_value) <= transitions(n), 1, 'last'));
%     end
%     opt_gap_bipass(2,n) = mean(systems_in_set);
% end
% 
% opt_gap_bipass(2,:) = means(end) - opt_gap_bipass(2,:);
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%
% % PLOTTING
% 
% nexttile
% hold on
% stairs(opt_gap_sttb(1,:), opt_gap_sttb(2,:))
% stairs(opt_gap_bipass(1,:), opt_gap_bipass(2,:))
% plot(n0*k, means(end) - means(1), 'r*')
% ylabel('$\max_{i \in S} \{ \mu_k - \mu_i \}$', 'interpreter', 'latex')
% hold off

% Average Optimality Gap Plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AVERAGE OPTIMALITY GAP STTB

avg_opt_gap_sttb = [(1:max_budget)*k ; zeros(1, max_budget)];
avg_opt_gap_sttb(:,1) = [0 ; means(end) - mean(means)];
temp = zeros(max_budget-1, M);
for m = 1:M
    temp(:,m) = subsetSTTB(:,:,m) * means' ./ sum(subsetSTTB(:,:,m),2);
end
avg_opt_gap_sttb(2,2:max_budget) = means(end) - mean(temp,2)';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AVERAGE OPTIMALITY GAP BIPASS

avg_opt_gap_bipass = [transitions ; zeros(1,n_transitions)];

for n = 1:n_transitions
    macrorep = subsetBiPASS;
    macrorep(macrorep == 0) = Inf;
    macrorep = macrorep > transitions(n);
    avg_opt_gap_bipass(2,n) = means(end) - mean(macrorep * means' ./ sum(macrorep,2));
end

%%%%%%%%%%%%%%%%%
% PLOTTING

nexttile
hold on

stairs(avg_opt_gap_sttb(1,:),avg_opt_gap_sttb(2,:))
stairs(avg_opt_gap_bipass(1,:), avg_opt_gap_bipass(2,:))
plot(n0*k, means(end) - mean(means), 'r*')
% lgd = legend('STTB', 'BiPASS', 'End of BiPASS Sampling');
% lgd.Interpreter = 'latex';
ylabel('$ \frac{1}{|S|} \sum_{i \in S}( \mu_k - \mu_i )$', 'interpreter', 'latex')
hold off