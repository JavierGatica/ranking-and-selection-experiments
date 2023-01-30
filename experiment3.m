% distribution parameters
a = 5;
b = 2;

% experiment parameters
M = 100;
k = 1000;
variance = 0.1;
n0_range = [20, 40, 100];
n_budgets = size(n0_range,2);
alpha = 0.05;


% KNN parameters
h = 0.005;
disc = 0:h:1; % discretization
n_disc = size(disc,2);
knn_curveSTTB = zeros(1,n_disc);

% STTB parameters


% BiPASS parameters
c = 5;
n_init = 10;
dn = 2;



knn_curvesSTTB = zeros(n_budgets,n_disc);
knn_curvesBiPASS = zeros(n_budgets,n_disc);

j = 1;
resultsSTTB = zeros(2, k*M);
resultsBiPASS = zeros(2, k*M);
for n0 = n0_range
    
    % STTB variables dependent on n0
    t_value = tinv((1-alpha)^(1/(k-1)),n0-1);
    % BiPASS variables dependent on n0
    Nmax = n0*k;
    for m = 1:M
        means = betarnd(a,b,[1,k]); % the values must be between 0 and 1
        
        % for STTB
        outputSTTB = sim_output(means,variance,n0,1);
        mu = mean(outputSTTB);
        sigma = mean(var(outputSTTB));
        resultsSTTB(:,(k*(m-1)+1):(k*m)) = [means ; STTB(mu, sigma, t_value, n0)];
        
        
        % for BiPASS
        simulation = @(dn, active) sim_output(means(active), variance, dn, 1);
        resultsBiPASS(:,(k*(m-1)+1):k*m) = [means ; BiPASS(simulation,c, n_init, dn, Nmax, k) == 0];
    end
    
    % THIS IS THE KNN REGRESSION
    
    knn_curvesSTTB(j,:) = knn_reg(disc,resultsSTTB);
    knn_curvesBiPASS(j,:) = knn_reg(disc, resultsBiPASS);
    j = j+1;
end
%% Plotting

close all
figure
t = tiledlayout(2,1);

ylabel(t, 'Probability of Selection ($P( i \in S | \mu_i = \mu )$)', 'interpreter', 'latex')
xlabel(t, 'Sistem performance ($\mu$)','interpreter','latex')

p1 = nexttile;

lgd = legend('Location', 'NorthOutside', 'Orientation', 'Horizontal');
lgd.Layout.Tile = 'North';

title('STTB')
hold on
points = ~isnan(knn_curvesSTTB);
for i = 1:n_budgets 
    plot(disc(points(i,:)), knn_curvesSTTB(i,points(i,:)))
end
lg = legend(p1, '$n_0=20$', '$n_0=40$', '$n_0=100$');
lg.Interpreter = 'latex';
hold off

nexttile
title('BiPASS')
hold on
points = ~isnan(knn_curvesBiPASS);
for i = 1:n_budgets 
    plot(disc(points(i,:)), knn_curvesBiPASS(i,points(i,:)))
end
hold off


% given that there is no empirical dependence of k for the screening ratio,
% it can be considered a fixed variable of the problem. (k = 1000)