add_rm_paths('add')

clear

% distribution parameters
a = 4;
b = 2;

% experiment parameters
M = 400;
k = 1000;
variance = 0.5;
n0_range = [200, 400, 1000];
n_budgets = size(n0_range,2);
alpha = 0.05;
means = betainv(linspace(0,1,k),a,b);
n_disc = length(means);

% BiPASS parameters
c = 5;
n_init = 10;
dn = 10;


curvesSTTBm = zeros(n_budgets,n_disc-1);
curvesBiPASS = zeros(n_budgets,n_disc-1);

j = 1;
resultsSTTBm = zeros(2, k*M);
resultsBiPASS = zeros(2, k*M);
for n0 = n0_range
    
    % STTB variables dependent on n0
    t_value = tinv((1-alpha)^(1/(k-1)),n0-1);
    % BiPASS variables dependent on n0
    Nmax = n0*k;
    for m = 1:M
        %means = betarnd(a,b,[1,k]); % the values must be between 0 and 1
        % for STTB
        output1 = sim_output(means,variance,n0,1);
        mu = mean(output1);
        sigma = mean(var(output1));
        resultsSTTBm(:,(k*(m-1)+1):(k*m)) = [means ; STTB(mu, sigma, tinv((1-alpha)^(1/(k-1)),n0-1), n0)];
        
        
        % for BiPASS
        simulation = @(dn, active) sim_output(means(active), variance, dn, 1);
        resultsBiPASS(:,(k*(m-1)+1):k*m) = [means ; BiPASS(simulation,c, n_init, dn, Nmax, k) == 0];
    end
    
    % THIS IS THE KNN REGRESSION
    
    curvesSTTBm(j,:) = emp_cdf(means,resultsSTTBm);
    curvesBiPASS(j,:) = emp_cdf(means, resultsBiPASS);
    j = j+1;
end

add_rm_paths('remove')

cd ..
save("data/exp3_data.mat")

% given that there is no empirical dependence of k for the screening ratio,
% it can be considered a fixed variable of the problem. (k = 1000)