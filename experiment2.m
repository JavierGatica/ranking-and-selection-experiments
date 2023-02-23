close all
% beta parameters
a = 4;
b = 2;

% experiment parameters
kmin = 100;
step = 1;
kmax = 1000;
k_values = kmin:step:kmax;
n_iterations = size(k_values,2);
variance = 1;
n0 = 40;
M = 1;

% STTB parameters
alpha = 0.05;
t_valuesSTTBm = arrayfun(@(k) tinv((1-alpha)^(1/(k-1)),n0-1), kmin:step:kmax);
t_valuesSTTBu = arrayfun(@(k) tinv((1-alpha)^(1/((k-1)*k)),n0-1), kmin:step:kmax);
t_valuesODPm = arrayfun(@(k) tinv((1-alpha)^(1/k),n0-1), kmin:step:kmax);
t_valuesODPu = arrayfun(@(k) tinv(0.5*(1+(1-alpha)^(1/k)),n0-1), kmin:step:kmax);

% BiPASS parameters
c = 5;
n_init = 10;
dn = 10;

fraction_sttbm = zeros(M, n_iterations);
fraction_sttbu = zeros(M, n_iterations);
fraction_odpm = zeros(M, n_iterations);
fraction_odpu = zeros(M, n_iterations);
fraction_bipass = zeros(M, n_iterations);
i = 1;
for k = kmin:step:kmax
    N = n0 * k;
    means = betacdf(linspace(0,1,k),a,b);
    % for STTB
    output1 = sim_output(means,variance,n0,M);
    for m = 1:M
        mu = mean(output1(:,:,m));
        sigma = mean(var(output1(:,:,m)));
        fraction_sttbm(m,i) = (sum(STTB(mu, sigma, t_valuesSTTBm(i),n0))-1) / (k-1);
        fraction_sttbu(m,i) = (sum(STTB(mu, sigma, t_valuesSTTBu(i),n0))-1) / (k-1);
        fraction_odpm(m,i) = (sum(ODP(mu, sigma, t_valuesODPm(i),n0))-1) / (k-1);
        fraction_odpu(m,i) = (sum(ODP(mu, sigma, t_valuesODPu(i),n0))-1) / (k-1);
    end
    
    % for BiPASS
    simulation = @(dn, active) sim_output(means(active), variance, dn, 1);
    for m = 1:M
        bipass = BiPASS(simulation, c, n_init, dn, N, k);
        bipass = bipass == 0;
        fraction_bipass(m,i) = (sum(bipass)-1) / (k-1);
    end
    i = i+1;
end
fraction_sttbm = mean(fraction_sttbm,1);
fraction_sttbu = mean(fraction_sttbu,1);
fraction_odpm = mean(fraction_odpm,1);
fraction_odpu = mean(fraction_odpu,1);
fraction_bipass = mean(fraction_bipass,1);

%% Plotting

figure
hold on
plot(kmin:step:kmax, fraction_sttbm)
plot(kmin:step:kmax, fraction_sttbu)
plot(kmin:step:kmax, fraction_odpm)
plot(kmin:step:kmax, fraction_odpu)
plot(kmin:step:kmax, fraction_bipass)
xlabel('Number of systems ($k$)', 'interpreter', 'latex')
ylabel('Screening proportion $\bigg( \frac{|S|-1}{k-1}\bigg)$', 'interpreter', 'latex')
ylim([0,1])
legend('marginal STTB', 'uniform STTB', 'marginal ODP', 'uniform ODP','BiPASS')