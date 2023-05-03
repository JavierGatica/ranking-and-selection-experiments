add_rm_paths('add')

clear
close all
% beta parameters
a = 4;
b = 2;

% experiment parameters
k_values = linspace(200,2000,10);
n_iterations = size(k_values,2);
variance = 0.5;
n0 = 40;
M = 40;

% STTB parameters
alpha = 0.05;
t_valuesSTTBm = arrayfun(@(k) tinv((1-alpha)^(1/(k-1)),n0-1), k_values);

% BiPASS parameters
c = 5;
n_init = 10;
dn = 10;

fraction_sttbm = zeros(M, n_iterations);
fraction_bipass = zeros(M, n_iterations);

for i = 1:n_iterations
    N = n0 * k_values(i);
    means = betainv(linspace(0,1,k_values(i)),a,b);
    % for STTB
    output1 = sim_output(means,variance,n0,M);
    for m = 1:M
        mu = mean(output1(:,:,m));
        sigma = mean(var(output1(:,:,m)));
        fraction_sttbm(m,i) = (sum(STTB(mu, sigma, t_valuesSTTBm(i),n0))-1) / (k_values(i)-1);
    end
    
    % for BiPASS
    simulation = @(dn, active) sim_output(means(active), variance, dn, 1);
    for m = 1:M
        bipass = BiPASS(simulation, c, n_init, dn, N, k_values(i));
        bipass = bipass == 0;
        fraction_bipass(m,i) = (sum(bipass)-1) / (k_values(i)-1);
    end
    k_values(i)
end

add_rm_paths('remove')

cd ..
save("data/exp2_data.mat")
