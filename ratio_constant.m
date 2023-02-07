function [c_sttb, c_bipass] = ratio_constant(n0)

% beta parameters
a = 4;
b = 2;

% experiment parameters
kmin = 100;
step = 10;
kmax = 500;
k_values = kmin:step:kmax;
n_iterations = size(k_values,2);
variance = 1;
M = 10;

% STTB parameters
alpha = 0.05;
t_values = arrayfun(@(k) tinv((1-alpha)^(1/(k-1)),n0-1), kmin:step:kmax);

% BiPASS parameters
c = 5;
n_init = 10;
dn = 10;

fraction_sttb = zeros(M, n_iterations);
fraction_bipass = zeros(M, n_iterations);
i = 1;
for k = kmin:step:kmax
    N = n0 * k;

    means = [sort(betarnd(a,b,[1, k-1])), 1] * 2;
    % for STTB
    outputSTTB = sim_output(means,variance,n0,M);
    for m = 1:M
        mu = mean(outputSTTB(:,:,m));
        sigma = mean(var(outputSTTB(:,:,m)));
        fraction_sttb(m,i) = (sum(STTB(mu, sigma, t_values(i),n0))-1) / (k-1);
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
fraction_sttb = mean(fraction_sttb,1);
fraction_bipass = mean(fraction_bipass,1);

c_sttb = mean(fraction_sttb);
c_bipass = mean(fraction_bipass);
end

