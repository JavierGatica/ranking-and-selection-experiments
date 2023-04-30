k = 500;
means = arrayfun(@(i) sqrt(i), 1:k);
variance = 1;
alpha = 0.05;
B = 100;
delta = 1; % selects top 5 systems for this configuration
n0 = 200;
M = 20;
output = sim_output(means, variance, n0, M);

subBS = false(M,k);
subSTTB = false(M,k);

sample_means = reshape(mean(output),k,M)';
S2 = reshape(mean(var(output),2),1,M);
t_value = tinv((1-alpha)^(1/((k-1)*k)),n0-1) + delta*sqrt(n0/2)./sqrt(S2);
for m = 1:M
    subBS(m,:) = bs_delta(output(:,:,m), alpha, B, delta);
    subSTTB(m,:) = STTB(sample_means(m,:),S2(m), t_value(m), n0);
end

good_systems = means >= means(end) - delta;
est_pcsBS = sum(all(subBS(:,good_systems)')) / M;
est_pcsSTTB = sum(all(subSTTB(:,good_systems)')) / M;

close all

figure
bar([est_pcsBS, est_pcsSTTB])
set(gca,'xticklabel',{'Bootstrap','STTB'})

figure
hold on
plot(mean(subSTTB))
plot(mean(subBS))
xlabel('system (k)')
ylabel('probability of selection')
hold off

