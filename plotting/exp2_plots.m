%% Plotting

close all

cd ..
load("data/exp2_data.mat")

figure;
plot(k_values, mean(fraction_sttbm),"col","#0072BD")
hold on
plot(k_values, mean(fraction_bipass),"col","#D95319")

plot(k_values, quantile(fraction_sttbm,0.1),"--","col","#0072BD")

plot(k_values, quantile(fraction_bipass,0.1),"--","col","#D95319")

plot(k_values, quantile(fraction_sttbm,0.9),"--","col","#0072BD")
xlim([k_values(1) k_values(end)])


plot(k_values, quantile(fraction_bipass,0.9),"--","col","#D95319")
xlabel('Number of Systems ($k$)', 'interpreter', 'latex')
ylabel('Proportion of Additional Systems $\bigg( \frac{|S|-1}{k-1}\bigg)$', 'interpreter', 'latex')
ylim([0,1])
legend('STTB', 'Bi-PASS')
legend boxoff

plt = gca;

set(gcf,'position',[0,0,350,300])
set(plt,"XScale","log")
set(plt, 'FontSize', 14)

savefig("figures/exp2")