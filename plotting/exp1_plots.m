cd ..
load("data/exp1_data.mat")

Nmin = 10^4;

% Plot 1

close all

plot(x_sttbm, mean(card_sttbm),"col","#0072BD")


hold on

plot(x_bipass, mean(card_bipass),"col","#D95319")

plot(x_sttbm, quantile(card_sttbm,0.1),"--","col","#0072BD")

plot(x_bipass, quantile(card_bipass,0.1),"--","col","#D95319")


plot(x_sttbm, quantile(card_sttbm,0.9),"--","col","#0072BD")
plot(x_bipass, quantile(card_bipass,0.9),"--","col","#D95319")




xlim([Nmin Nmax])
xlabel('Total Sample Size ($N$)', 'interpreter', 'latex')
legend('STTB', 'Bi-PASS','interpreter','latex');
legend boxoff
ylabel('Subset Size $(|S|)$', 'interpreter', 'latex')


plt = gca;
set(gcf,'position',[0,0,280,200])
set(plt,"XScale","log")
set(plt, 'FontSize', 14)
hold off

savefig("figures/exp1_1")

% Plot 2


plot(x_sttbm, mean(avg_opt_gap_sttbm),"col","#0072BD")

hold on

plot(x_bipass, mean(avg_opt_gap_bipass),"col","#D95319")

plot(x_sttbm, quantile(avg_opt_gap_sttbm,0.1),"--","col","#0072BD")

plot(x_bipass, quantile(avg_opt_gap_bipass,0.1),"--","col","#D95319")


plot(x_sttbm, quantile(avg_opt_gap_sttbm,0.9),"--","col","#0072BD")
plot(x_bipass, quantile(avg_opt_gap_bipass,0.9),"--","col","#D95319")



xlim([Nmin Nmax])
legend('STTB','Bi-PASS','interpreter','latex');
legend boxoff
xlabel('Total Sample Size ($N$)', 'interpreter', 'latex')
ylabel('AOG $(\frac{1}{|S|} \sum_{i \in S} (\mu_k - \mu_i))$', 'interpreter', 'latex')

plt = gca;

set(gcf,'position',[0,0,280,200])
set(plt,"XScale","log")
set(plt, 'FontSize', 14)

hold off
savefig("figures/exp1_2")

close all
