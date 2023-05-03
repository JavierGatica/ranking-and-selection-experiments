%% TODO

% put more labels in the YAxis

cd ..
clear
load("data/exp1b_data.mat")

Nmin = 10^4;

% Plot 1

close all

hold on

plot(x, mean(card_sttbm),"col","#0072BD")
plot(x, mean(card_sttbu),"col","#EDB120")
plot(x, mean(card_dsttbm),"col","#7E2F8E")
plot(x, mean(card_dsttbu),"col","#77AC30")


% plot(x, quantile(card_sttbm,0.1),"--","col","#0072BD")
% 
% plot(x, quantile(card_sttbm,0.9),"--","col","#0072BD")
% 
% 
% 
% 
% plot(x, quantile(card_sttbu,0.1),"--","col","#EDB120")
% 
% plot(x, quantile(card_sttbu,0.9),"--","col","#EDB120")
% 
% 
% 
% 
% plot(x, quantile(card_dsttbm,0.1),"--","col","#7E2F8E")
% 
% plot(x, quantile(card_dsttbm,0.9),"--","col","#7E2F8E")
% 
% 
% 
% 
% plot(x, quantile(card_dsttbu,0.1),"--","col","#77AC30")
% 
% plot(x, quantile(card_dsttbu,0.9),"--","col","#77AC30")



xlim([Nmin Nmax])
xlabel('Total Sample Size ($N$)', 'interpreter', 'latex')
legend('STTB marginal', 'STTB uniform','DSTTB marginal','DSTTB uniform','interpreter','latex');
legend boxoff
ylabel('Subset Size $(|S|)$', 'interpreter', 'latex')


plt = gca;
set(gcf,'position',[0,0,280,200])
set(plt,"YScale","log")
set(plt, 'FontSize', 14)
hold off


savefig("figures/exp1b")

