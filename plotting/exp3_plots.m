%% Plotting

clear

cd ..
load("data/exp3_data.mat")

close all

figure;

hold on

plot(0,0,':','col',[0,0,0])
plot(0,0,'--','col',[0,0,0])
plot(0,0,'-','col',[0,0,0])

line_style = [":","--","-"];

points = ~isnan(curvesSTTBm);
for i = 1:n_budgets 
    plot(means(points(i,:)), smooth(curvesSTTBm(i,points(i,:))),"Color","#0072BD",'LineStyle',line_style(i))
end

points = ~isnan(curvesBiPASS);
for i = 1:n_budgets 
    plot(means(points(i,:)), smooth(curvesBiPASS(i,points(i,:))),"Color","#D95319",'LineStyle',line_style(i))
end

xlim([0.5 1]);
ylabel('Probability of Selection ($P( i \in S | \mu_i = \mu )$)', 'interpreter', 'latex')
xlabel('System Performance ($\mu$)','interpreter','latex')

lgd = legend('$N=2\cdot 10^5$', '$N=4 \cdot 10^5$', '$N=10^6$','Location', 'NorthOutside', 'Orientation', 'Horizontal','interpreter','latex');
legend boxoff

plt = gca;

set(gcf,'position',[0,0,350,300])
set(plt, 'FontSize', 14)


hold off

savefig("figures/exp3")
