cd ..
load("data/exp4a_data.mat")

close all

[N, S] = meshgrid((n_min:n_max)*k, variances);

SampleSize = reshape(mean(budget_used,1),n_diff,n_sigma)';

s = pcolor(N,S,SampleSize);
s.FaceColor = 'interp';
s.LineStyle = 'none';


xlabel('Screening Total Sample Size ($N_1$)', 'interpreter', 'latex')
ylabel('System Variance ($\sigma^2$)', 'interpreter', 'latex')
c = colorbar;
c.Label.String = 'Total Sample Size ($N$)';
c.Label.Interpreter = 'latex';


plt = gca;
set(gcf,'position',[0,0,280,200])
set(plt,"XScale","log")
set(plt, 'FontSize', 14)
hold off

savefig("figures/exp4a_1")



SubsetSize = reshape(mean(subset_size,1),n_diff,n_sigma)';

s = pcolor(N,S,SubsetSize);
s.FaceColor = 'interp';
s.LineStyle = 'none';


xlabel('Screening Total Sample Size ($N_1$)', 'interpreter', 'latex')
ylabel('System Variance ($\sigma^2$)', 'interpreter', 'latex')
c = colorbar;
c.Label.String = 'Screening Subset Size ($|S|$)';
c.Label.Interpreter = 'latex';


plt = gca;
set(gcf,'position',[0,0,280,200])
set(plt,"XScale","log")
set(plt, 'FontSize', 14)
hold off

savefig("figures/exp4a_2")


close all