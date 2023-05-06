cd ..
load("data/exp4a_data.mat")

close all

[N, S] = meshgrid((n_min:n_max), variances);

SampleSize = reshape(mean(budget_used,1),n_diff,n_sigma)';

s = pcolor(N,S,SampleSize);
s.FaceColor = 'interp';
s.LineStyle = 'none';


xlabel('Screening Sample Size ($n_0$)', 'interpreter', 'latex')
ylabel('System Variance ($\sigma^2$)', 'interpreter', 'latex')
c = colorbar;
c.Label.String = 'Total Sample Size ($N$)';
c.Label.Interpreter = 'latex';
c.Ruler.Exponent = 3;



plt = gca;
set(gcf,'position',[0,0,280,205])
set(plt, 'FontSize', 14)
hold off

savefig("figures/exp4a_1")



SubsetSize = reshape(mean(subset_size,1),n_diff,n_sigma)';

s = pcolor(N,S,SubsetSize);
s.FaceColor = 'interp';
s.LineStyle = 'none';


xlabel('Screening Sample Size ($n_0$)', 'interpreter', 'latex')
ylabel('System Variance ($\sigma^2$)', 'interpreter', 'latex')
c = colorbar;
c.Label.String = 'Screening Subset Size ($|S|$)';
c.Label.Interpreter = 'latex';


plt = gca;
set(gcf,'position',[0,0,280,205])
set(plt, 'FontSize', 14)
hold off

savefig("figures/exp4a_2")


close all