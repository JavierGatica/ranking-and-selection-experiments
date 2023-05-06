cd ..
load("data/exp4b_data.mat")

close all

[A, G] = meshgrid(alpha_screening, delta_screening);

s = pcolor(A,G,SampleSize);
s.FaceColor = 'interp';
s.LineStyle = 'none';

hold on



contour(A,G,SampleSize,'Color','black')

plot((1 - sqrt(1-alpha))*ones(1,2), [min(delta_screening), max(delta_screening)], '--r', 'LineWidth',2)

[~,I] = min(SampleSize,[], "all");
[r,c] = ind2sub(size(SampleSize), I);
scatter(A(r,c), G(r,c),"gx","LineWidth",5)


xlabel('Screening Significance Level ($\alpha_1$)','interpreter','latex')
ylabel('Screening Tolerance ($\delta_1$)', 'interpreter','latex')
c = colorbar;
c.Label.String = 'Total Sample Size ($N$)';
c.Label.Interpreter = 'latex';

plt = gca;
set(gcf,'position',[0,0,290,230])
set(plt, 'FontSize', 14)

ylim([min(delta_screening), max(delta_screening)])

hold off

savefig("figures/exp4b")