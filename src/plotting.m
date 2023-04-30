function plotting(subsets, means)
close all
[max_budget, k, M] = size(subsets);
max_budget = max_budget + 1;

% for the |S| plot with respect to the budget
sizes = sum(subsets, 2);
card_means = mean(sizes, 3);
q05 = quantile(sizes, 0.05, 3);
q95 = quantile(sizes, 0.95, 3);
maxs = max(sizes,[],3);
mins = min(sizes,[],3);


figure
plot(2:max_budget, [card_means , q05, q95, maxs, mins])
title('Cardinality of Subset')
xlabel('Budget per system')
ylabel('Cardinality |S|')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

blue = [0.3010 0.7450 0.9330];
orange = [0.9290 0.6940 0.1250];
n_lines = max_budget - 1;

colors = zeros(n_lines, 3);
for i = 1:n_lines
    colors(i,:) = (i / n_lines) * orange + (1 - i / n_lines) * blue; 
end

% for the empirical probability plot evolving with the budget
probabilities = mean(subsets, 3);
figure
hold on
for i = 1:n_lines
    plot(means, probabilities(i,:), 'b--o', 'Color', colors(i,:), 'MarkerSize', 3)
end
colormap(colors);
caxis([2 max_budget]);
cbar = colorbar;
ylabel(cbar, 'Budget per system ($n_0$)', 'interpreter', 'latex', 'FontSize', 14)
title('\textbf{Empirical Probability of Selection}', 'interpreter', 'latex', 'FontSize', 12)
xlabel('System Index ($i$)', 'interpreter', 'latex', 'FontSize', 14)
ylabel('Pr$\{ i \in S \}$', 'interpreter', 'latex', 'FontSize', 14)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% for the optimality gap plot
opt_gap = zeros(max_budget-1, M); % optimality gaps
for n = 1:(max_budget-1)
    for m = 1:M
        opt_gap(n,m) = means(find(subsets(n,:,m),1,'first'));
    end
end
avg_opt_gap = means(end) - mean(opt_gap,2); % average optimality gap for any given budget
q95_opt_gap = means(end) - quantile(opt_gap, 0.95, 2);
q05_opt_gap = means(end) - quantile(opt_gap, 0.05, 2);

figure
plot(2:max_budget, [avg_opt_gap , q95_opt_gap , q05_opt_gap])
title('Optimality gap for subset selection')
xlabel('Budget per system ($n_0$)', 'interpreter', 'latex')
ylabel('$\max_{i \in S} \{\mu_k - \mu_i\}$', 'interpreter', 'latex')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% for the subset size histogram plot
figure
hold on
for j = 5:5:(max_budget - 1)
    [GC, GR] = groupcounts(reshape(sizes(j,1,:),[1, M])'); % average cari
    groupcount = zeros(1,k);
    for i = 1:size(GR,1)
        groupcount(GR(i)) = GC(i);
    end
    area(groupcount, 'FaceColor', colors(j,:), 'LineWidth', 0.01)
    
end
colormap(colors);
caxis([5 max_budget]);
cbar = colorbar;
ylabel(cbar, 'Budget per system ($n_0$)', 'interpreter', 'latex', 'FontSize', 14)
title('\textbf{Cardinality of Selections Histogram}', 'interpreter', 'latex', 'FontSize', 12)
xlabel('Number of systems selected ($|S|$)', 'interpreter', 'latex', 'FontSize', 14)
ylabel('Proportion', 'interpreter', 'latex', 'FontSize', 14)

end



