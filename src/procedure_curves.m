function [card_plot, avg_opt_gap] = procedure_curves(subset, max_budget, k, means, M)

card_plot = [(0:max_budget) * k ;[k, k, mean(sum(subset,2),3)']];


%%
avg_opt_gap = [(1:max_budget)*k ; zeros(1, max_budget)];
avg_opt_gap(:,1) = [0 ; means(end) - mean(means)];
temp = zeros(max_budget-1, M);
for m = 1:M
    temp(:,m) = subset(:,:,m) * means' ./ sum(subset(:,:,m),2);
end
avg_opt_gap(2,2:max_budget) = means(end) - mean(temp,2)';


end

