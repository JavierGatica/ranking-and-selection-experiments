function output = sim_output(means, variance, max_budget, M)
k = size(means,2);
output = zeros(max_budget, k, M);
for i = 1:k
    output(:,i,:) = normrnd(means(i), variance, [max_budget 1 M]);
end
end