function subset = ODP(sample_means, S2, t_value, n0)
k = size(sample_means,2);
W = 2 * t_value * sqrt(S2 / n0);
subset = false(1,k);
for i = 1:k
    subset(i) = ~any(sample_means(i) < sample_means - W);
end
end