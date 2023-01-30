function subset = STTB(sample_means, S2, t_value, n0)
k = size(sample_means,2);
W = t_value * sqrt(2 * S2 / n0);
subset = false(1,k);
for i = 1:k
    subset(i) = ~any(sample_means(i) < sample_means - W);
end
end

