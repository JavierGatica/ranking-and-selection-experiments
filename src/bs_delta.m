function [subset, G] = bs_delta(output, alpha, B, delta)

[n0, k] = size(output);

bs_sample = zeros(B,k,n0);
for i = 1:k
    bs_sample(:,i,:) = reshape(randsample(output(:,i),n0*B,true),B,1,n0);
end

G = mean(bs_sample,3); % element b,i gives the b-th bootstrap estimation for system i
% now we define the G matrix from my notes
theta_max = max(G,[],2);
for b = 1:B
    G(b,:) = arrayfun(@(theta_bi) theta_bi >= theta_max(b) - delta, G(b,:));
end
G = logical(G);

subset = false(1,k);
i_star = sortrows([(1:k) ; sum(G)]',2,'descend');
for i = 1:k
    subset(i_star(i,1)) = true;
    pcs = sum(arrayfun(@(b) all(subset(G(b,:))), 1:B)) / B;
    if pcs > 1-alpha
        break
    end
end
end