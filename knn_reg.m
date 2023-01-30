function curve = knn_reg(disc,results)
h = disc(2) - disc(1);
n_disc = size(disc,2);
curve = zeros(1,n_disc);
for i = 1:n_disc % change limits if the domains of the elements mu_i change
    local_results = results(:,disc(i) <= results(1,:) & results(1,:) < disc(i)+h);
    if isempty(local_results)
        curve(i) = NaN;
    else
        curve(i) = mean(local_results(2,:));
    end
end
end

