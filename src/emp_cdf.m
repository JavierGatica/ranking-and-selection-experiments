function curve = emp_cdf(disc,results)
n_disc = length(disc);
curve = zeros(1,n_disc-1);
for i = 1:(n_disc-1) % change limits if the domains of the elements mu_i change
    local_results = results(:,disc(i) <= results(1,:) & results(1,:) <= disc(i+1));
    if isempty(local_results)
        curve(i) = NaN;
    else
        curve(i) = mean(local_results(2,:));
    end
end
end

    