%% TODO:

% separate calculations from plotting
% 
% make functions for more concise code
% 
% adapt section 2 code to reach the linear growth portion
%
% adapt the code to work for the new folder structure

add_rm_paths('add')

k = 50;
means = 1:k;



M = 100;

n_min = 5;
n_max = 600;
n_diff = n_max-n_min+1;

alpha = 0.05;
alpha_screening = 1 - sqrt(1 - alpha);
alpha_selection = 1 - sqrt(1 - alpha);

delta = 0.05;
n0_selection = 20;

var_min = 0.5;
var_max = 0.7;
n_sigma = 20;
variances = linspace(var_min,var_max,n_sigma);

budget_used = zeros(M,n_max-n_min+1,n_sigma);
for n_var = 1:n_sigma
    variance = variances(n_var);
    parfor m = 1:M
        for i = 1:n_diff
            sample_means = means_simulation(true(1,k),i+n_min-1, means, variance);
            sample_var = var_simulation(true(1,k),i+n_min-1,variance);
            subset = screening(sample_means,sample_var,i+n_min-1,alpha_screening,0);
            k2 = sum(subset);
            if k2 >= 2
                [h,~] = Rinott_Number(sum(subset),n0_selection,1-alpha_selection,0.99,1000);
                budget_used(m,i,n_var) = selection(subset, h, n0_selection, delta,variance) + (n_min + i - 1) * k;
            else
                budget_used(m,i,n_var) = (n_min + i - 1) * k;
            end
        end
    end
end



[N, S] = meshgrid(n_min:n_max, variances);

SampleSize = reshape(mean(budget_used,1),n_diff,n_sigma)';

s = pcolor(N,S,SampleSize);
s.FaceColor = 'interp';
s.LineStyle = 'none';

set(gca,"XScale","log")

xlabel('Screening sample size ($n_0$)', 'interpreter', 'latex')
ylabel('System variance ($\sigma^2$)', 'interpreter', 'latex')
c = colorbar;
c.Label.String = 'Total Sample Size';
c.Label.Interpreter = 'latex';



add_rm_paths('remove')

function sample_means = means_simulation(systems, n0, means, variance)

    if length(systems) ~= length(means)
        error("System indexing does not match length")
    end
    
    sample_means = means(systems) + randn([1,sum(systems)]) * sqrt(variance / n0);
end

function sample_var = var_simulation(systems, n0, variance)

    sample_var = variance * chi2rnd(sum(systems)*n0) / (sum(systems)*n0 - 1);

end


function subset = screening(sample_means,S2,n0,alpha,delta)
    % the screening procedure used is DSTTB with marginal PGS
    
    k = size(sample_means,2);
    t_value = tinv((1-alpha)^(1/k),n0-1);
    W = 2 * t_value * sqrt(S2 / n0) + delta;
    subset = false(1,k);
    for i = 1:k
        subset(i) = ~any(sample_means(i) < sample_means - W);
    end
end


function budget_used = selection(systems, h, n0, delta,variance)
    % the procedure used is Rinott's procedure
    
    sample_var = var_simulation(systems, n0,variance);
    N = ceil(((h / delta) ^ 2) * sample_var);
    budget_used = sum(systems)*N;
end


    