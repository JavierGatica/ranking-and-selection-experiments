add_rm_paths('add')

k = 100;
means = 1:k;


M = 100;

n_min = 5;
n_max = 15;
n_diff = n_max-n_min+1;

alpha = 0.05;
alpha_screening = 1 - sqrt(1 - alpha);
alpha_selection = 1 - sqrt(1 - alpha);

delta = 5.5;
n0_selection = 20;

var_min = 30;
var_max = 35;
n_sigma = 20;
variances = linspace(var_min,var_max,n_sigma);

subset_size = ones(M,n_max-n_min+1,n_sigma);
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
                subset_size(m,i,n_var) = k2;
                [h,~] = Rinott_Number(sum(subset),n0_selection,1-alpha_selection,0.99,1000);
                budget_used(m,i,n_var) = selection(subset, h, n0_selection, delta,variance) + (n_min + i - 1) * k;
            else
                budget_used(m,i,n_var) = (n_min + i - 1) * k;
            end
        end
    end
    n_var
end

add_rm_paths('remove')

cd ..
save("data/exp4a_data.mat")




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
    % the screening procedure used is STTB with marginal PGS
    
    k = size(sample_means,2);
    t_value = tinv((1-alpha)^(1/(k-1)),n0-1);
    W = t_value * sqrt(2 * S2 / n0) + delta;
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


    