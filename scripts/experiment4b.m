%% TODO:

% change to have more than one good system

% start alpha_1 from 0

% reduce n_screening for bigger subset size


add_rm_paths('add')

k = 1000;
means = 1:k; % 


n_selection = 20;
rep = 1000;
Z = randn(rep,k-1);
Y = chi2rnd(n_selection-1,rep,k-1);
C = chi2rnd(n_selection-1,rep,1);


variance = 1;

M = 200;
n_screening = 30;
alpha = 0.05;


n_alpha = 100;
alpha_screening = linspace(0.005,0.03,n_alpha);
alpha_selection = 1 - (1 - alpha) * (1 - alpha_screening) .^ -1;

n_delta = 100;
delta = 0.1;
delta_screening = linspace(0,delta/2,n_delta);
delta_selection = delta - delta_screening;

h = zeros(k,n_alpha);
for i = 2:200
    for j = 1:n_alpha
        h(i,j) = CRN_rinott_number(Z,Y,C,i,1 - alpha_selection(j),n_selection);
    end
end


budget_used = ones(M,n_alpha,n_delta) * n_screening * k;
parfor m = 1:M
    sample_means = means_simulation(true(1,k), n_screening,means,variance);
    sample_var = var_simulation(true(1,k), n_screening,variance);
    for j = 1:n_delta
        for i = 1:n_alpha
            subset = screening(sample_means,sample_var,n_screening,alpha_screening(i),delta_screening(j));
            k2 = sum(subset);
            if k2 >= 2
                budget_used(m,i,j) = selection(subset, h(k2,i), n_selection, delta_selection(j),variance);
            end
        end
    end
    m
end

SampleSize = reshape(mean(budget_used,1),n_alpha,n_delta)';
SampleSize = smooth2a(SampleSize, 5,5);


add_rm_paths('remove')

cd ..
save("data/exp4b_data.mat")
%% Function definitions

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


function Hstar = CRN_rinott_number(Z,Y,C,k,pstar,n0)
    % calculate Rinott's number with CRN from a common sample
    Cmat=repmat(C,1,k-1);
    denom = sqrt((Y(:,1:(k-1)) + Cmat) .* (1/(n0-1)));
    H= max(Z(:,1:(k-1))./ denom,[],2);
    Hstar=quantile(H,pstar);
end
    
    