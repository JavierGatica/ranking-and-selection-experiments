%% Global variables
k = 1000;
variance = 1;
means = sqrt(1:k);

%% CRN Rinott's number
n_selection = 20;
rep = 12000;
Z = randn(rep,k-1);
Y = chi2rnd(n_selection-1,rep,k-1);
C = chi2rnd(n_selection-1,rep,1);




%% Section 1:

M = 20;

n_min = 10;
n_max = 30;
n_diff = n_max-n_min+1;

alpha = 0.05;
alpha_screening = 1 - sqrt(1 - alpha);
alpha_selection = 1 - sqrt(1 - alpha);

delta = 0.05;
n0_selection = 20;

var_min = 0.5;
var_max = 1.5;
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
                budget_used(m,i,n_var) = selection(subset, h, n0_selection, delta,variance) + k*(n_min+i-1);
            else
                budget_used(m,i,n_var) = k*(n_min+i-1)
            end
        end
        sum(subset)
    end
end


for n_screening = n_min:n_max
    budget_used(:,n_screening-n_min+1) = budget_used(:,n_screening-n_min+1) + n_screening * k;
end

[N, S] = meshgrid(n_min:n_max, variances);

SampleSize = reshape(mean(budget_used,1),n_diff,n_sigma)';

s = pcolor(N,S,SampleSize);
s.FaceColor = 'interp';
s.LineStyle = 'none';

xlabel('Screening sample size ($n_0$)', 'interpreter', 'latex')
ylabel('System variance ($\sigma^2$)', 'interpreter', 'latex')
c = colorbar;
c.Label.String = 'Total Sample Size';
c.Label.Interpreter = 'latex';
%% Section 2:

variance = 1;
M = 50;
alpha = 0.05;
alpha_screening = 1 - sqrt(1 - alpha);
alpha_selection = 1 - sqrt(1 - alpha);
n_screening = 10;
n_selection = 20;
deltamax = 0.1; % max value of delta for screening (the min value is 0) change this for resolution
delta = 0.16;
Drange = linspace(0,deltamax,100);

budget_delta = zeros(M,100);
parfor m = 1:M
    for i = 1:100
        sample_means = means_simulation(true(1,k),n_screening, means, variance);
        S2 = var_simulation(true(1,k),n_screening,variance);
        
        subset = screening(sample_means,S2,n_screening,alpha_screening,deltamax * (i/100));
        k2 = sum(subset);
        if k2 >= 2
            [h,~] = Rinott_Number(k2,n_screening,1-alpha_selection,0.99,1000);
            temp = selection(subset, h, n_selection, delta - deltamax * (i/100),variance);
            budget_delta(m,i) = budget_delta(m,i) + temp;
        end
        budget_delta(m,i) = budget_delta(m,i) + k*n_screening;
    end
end


plot(Drange, mean(budget_delta))
xlabel('Screening tolerance $\delta_1$','interpreter','latex')
ylabel('Total sample size (with $\delta_2 = \delta - \delta_1)$', 'interpreter', 'latex')

%% Section 3: (this experiment combines section 2 and the original section 3)



variance = 1;

M = 50;
n_screening = 30;
alpha = 0.05;


n_alpha = 100;
alpha_screening = linspace(0.008,0.03,n_alpha);
alpha_selection = 1 - (1 - alpha) * (1 - alpha_screening) .^ -1;

n_delta = 100;
delta = 0.1;
delta_screening = linspace(0,delta/2,n_delta);
delta_selection = delta - delta_screening;
%%

h = zeros(k,n_alpha);
for i = 2:200
    for j = 1:n_alpha
        h(i,j) = CRN_rinott_number(Z,Y,C,i,1 - alpha_selection(j),n_selection);
    end
end
%%

budget_used = zeros(M,n_alpha,n_delta);
subset_size = zeros(M,n_alpha,n_delta);
parfor m = 1:M
    sample_means = means_simulation(true(1,k), n_screening,means,variance);
    sample_var = var_simulation(true(1,k), n_screening,variance);
    for j = 1:n_delta
        for i = 1:n_alpha
            subset = screening(sample_means,sample_var,n_screening,alpha_screening(i),delta_screening(j));
            k2 = sum(subset);
            subset_size(m,i,j) = k2;
            if k2 >= 2
                budget_used(m,i,j) = selection(subset, h(k2,i), n_selection, delta_selection(j),variance);
            end
        end
    end
    m
end

%% Plotting
close all

[A, G] = meshgrid(alpha_screening, delta_screening);

SampleSize = reshape(mean(budget_used,1),n_alpha,n_delta)';
SubsetSize = reshape(mean(subset_size,1),n_alpha,n_delta)';

f1 = figure;

% contour plot
% https://www.mathworks.com/help/matlab/ref/contour.html


s = contour(A,G,SampleSize);
% s.FaceColor = 'interp';
% s.LineStyle = 'none';

xlabel('screening significance level ($\alpha_1$)','interpreter','latex')
ylabel('screening selection tolerance ($\delta_2$)', 'interpreter','latex')
c = colorbar;
c.Label.String = 'Average Total Sample Size';
c.Label.Interpreter = 'latex';

%%
f2 = figure;

s = pcolor(A,G,SubsetSize);
s.FaceColor = 'interp';
s.LineStyle = 'none';

xlabel('screening significance level ($\alpha_1$)','interpreter','latex')
ylabel('screening selection tolerance ($\delta_2$)', 'interpreter','latex')
c = colorbar;
c.Label.String = 'Average Subset Size';
c.Label.Interpreter = 'latex';

% plot idea:

% x axis: delta 1 (constant overall delta)
% y axis: alpha 1 (constant overall alpha)
% z axis: average sample size used for selection

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
% the screening procedure used is ODP with marginal PGS

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


function Hstar = CRN_rinott_number(Z,Y,C,k,pstar,n0)
% calculate Rinott's number with CRN from a common sample
Cmat=repmat(C,1,k-1);
denom = sqrt((Y(:,1:(k-1)) + Cmat) .* (1/(n0-1)));
H= max(Z(:,1:(k-1))./ denom,[],2);
Hstar=quantile(H,pstar);
end
    
    