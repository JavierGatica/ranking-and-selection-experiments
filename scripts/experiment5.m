%% Paramters

add_rm_paths('add');

n0 = 20; % samples per system
alpha = 0.05;
nproc = 10; % number of abstract procesors
k = 10000;


% system index assigned to each processor in the first screening stage
partition = zeros(1,nproc+1);
partition(1) = 1;
partition(nproc+1) = k+1;
for i = 2:nproc
    partition(i) = partition(i-1) + ceil(k/nproc);
end


variance = 0.01;
means = betainv(linspace(0,1,k),4,2);
M = 100;

total_sample = sim_output(means,variance,n0,1);

add_rm_paths('remove');

%% DSTTB + STTB:

proc_times = zeros(1,nproc);
subset = false(1,k);
t1 = zeros(1,M);
t2 = zeros(1,M);
for m = 1:M

    % stage 1: DSTTB over partition of systems
    for i = 1:nproc
        [proc_times(i), proc_subset] = timeDSTTB(total_sample(:,partition(i):(partition(i+1)-1)),alpha);
        subset(partition(i):(partition(i+1)-1)) = proc_subset;
        sum(proc_subset)/ceil(k/nproc);
    end
    t1(m) = max(proc_times); % representing the last batch of the partition to finish

    % stage 2: STTB over surviving systems
    t2(m) = timeSTTB(total_sample(:,subset),alpha);
end

t_DS_mean = mean(t1+t2);
t_DS_q10 = quantile(t1+t2,0.1);
t_DS_q90 = quantile(t1+t2,0.9);

[t_DS_q10, t_DS_mean, t_DS_q90]

%% STTB + STTB

proc_times = zeros(1,nproc);
subset = false(1,k);
t1 = zeros(1,M);
t2 = zeros(1,M);
for m = 1:M

    % stage 1: DSTTB over partition of systems
    for i = 1:nproc
        [proc_times(i), proc_subset] = timeSTTB(total_sample(:,partition(i):(partition(i+1)-1)),alpha);
        subset(partition(i):(partition(i+1)-1)) = proc_subset;
        sum(proc_subset)/ceil(k/nproc);
    end
    t1(m) = max(proc_times); % representing the last batch of the partition to finish

    % stage 2: STTB over surviving systems
    t2(m) = timeSTTB(total_sample(:,subset),alpha);
end

t_SS_mean = mean(t1+t2);
t_SS_q10 = quantile(t1+t2,0.1);
t_SS_q90 = quantile(t1+t2,0.9);

[t_SS_q10, t_SS_mean, t_SS_q90]


%% Procedure definitions

function [t,subset] = timeSTTB(sample,alpha)

tic;

[n0,k] = size(sample);

sample_means = mean(sample);
sample_variances = var(sample);

W = zeros(k,k);
t_value = tinv((1-alpha)^(1/(k-1)),n0-1);
for i = 1:k
    for j = (i+1):k
        temp = t_value * sqrt((sample_variances(i)+sample_variances(j))/n0);
        W(i,j) = temp;
        W(j,i) = temp;
    end
end

subset = false(1,k);
for i = 1:k
    subset(i) = ~any(sample_means(i) < sample_means - W(i,:));
end

t = toc;

end


function [t,subset] = timeDSTTB(sample,alpha)

tic;

[n0,k] = size(sample);

sample_means = mean(sample);
sample_variances = var(sample);

t_value = tinv((1-alpha)^(1/k),n0-1);
W = t_value * sqrt(sample_variances / n0);

subset = false(1,k);
temp = sample_means - W;
for i = 1:k
    subset(i) = ~any(sample_means(i) < temp - W(i));
end

t = toc;

end