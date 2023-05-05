%% Paramters

add_rm_paths('add');

n0 = 20; 
alpha = 0.05;
nproc = [10, 20, 50, 100]; 
k = 100000;


variance = 0.1;
means = betainv(linspace(0,1,k),2,4);
means = means(randperm(k));
delta = 1 - betainv(0.9,2,4);
M = 20;


subset = false(1,k);
t1 = zeros(1,M);
t2 = zeros(1,M);


t_mS = zeros(1,M);
t_mDS = zeros(1,M);
t_SmS = zeros(length(nproc),M);
t_DmS = zeros(length(nproc),M);
t_DmDS = zeros(length(nproc),M);
t_DSmDS = zeros(length(nproc),M);


subset_size = zeros(length(nproc), 6);

rng(1);

proc_times = zeros(1,nproc(p));

% %% | STTB:
% 
% "| STTB"
% parfor m = 1:M
%     total_sample = sim_output(means,variance,n0,1);
%     [t2(m),subset] = timeSTTB(total_sample, alpha,delta);
%     m
% end
% 
% t_mS = t2;
% 
% 
% %% | DSTTB + STTB:
% 
% "| DSTTB + STTB"
% parfor m = 1:M
%     total_sample = sim_output(means,variance,n0,1);
%     [t2(m), subset] = timeDSTTB(total_sample, alpha, delta);
%     [temp, subset2] = timeSTTB(total_sample(:,subset),alpha,delta);
%     t2(m) = t2(m) + temp;
%     m
% end
% 
% t_mDS = t2;


for p = 1:length(nproc)

partition = partition_indices(k, nproc(p));

%% STTB | STTB:

"STTB | STTB"
for m = 1:M
    total_sample = sim_output(means,variance,n0,1);
    for i = 1:nproc(p)
	    indices = partition(i):(partition(i+1)-1);
        [proc_times(i), proc_subset] = timeSTTB(total_sample(:,indices),alpha,delta);
        subset(indices) = proc_subset;
    end
    t1(m) = max(proc_times);
    [t2(m), subset2] = timeSTTB(total_sample(:,subset),alpha,delta);
    subset(subset) = subset2;
    m
end

t_SmS(p,:) = t1 + t2;


%% DSTTB | STTB:

"DSTTB | STTB"
for m = 1:M
    total_sample = sim_output(means,variance,n0,1);
    for i = 1:nproc(p)
	    indices = partition(i):(partition(i+1)-1);
        [proc_times(i), proc_subset] = timeDSTTB(total_sample(:,indices),alpha,delta);
        subset(indices) = proc_subset;
    end
    t1(m) = max(proc_times);
    [t2(m), subset2] = timeSTTB(total_sample(:,subset),alpha,delta);
    subset(subset) = subset2;
    m
end

t_DmS(p,:) = t1 + t2;


%% DSTTB | DSTTB + STTB:

"DSTTB | DSTTB + STTB"
for m = 1:M
    total_sample = sim_output(means,variance,n0,1);
    for i = 1:nproc(p)
	    indices = partition(i):(partition(i+1)-1);
        [proc_times(i), proc_subset] = timeDSTTB(total_sample(:,indices),alpha,delta);
        subset(indices) = proc_subset;
    end
    t1(m) = max(proc_times); 
    [t2(m), subset2] = timeDSTTB(total_sample(:,subset), alpha, delta);
    subset(subset) = subset2;
    [temp, subset2] = timeSTTB(total_sample(:,subset),alpha,delta);
    subset(subset) = subset2;
    t2(m) = t2(m) + temp;
    m
end

t_DmDS(p,:) = t1 + t2;



%% DSTTB + STTB | DSTTB + STTB:

"DSTTB + STTB | DSTTB + STTB"
for m = 1:M
    total_sample = sim_output(means,variance,n0,1);
    for i = 1:nproc(p)
	    indices = partition(i):(partition(i+1)-1);
        [proc_times(i), proc_subset] = timeDSTTB(total_sample(:,indices),alpha,delta);
        subset(indices) = proc_subset;
	    [temp, proc_subset2] = timeSTTB(total_sample(:,indices(proc_subset)),alpha,delta);
        proc_times(i) = proc_times(i) + temp;
	    subset(indices(proc_subset)) = proc_subset2;
    end
    t1(m) = max(proc_times);
    [t2(m), subset2] = timeDSTTB(total_sample(:,subset),alpha,delta);
    subset(subset) = subset2;
    [temp, subset2] = timeSTTB(total_sample(:,subset),alpha,delta);
    subset(subset) = subset2;
    t2(m) = t2(m) + temp;
    m
end

t_DSmDS(p,:) = t1 + t2;

end


add_rm_paths('remove');

cd ..
save("data/exp5.mat")

%% Procedure definitions

function [t,subset] = timeSTTB(sample,alpha,delta)

tic;

[n0,k] = size(sample);

sample_means = mean(sample);
sample_variances = var(sample);

t_value = tinv((1-alpha)^(1/(k-1)),n0-1);

subset = false(1,k);
for i = 1:k
    subset(i) = ~any(sample_means(i) < sample_means - t_value * sqrt((sample_variances(i) + sample_variances)/n0) - delta);
end

t = toc;

end


function [t,subset] = timeDSTTB(sample,alpha,delta)

tic;

[n0,k] = size(sample);

sample_means = mean(sample);
sample_variances = var(sample);

t_value = tinv((1-alpha)^(1/k),n0-1);
W = t_value * sqrt(sample_variances / n0);

subset = false(1,k);
temp = sample_means - W;
for i = 1:k
    subset(i) = ~any(sample_means(i) < temp - W(i) - delta);
end

t = toc;

end

function ranges = partition_indices(k, p)

    ranges = zeros(1,p+1);
    ranges(1) = 1;
    ranges(end) = k + 1;
    for i = 2:p
        ranges(i) = ranges(i-1) + ceil(k/p);
    end


end
