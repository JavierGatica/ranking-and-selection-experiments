function max_loc= Rinott_Selection(k, alpha, n0, delta)
%  implements Rinott's procedure
%  k = number of systems
%  1-alpha = desired PCS
%  n0 = first-stage sample size
%  delta = indifference-zone parameter
%  note: uses 99% UCB for Rinott's h
sys_mean=inf*ones(k,1);
sys_var=inf*ones(k,1);
sys_sample_size=inf*ones(k,1);
[~,h]= Rinott_Number(k, n0, 1-alpha, 0.99, 10000);%using UCB

for rep=1:k
    input=rep;%For special Mysim
    Simulation_output=Mysim(input,n0);
    Sample_Var = var(Simulation_output);
    N =ceil((h^2)*Sample_Var/(delta^2));
    if N>n0
       Simulation_output=[Simulation_output;Mysim(input,N-n0)];
    end
    sys_mean(rep,1)=mean(Simulation_output);
    sys_var(rep,1)=var(Simulation_output);
    sys_sample_size(rep,1)=max([N,n0]);
    [~,max_loc] = max(sys_mean);
end

