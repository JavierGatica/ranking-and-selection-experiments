function elim = BiPASS(simulation,c,n0,dn,Nmax,k)
g = @(t) sqrt((c + log(t+1))*(t+1));
II = 1:k;
active = true(1,k);
elim = zeros(1,k);
Yn0 = simulation(n0, active);
S2 = mean(var(Yn0));
Ysum = sum(Yn0);
r = n0;
N = n0*k;
while sum(active) > 1 && N < Nmax
    r = r + dn;
    N = N + dn*sum(active);
    Ynew = simulation(dn, active);
    Ysum(active) = Ysum(active) + sum(Ynew,1);
    rmuhat = mean(Ysum(active));
    for l = II(active)
        if Ysum(l) - rmuhat < -g(r/S2)*S2
            active(l) = false;
            elim(l) = N;
        end
    end
end
end

