macro_rep=1000;
record=zeros(macro_rep,1);
avg_reocrd=zeros(macro_rep,1);
for rep=1:macro_rep
    rep
    max_loc= Rinott_Selection(12, 0.05, 2, 0.5);
    if max_loc<3
        record(rep)=1;
    end

    avg_reocrd(rep)=mean(record(1:rep));
end
plot(avg_reocrd)
hold on
line([0,macro_rep],[0.95,0.95],'color','r')