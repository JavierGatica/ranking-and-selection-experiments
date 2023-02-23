function rand_nums=Mysim(input,num)
    k=-[[0 0.5],ones(1,10)];
    rand_nums=randn(num,1)+k(input);%Generating normal dustribued random number with mean=k(input), var=1
end
