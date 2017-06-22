pos1 = [];

for k=1:100
    [pos1(k),~] = calcX5_4(actLimit2,initState,finalState,k*0.1)
end
