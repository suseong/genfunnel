function [time_seq,pos,acc,iter,act,num,isOkay] = calc_mintime_traj(initState,finalState,actLimit)

actLimit1 = actLimit; actLimit2 = -actLimit;
numValidSol = 0;
%% case 1 t1 = t2, t3 = t4
infos = [];
tfList = [];
% pos / time sequence
[chad1,chad2] = bisection_1(actLimit1,initState,finalState); chad3 = [];
for k=1:length(chad1)
    [solution,chad3{k}] = check_sol(chad1{k},initState,finalState,chad2{k},actLimit1);
    if solution
        numValidSol = numValidSol + 1;
        infos{numValidSol,1} = chad1{k}; infos{numValidSol,2} = chad2{k}; infos{numValidSol,3} = chad3{k}; infos{numValidSol,4} = 10;
        tfList = [tfList;chad2{k}(end)];
    end
end

[chad1,chad2] = bisection_1(actLimit2,initState,finalState); chad3 = [];
for k=1:length(chad1)
    [solution,chad3{k}] = check_sol(chad1{k},initState,finalState,chad2{k},actLimit2);
    if solution
        numValidSol = numValidSol + 1;
        infos{numValidSol,1} = chad1{k}; infos{numValidSol,2} = chad2{k}; infos{numValidSol,3} = chad3{k}; infos{numValidSol,4} = 15;
        tfList = [tfList;chad2{k}(end)];
    end
end

%% case 2 t3 = t4
[chad1,chad2] = bisection_2(actLimit1,initState,finalState,1); chad3 = [];
for k=1:length(chad1)
    [solution,chad3{k}] = check_sol(chad1{k},initState,finalState,chad2{k},actLimit1);
    if solution
        numValidSol = numValidSol + 1;
        infos{numValidSol,1} = chad1{k}; infos{numValidSol,2} = chad2{k}; infos{numValidSol,3} = chad3{k}; infos{numValidSol,4} = 20;
        tfList = [tfList;chad2{k}(end)];        
    end
end

[chad1,chad2] = bisection_2(actLimit1,initState,finalState,2); chad3 = [];
for k=1:length(chad1)
    [solution,chad3{k}] = check_sol(chad1{k},initState,finalState,chad2{k},actLimit1);
    if solution
        numValidSol = numValidSol + 1;
        infos{numValidSol,1} = chad1{k}; infos{numValidSol,2} = chad2{k}; infos{numValidSol,3} = chad3{k}; infos{numValidSol,4} = 22.5;
        tfList = [tfList;chad2{k}(end)];        
    end
end

[chad1,chad2] = bisection_2(actLimit2,initState,finalState,1); chad3 = [];
for k=1:length(chad1)
    [solution,chad3{k}] = check_sol(chad1{k},initState,finalState,chad2{k},actLimit2);
    if solution
        numValidSol = numValidSol + 1;
        infos{numValidSol,1} = chad1{k}; infos{numValidSol,2} = chad2{k}; infos{numValidSol,3} = chad3{k}; infos{numValidSol,4} = 25;
        tfList = [tfList;chad2{k}(end)];        
    end
end

[chad1,chad2] = bisection_2(actLimit2,initState,finalState,2); chad3 = [];
for k=1:length(chad1)
    [solution,chad3{k}] = check_sol(chad1{k},initState,finalState,chad2{k},actLimit2);
    if solution
        numValidSol = numValidSol + 1;
        infos{numValidSol,1} = chad1{k}; infos{numValidSol,2} = chad2{k}; infos{numValidSol,3} = chad3{k}; infos{numValidSol,4} = 27.5;
        tfList = [tfList;chad2{k}(end)];        
    end
end

%% case 3 t1 = t2
[chad1,chad2] = bisection_3(actLimit1,initState,finalState,1); chad3 = [];
for k=1:length(chad1)
    [solution,chad3{k}] = check_sol(chad1{k},initState,finalState,chad2{k},actLimit1);
    if solution
        numValidSol = numValidSol + 1;
        infos{numValidSol,1} = chad1{k}; infos{numValidSol,2} = chad2{k}; infos{numValidSol,3} = chad3{k}; infos{numValidSol,4} = 30;
        tfList = [tfList;chad2{k}(end)];        
    end
end

[chad1,chad2] = bisection_3(actLimit1,initState,finalState,2); chad3 = [];
for k=1:length(chad1)
    [solution,chad3{k}] = check_sol(chad1{k},initState,finalState,chad2{k},actLimit1);
    if solution
        numValidSol = numValidSol + 1;
        infos{numValidSol,1} = chad1{k}; infos{numValidSol,2} = chad2{k}; infos{numValidSol,3} = chad3{k}; infos{numValidSol,4} = 32.5;
        tfList = [tfList;chad2{k}(end)];        
    end
end

[chad1,chad2] = bisection_3(actLimit2,initState,finalState,1); chad3 = [];
for k=1:length(chad1)
    [solution,chad3{k}] = check_sol(chad1{k},initState,finalState,chad2{k},actLimit2);
    if solution
        numValidSol = numValidSol + 1;
        infos{numValidSol,1} = chad1{k}; infos{numValidSol,2} = chad2{k}; infos{numValidSol,3} = chad3{k}; infos{numValidSol,4} = 35;
        tfList = [tfList;chad2{k}(end)];        
    end
end

[chad1,chad2] = bisection_3(actLimit2,initState,finalState,2); chad3 = [];
for k=1:length(chad1)
    [solution,chad3{k}] = check_sol(chad1{k},initState,finalState,chad2{k},actLimit2);
    if solution
        numValidSol = numValidSol + 1;
        infos{numValidSol,1} = chad1{k}; infos{numValidSol,2} = chad2{k}; infos{numValidSol,3} = chad3{k}; infos{numValidSol,4} = 37.5;
        tfList = [tfList;chad2{k}(end)];        
    end
end

%% case 4
[chad1,chad2] = bisection_4(actLimit1,initState,finalState); chad3 = [];
for k=1:length(chad1)
    [solution,chad3{k}] = check_sol(chad1{k},initState,finalState,chad2{k},actLimit1);
    if solution
        numValidSol = numValidSol + 1;
        infos{numValidSol,1} = chad1{k}; infos{numValidSol,2} = chad2{k}; infos{numValidSol,3} = chad3{k}; infos{numValidSol,4} = 40;
        tfList = [tfList;chad2{k}(end)];        
    end
end

[chad1,chad2] = bisection_4(actLimit2,initState,finalState); chad3 = [];
for k=1:length(chad1)
    [solution,chad3{k}] = check_sol(chad1{k},initState,finalState,chad2{k},actLimit2);
    if solution
        numValidSol = numValidSol + 1;
        infos{numValidSol,1} = chad1{k}; infos{numValidSol,2} = chad2{k}; infos{numValidSol,3} = chad3{k}; infos{numValidSol,4} = 45;
        tfList = [tfList;chad2{k}(end)];        
    end
end

%%
% num = 1;
if size(infos,1) >= 1
%     num = size(infos,1);
    isOkay = 1;
else
% elseif size(infos,1) == 0
    isOkay = 0;
%     keyboard;
end

if isOkay == 1
    [~,k] = min(tfList);
    pos = infos{k,1};
    iter = 0;
    acc = infos{k,3};
    act = infos{k,4};
    num = 2*(floor(mod(act,10)/5)-0.5);
    time_seq = infos{k,2};
else
    pos = [];
    iter = 0;
    acc = [];
    act = [];
    num = 3;
    time_seq = [];    
end
    
end