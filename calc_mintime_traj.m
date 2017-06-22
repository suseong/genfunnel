function [time_seq,pos,acc,iter,act] = calc_mintime_traj(initState,finalState,actLimit)

actLimit1 = actLimit; actLimit2 = -actLimit;

%% case 1 t1 = t2, t3 = t4
initTime = [0.1 20];

infos = [];
% pos / time sequence / interation number
k=1;[infos{k,1},infos{k,2},infos{k,3}] = bisection_1(actLimit1,initState,finalState,initTime);

initTime = [0.1 20];
k=2;[infos{k,1},infos{k,2},infos{k,3}] = bisection_1(actLimit2,initState,finalState,initTime);

%% case 2 t3 = t4
initTime = [0.1 20];
k=3;[infos{k,1},infos{k,2},infos{k,3}] = bisection_2(actLimit1,initState,finalState,initTime,1);
k=4;[infos{k,1},infos{k,2},infos{k,3}] = bisection_2(actLimit1,initState,finalState,initTime,2);

initTime = [0.1 20];
k=5;[infos{k,1},infos{k,2},infos{k,3}] = bisection_2(actLimit2,initState,finalState,initTime,1);
k=6;[infos{k,1},infos{k,2},infos{k,3}] = bisection_2(actLimit2,initState,finalState,initTime,2);

%% case 3 t1 = t2
initTime = [0.1 20];
k=7;[infos{k,1},infos{k,2},infos{k,3}] = bisection_3(actLimit1,initState,finalState,initTime,1);
k=8;[infos{k,1},infos{k,2},infos{k,3}] = bisection_3(actLimit1,initState,finalState,initTime,2);

initTime = [0.1 20];
k=9;[infos{k,1},infos{k,2},infos{k,3}] = bisection_3(actLimit2,initState,finalState,initTime,1);
k=10;[infos{k,1},infos{k,2},infos{k,3}] = bisection_3(actLimit2,initState,finalState,initTime,2);

%% case 4
initTime = [0.1 20];
k=11;[infos{k,1},infos{k,2},infos{k,3}] = bisection_4(actLimit1,initState,finalState,initTime);

initTime = [0.1 20];
k=12;[infos{k,1},infos{k,2},infos{k,3}] = bisection_4(actLimit2,initState,finalState,initTime);

%%
solution = zeros(1,12);

k=1;  [solution(k),infos{k,4}] = check_sol(infos{k,1},initState,finalState,infos{k,2},actLimit1);
k=2;  [solution(k),infos{k,4}] = check_sol(infos{k,1},initState,finalState,infos{k,2},actLimit2);

k=3;  [solution(k),infos{k,4}] = check_sol(infos{k,1},initState,finalState,infos{k,2},actLimit1);
k=4;  [solution(k),infos{k,4}] = check_sol(infos{k,1},initState,finalState,infos{k,2},actLimit1);
k=5;  [solution(k),infos{k,4}] = check_sol(infos{k,1},initState,finalState,infos{k,2},actLimit2);
k=6;  [solution(k),infos{k,4}] = check_sol(infos{k,1},initState,finalState,infos{k,2},actLimit2);

k=7;  [solution(k),infos{k,4}] = check_sol(infos{k,1},initState,finalState,infos{k,2},actLimit1);
k=8;  [solution(k),infos{k,4}] = check_sol(infos{k,1},initState,finalState,infos{k,2},actLimit1);
k=9;  [solution(k),infos{k,4}] = check_sol(infos{k,1},initState,finalState,infos{k,2},actLimit2);
k=10; [solution(k),infos{k,4}] = check_sol(infos{k,1},initState,finalState,infos{k,2},actLimit2);

k=11; [solution(k),infos{k,4}] = check_sol(infos{k,1},initState,finalState,infos{k,2},actLimit1);
k=12; [solution(k),infos{k,4}] = check_sol(infos{k,1},initState,finalState,infos{k,2},actLimit2);

%%
k = find(solution == 1);

if length(k) > 1
    if infos{k(1),2}(end) < infos{k(2),2}(end)
        k = k(1);
    else
        k = k(2);
    end
end

if length(k) > 2
    keyboard
end

try
    time_seq = infos{k,2};
catch
    keyboard
end

pos = infos{k,1};
iter = infos{k,3};
acc = infos{k,4};

try
    if sum(k == [1 3 4 7 8 11])
        act = 1;
    else
        act = -1;
    end 
catch
    keyboard
end

end