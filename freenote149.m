clear all
close all
% clc

initState = [1;-1;0]; finalState = [-1;2;0];
actLimit1 = [10;6];  actLimit2 = [-10;-6];

%% case 1 t1 = t2, t3 = t4
initTime = [0.01 30];

infos = [];
% pos / time sequence / interation number
k=1;[infos{k,1},infos{k,2},infos{k,3}] = bisection_1(actLimit1,initState,finalState,initTime);

initTime = [0.01 30];
k=2;[infos{k,1},infos{k,2},infos{k,3}] = bisection_1(actLimit2,initState,finalState,initTime);

%% case 2 t3 = t4
initTime = [0.01 30];
k=3;[infos{k,1},infos{k,2},infos{k,3}] = bisection_2(actLimit1,initState,finalState,initTime,1);
k=4;[infos{k,1},infos{k,2},infos{k,3}] = bisection_2(actLimit1,initState,finalState,initTime,2);

initTime = [0.01 30];
k=5;[infos{k,1},infos{k,2},infos{k,3}] = bisection_2(actLimit2,initState,finalState,initTime,1);
k=6;[infos{k,1},infos{k,2},infos{k,3}] = bisection_2(actLimit2,initState,finalState,initTime,2);

%% case 3 t1 = t2
initTime = [0.01 30];
k=7;[infos{k,1},infos{k,2},infos{k,3}] = bisection_3(actLimit1,initState,finalState,initTime,1);
k=8;[infos{k,1},infos{k,2},infos{k,3}] = bisection_3(actLimit1,initState,finalState,initTime,2);

initTime = [0.01 30];
k=9;[infos{k,1},infos{k,2},infos{k,3}] = bisection_3(actLimit2,initState,finalState,initTime,1);
k=10;[infos{k,1},infos{k,2},infos{k,3}] = bisection_3(actLimit2,initState,finalState,initTime,2);

%% case 4
initTime = [0.01 30];
k=11;[infos{k,1},infos{k,2},infos{k,3}] = bisection_4(actLimit1,initState,finalState,initTime);

initTime = [0.01 30];
k=12;[infos{k,1},infos{k,2},infos{k,3}] = bisection_4(actLimit2,initState,finalState,initTime);

%%
solution = zeros(1,12);

k=1;  solution(k) = check_sol(infos{k,1},initState,finalState,infos{k,2},actLimit1);
k=2;  solution(k) = check_sol(infos{k,1},initState,finalState,infos{k,2},actLimit2);

k=3;  solution(k) = check_sol(infos{k,1},initState,finalState,infos{k,2},actLimit1);
k=4;  solution(k) = check_sol(infos{k,1},initState,finalState,infos{k,2},actLimit1);
k=5;  solution(k) = check_sol(infos{k,1},initState,finalState,infos{k,2},actLimit2);
k=6;  solution(k) = check_sol(infos{k,1},initState,finalState,infos{k,2},actLimit2);

k=7;  solution(k) = check_sol(infos{k,1},initState,finalState,infos{k,2},actLimit1);
k=8;  solution(k) = check_sol(infos{k,1},initState,finalState,infos{k,2},actLimit1);
k=9;  solution(k) = check_sol(infos{k,1},initState,finalState,infos{k,2},actLimit2);
k=10; solution(k) = check_sol(infos{k,1},initState,finalState,infos{k,2},actLimit2);

k=11; solution(k) = check_sol(infos{k,1},initState,finalState,infos{k,2},actLimit1);
k=12; solution(k) = check_sol(infos{k,1},initState,finalState,infos{k,2},actLimit2);

