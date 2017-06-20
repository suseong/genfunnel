clear all
close all
clc

initState = [0;0;0]; finalState = [5;2;0];
actLimit1 = [10;5];  actLimit2 = [-10;-5];

%% case 1 t1 = t2, t3 = t4
initTime = [0.01 30]; 
[pos11,tsq11,k11] = bisection_1(actLimit1,initState,finalState,initTime);

initTime = [0.01 30];
[pos12,tsq12,k12] = bisection_1(actLimit2,initState,finalState,initTime);

%% case 2 t3 = t4
initTime = [0.01 3];
[pos21,tsq21,k21] = bisection_2(actLimit1,initState,finalState,initTime);

initTime = [0.01 3];
[pos22,tsq22,k22] = bisection_2(actLimit2,initState,finalState,initTime);

%% case 3 t1 = t2
initTime = [0.01 30];
[pos31,tsq31,k31] = bisection_3(actLimit1,initState,finalState,initTime);

initTime = [0.01 30];
[pos32,tsq32,k32] = bisection_3(actLimit2,initState,finalState,initTime);

%% case 4
initTime = [0.01 30];
[pos41,tsq41,k41] = bisection_4(actLimit1,initState,finalState,initTime);

initTime = [0.01 30];
[pos42,tsq42,k42] = bisection_4(actLimit2,initState,finalState,initTime);

%%
