function [time_seq,pos,acc,iter,act,num,isOkay] = calc_mintime_traj_(initState,finalState,actLimit)

%%        1-1   1-2   1-3 | 2-1         2-2   2-3       | 3-1   3-2   3-3         3-4           | 4-1     4-2     4-3     4-4 
input_ = [1  1  2  2  3  3  4  4  4  4  5  5  6  6  6  6  7  7  8  8  9  9  9  9  10  10  10  10  11  11  12  12  13  13  14  14;
          1 -1  1 -1  1 -1  1  1 -1 -1  1 -1  1  1 -1 -1  1 -1  1 -1  1  1 -1 -1   1   1  -1  -1   1  -1   1  -1   1  -1   1  -1; 
          1  1  1  1  1  1  1  2  1  2  1  1  1  2  1  2  1  1  1  1  1  2  1  2   1   2   1   2   1   1   1   1   1   1   1   1];
%         1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24  25  26  27  28  29  30  31  32  33  34  35  36
        
infos = [];
tfList = [];
% pos / time sequence

actLimit = [20 1];

clc
% init = [10*(rand(1,3)-0.5)] 
% final = [10*(rand(1,3)-0.5)]

init = [ 3.0033   -0.4620   -0.6761];
final = [ 3.2531   -4.1653   -3.6683];

numValidSol = 0;

for j = 1:size(input_,2)
    [chad1,chad2,chad3] = bisection_(input_(2,j)*actLimit,init,final,input_(3,j),input_(1,j));
    if ~isempty(chad1)
        j
        chad1{1}
        chad2{1}
        chad3{1}
        disp('-----------')
    end
    
    for k=1:length(chad1)
%         chad1{k}
%         chad2{k}
%         chad3{k}
        if min(diff(chad2{k})) > -1e-3 && chad2{k}(1) > 0 && max(abs(chad3{k}(1:4))) < actLimit(2)*1.01
            numValidSol = numValidSol + 1;
            infos{numValidSol,1} = chad1{k}; infos{numValidSol,2} = chad2{k}; infos{numValidSol,3} = chad3{k}; infos{numValidSol,4} = j;
            tfList = [tfList;chad2{k}(end)];
            chad1{k}
            chad2{k}
            chad3{k}
            disp('-----------')
        end
    end
    
end
%%

% [chad1,chad2] = bisection_1(actLimit1,initState,finalState); chad3 = [];
% for k=1:length(chad1)
%     [solution,chad3{k}] = check_sol(chad1{k},initState,finalState,chad2{k},actLimit1);
%     if solution
%         numValidSol = numValidSol + 1;
%         infos{numValidSol,1} = chad1{k}; infos{numValidSol,2} = chad2{k}; infos{numValidSol,3} = chad3{k}; infos{numValidSol,4} = 10;
%         tfList = [tfList;chad2{k}(end)];
%     end
% end
% 
% [chad1,chad2] = bisection_1(actLimit2,initState,finalState); chad3 = [];
% for k=1:length(chad1)
%     [solution,chad3{k}] = check_sol(chad1{k},initState,finalState,chad2{k},actLimit2);
%     if solution
%         numValidSol = numValidSol + 1;
%         infos{numValidSol,1} = chad1{k}; infos{numValidSol,2} = chad2{k}; infos{numValidSol,3} = chad3{k}; infos{numValidSol,4} = 15;
%         tfList = [tfList;chad2{k}(end)];
%     end
% end
% 
% %% case 2 t3 = t4
% [chad1,chad2] = bisection_2(actLimit1,initState,finalState,1); chad3 = [];
% for k=1:length(chad1)
%     [solution,chad3{k}] = check_sol(chad1{k},initState,finalState,chad2{k},actLimit1);
%     if solution
%         numValidSol = numValidSol + 1;
%         infos{numValidSol,1} = chad1{k}; infos{numValidSol,2} = chad2{k}; infos{numValidSol,3} = chad3{k}; infos{numValidSol,4} = 20;
%         tfList = [tfList;chad2{k}(end)];        
%     end
% end
% 
% [chad1,chad2] = bisection_2(actLimit1,initState,finalState,2); chad3 = [];
% for k=1:length(chad1)
%     [solution,chad3{k}] = check_sol(chad1{k},initState,finalState,chad2{k},actLimit1);
%     if solution
%         numValidSol = numValidSol + 1;
%         infos{numValidSol,1} = chad1{k}; infos{numValidSol,2} = chad2{k}; infos{numValidSol,3} = chad3{k}; infos{numValidSol,4} = 22.5;
%         tfList = [tfList;chad2{k}(end)];        
%     end
% end
% 
% [chad1,chad2] = bisection_2(actLimit2,initState,finalState,1); chad3 = [];
% for k=1:length(chad1)
%     [solution,chad3{k}] = check_sol(chad1{k},initState,finalState,chad2{k},actLimit2);
%     if solution
%         numValidSol = numValidSol + 1;
%         infos{numValidSol,1} = chad1{k}; infos{numValidSol,2} = chad2{k}; infos{numValidSol,3} = chad3{k}; infos{numValidSol,4} = 25;
%         tfList = [tfList;chad2{k}(end)];        
%     end
% end
% 
% [chad1,chad2] = bisection_2(actLimit2,initState,finalState,2); chad3 = [];
% for k=1:length(chad1)
%     [solution,chad3{k}] = check_sol(chad1{k},initState,finalState,chad2{k},actLimit2);
%     if solution
%         numValidSol = numValidSol + 1;
%         infos{numValidSol,1} = chad1{k}; infos{numValidSol,2} = chad2{k}; infos{numValidSol,3} = chad3{k}; infos{numValidSol,4} = 27.5;
%         tfList = [tfList;chad2{k}(end)];        
%     end
% end
% 
% %% case 3 t1 = t2
% [chad1,chad2] = bisection_3(actLimit1,initState,finalState,1); chad3 = [];
% for k=1:length(chad1)
%     [solution,chad3{k}] = check_sol(chad1{k},initState,finalState,chad2{k},actLimit1);
%     if solution
%         numValidSol = numValidSol + 1;
%         infos{numValidSol,1} = chad1{k}; infos{numValidSol,2} = chad2{k}; infos{numValidSol,3} = chad3{k}; infos{numValidSol,4} = 30;
%         tfList = [tfList;chad2{k}(end)];        
%     end
% end
% 
% [chad1,chad2] = bisection_3(actLimit1,initState,finalState,2); chad3 = [];
% for k=1:length(chad1)
%     [solution,chad3{k}] = check_sol(chad1{k},initState,finalState,chad2{k},actLimit1);
%     if solution
%         numValidSol = numValidSol + 1;
%         infos{numValidSol,1} = chad1{k}; infos{numValidSol,2} = chad2{k}; infos{numValidSol,3} = chad3{k}; infos{numValidSol,4} = 32.5;
%         tfList = [tfList;chad2{k}(end)];        
%     end
% end
% 
% [chad1,chad2] = bisection_3(actLimit2,initState,finalState,1); chad3 = [];
% for k=1:length(chad1)
%     [solution,chad3{k}] = check_sol(chad1{k},initState,finalState,chad2{k},actLimit2);
%     if solution
%         numValidSol = numValidSol + 1;
%         infos{numValidSol,1} = chad1{k}; infos{numValidSol,2} = chad2{k}; infos{numValidSol,3} = chad3{k}; infos{numValidSol,4} = 35;
%         tfList = [tfList;chad2{k}(end)];        
%     end
% end
% 
% [chad1,chad2] = bisection_3(actLimit2,initState,finalState,2); chad3 = [];
% for k=1:length(chad1)
%     [solution,chad3{k}] = check_sol(chad1{k},initState,finalState,chad2{k},actLimit2);
%     if solution
%         numValidSol = numValidSol + 1;
%         infos{numValidSol,1} = chad1{k}; infos{numValidSol,2} = chad2{k}; infos{numValidSol,3} = chad3{k}; infos{numValidSol,4} = 37.5;
%         tfList = [tfList;chad2{k}(end)];        
%     end
% end
% 
% %% case 4
% [chad1,chad2] = bisection_4(actLimit1,initState,finalState); chad3 = [];
% for k=1:length(chad1)
%     [solution,chad3{k}] = check_sol(chad1{k},initState,finalState,chad2{k},actLimit1);
%     if solution
%         numValidSol = numValidSol + 1;
%         infos{numValidSol,1} = chad1{k}; infos{numValidSol,2} = chad2{k}; infos{numValidSol,3} = chad3{k}; infos{numValidSol,4} = 40;
%         tfList = [tfList;chad2{k}(end)];        
%     end
% end
% 
% [chad1,chad2] = bisection_4(actLimit2,initState,finalState); chad3 = [];
% for k=1:length(chad1)
%     [solution,chad3{k}] = check_sol(chad1{k},initState,finalState,chad2{k},actLimit2);
%     if solution
%         numValidSol = numValidSol + 1;
%         infos{numValidSol,1} = chad1{k}; infos{numValidSol,2} = chad2{k}; infos{numValidSol,3} = chad3{k}; infos{numValidSol,4} = 45;
%         tfList = [tfList;chad2{k}(end)];        
%     end
% end
% 
% %%
% % num = 1;
% if size(infos,1) >= 1
% %     num = size(infos,1);
%     isOkay = 1;
% else
% % elseif size(infos,1) == 0
%     isOkay = 0;
% %     keyboard;
% end
% 
% if isOkay == 1
%     [~,k] = min(tfList);
%     pos = infos{k,1};
%     iter = 0;
%     acc = infos{k,3};
%     act = infos{k,4};
%     num = 2*(floor(mod(act,10)/5)-0.5);
%     time_seq = infos{k,2};
% else
%     pos = [];
%     iter = 0;
%     acc = [];
%     act = [];
%     num = 3;
%     time_seq = [];    
% end
    
end