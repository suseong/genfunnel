function [pos_,tsq_,k] = bisection_2(limit,initState,finalState,initTime,idx)

smallCont = [];
bigCont = [];
pos_ = [];
tsq_ = [];
k = 0;

goal = finalState(1);

[pos,tsq] = calcX5_2(limit,initState,finalState,0.1);

if pos(idx) < goal
    smallCont = [initTime(1);pos(idx)];
else
    bigCont = [initTime(1);pos(idx)];
end
% 
% [pos,tsq] = calcX5_2(limit,initState,finalState,initTime(2));
% if pos(idx) < goal
%     smallCont = [initTime(2);pos(idx)];
% else
%     bigCont = [initTime(2);pos(idx)];
% end
% 
% if and(~isempty(smallCont),~isempty(bigCont))
%     time = (initTime(1) + initTime(2))/2;

for k=1:30
   initTime = 10*rand(1,1);
   [pos,tsq] = calcX5_2(limit,initState,finalState,initTime);
   if pos < goal
       smallCont = [initTime;pos(idx)];
   else
       bigCont = [initTime;pos(idx)];
   end
   
   if and(~isempty(smallCont),~isempty(bigCont)) 
       break;
   end
end

if and(~isempty(smallCont),~isempty(bigCont))
    time = (smallCont(1) + bigCont(1))/2;
    [pos,tsq] = calcX5_2(limit,initState,finalState,time);
   
    if pos(idx) < goal
        smallCont = [smallCont [time;pos(idx)]];
    else
        bigCont = [bigCont [time;pos(idx)]];
    end
    
    for k = 1:100
        if (pos(idx) >= goal)
            time = (time + smallCont(1,end)) / 2;
        else
            time = (time + bigCont(1,end)) / 2;
        end
        [pos,tsq] = calcX5_2(limit,initState,finalState,time);
        
        if pos(idx) < goal
            smallCont = [smallCont [time;pos(idx)]];
        else
            bigCont = [bigCont [time;pos(idx)]];
        end

        if (norm(pos(idx) - goal) < 1e-3)
            break;
        end
    end
    
    pos_ = pos(idx);
    tsq_ = tsq(idx,:);
end

end
