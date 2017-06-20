function [pos,tsq,k] = bisection_2(limit,initState,finalState,initTime)

smallCont = [];
bigCont = [];
k = 0;

goal = finalState(1);

[pos,tsq] = calcX5_2(limit,initState,finalState,initTime(1));
if pos < goal
    smallCont = [initTime(1);pos];
else
    bigCont = [initTime(1);pos];
end

[pos,tsq] = calcX5_2(limit,initState,finalState,initTime(2));
if pos < goal
    smallCont = [initTime(2);pos];
else
    bigCont = [initTime(2);pos];
end

if and(~isempty(smallCont),~isempty(bigCont))
    time = (initTime(1) + initTime(2))/2;
    [pos,tsq] = calcX5_2(limit,initState,finalState,time);
   
    if pos < goal
        smallCont = [smallCont [time;pos]];
    else
        bigCont = [bigCont [time;pos]];
    end
    
    for k = 1:100
        if (pos >= goal)
            time = (time + smallCont(1,end)) / 2;
        else
            time = (time + bigCont(1,end)) / 2;
        end
        [pos,tsq] = calcX5_2(limit,initState,finalState,time);
        
        if pos < goal
            smallCont = [smallCont [time;pos]];
        else
            bigCont = [bigCont [time;pos]];
        end

        if (norm(pos - goal) < 1e-3)
            break;
        end
        
        pos
    end
    
end
