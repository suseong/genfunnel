function [pos,tsq,k] = bisection_4(limit,initState,finalState,initTime)

smallCont = [];
bigCont = [];
k = 0;

goal = finalState(1);

timeTrials = 0.2:0.5:10;
posTrials = zeros(1,length(timeTrials));

for k=1:length(timeTrials)
   [posTrials(k),~] = calcX5_4(limit,initState,finalState,timeTrials(k)); 
end

chad = 0;
for k=1:length(timeTrials)-1
   if (posTrials(k) - goal)*(posTrials(k+1) - goal) <= 0
       chad = chad + 1;
      if posTrials(k) < goal
          smallCont{chad} = [timeTrials(k);posTrials(k)];
          bigCont{chad} = [timeTrials(k+1);posTrials(k+1)];
      else
          smallCont{chad} = [timeTrials(k);posTrials(k)];
          bigCont{chad} = [timeTrials(k+1);posTrials(k+1)];
      end
   end
end

[pos,tsq] = calcX5_4(limit,initState,finalState,0.1);

if pos < goal
    smallCont = [initTime(1);pos];
else
    bigCont = [initTime(1);pos];
end

for k=1:30
   initTime = 10*rand(1,1);
   [pos,tsq] = calcX5_4(limit,initState,finalState,initTime);
   if pos < goal
       smallCont = [initTime;pos];
   else
       bigCont = [initTime;pos];
   end
   
   if and(~isempty(smallCont),~isempty(bigCont)) 
       break;
   end
end

if and(~isempty(smallCont),~isempty(bigCont))
    time = (smallCont(1) + bigCont(1))/2;
    [pos,tsq] = calcX5_4(limit,initState,finalState,time);
   
    if pos < goal
        smallCont = [time;pos];
    else
        bigCont = [time;pos];
    end
    
    for k = 1:100
        if (pos >= goal)
            time = (time + smallCont(1,end)) / 2;
        else
            time = (time + bigCont(1,end)) / 2;
        end
        [pos,tsq] = calcX5_4(limit,initState,finalState,time);
        
        if pos < goal
            smallCont = [time;pos];
        else
            bigCont = [time;pos];
        end
        if (norm(pos - goal) < 1e-3)
            break;
        end
    end
    
end
