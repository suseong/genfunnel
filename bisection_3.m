function [posOut,tsqOut] = bisection_3(limit,initState,finalState,idx)

smallCont = [];
bigCont = [];
k = 0;

goal = finalState(1);

timeTrials = 0.02:1:20;
posTrials = zeros(length(timeTrials),2);

for k=1:length(timeTrials)
   [posTrials(k,:),~] = calcX5_3(limit,initState,finalState,timeTrials(k)); 
end

chad = 0;
for k=1:length(timeTrials)-1
   if (posTrials(k,idx) - goal)*(posTrials(k+1,idx) - goal) <= 0
       chad = chad + 1;
      if posTrials(k,idx) < goal
          smallCont{chad} = [timeTrials(k);posTrials(k,idx)];
          bigCont{chad} = [timeTrials(k+1);posTrials(k+1,idx)];
      else
          bigCont{chad} = [timeTrials(k);posTrials(k,idx)];
          smallCont{chad} = [timeTrials(k+1);posTrials(k+1,idx)];
      end
   end
end

posOut = [];
tsqOut = [];

if and(~isempty(smallCont),~isempty(bigCont))
    for j = 1:length(smallCont)
        time = (smallCont{j}(1) + bigCont{j}(1))/2;
        [pos,~] = calcX5_3(limit,initState,finalState,time);
        
        if pos(idx) < goal
            smallCont{j} = [time;pos(idx)];
        else
            bigCont{j} = [time;pos(idx)];
        end
        
        for k = 1:200
            if (pos(idx) >= goal)
                time = (time + smallCont{j}(1)) / 2;
            else
                time = (time + bigCont{j}(1)) / 2;
            end
            [pos,tsq] = calcX5_3(limit,initState,finalState,time);
            
            if pos(idx) < goal
                smallCont{j} = [time;pos(idx)];
            else
                bigCont{j} = [time;pos(idx)];
            end
            if (norm(pos(idx) - goal) < 1e-4)
                posOut{j} = pos(idx);
                tsqOut{j} = tsq(idx,:);
                break;
            end
        end
    end
end

end
