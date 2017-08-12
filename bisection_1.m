function [posOut,tsqOut] = bisection_1(limit,initState,finalState)

smallCont = [];
bigCont = [];
k = 0;

goal = finalState(1);

timeTrials = 0.02:1:20;
posTrials = zeros(1,length(timeTrials));

for k=1:length(timeTrials)
   [posTrials(k),~] = calcX5_1(limit,initState,finalState,timeTrials(k)); 
end

chad = 0;
for k=1:length(timeTrials)-1
   if (posTrials(k) - goal)*(posTrials(k+1) - goal) <= 0
       chad = chad + 1;
      if posTrials(k) < goal
          smallCont{chad} = [timeTrials(k);posTrials(k)];
          bigCont{chad} = [timeTrials(k+1);posTrials(k+1)];
      else
          bigCont{chad} = [timeTrials(k);posTrials(k)];
          smallCont{chad} = [timeTrials(k+1);posTrials(k+1)];
      end
   end
end

posOut = [];
tsqOut = [];

if and(~isempty(smallCont),~isempty(bigCont))
    for j = 1:length(smallCont)
        time = (smallCont{j}(1) + bigCont{j}(1))/2;
        [pos,tsq] = calcX5_1(limit,initState,finalState,time);
        
        if pos < goal
            smallCont{j} = [time;pos];
        else
            bigCont{j} = [time;pos];
        end
        
        for k = 1:200
            if (pos >= goal)
                time = (time + smallCont{j}(1)) / 2;
            else
                time = (time + bigCont{j}(1)) / 2;
            end
            [pos,tsq] = calcX5_1(limit,initState,finalState,time);
            
            if pos < goal
                smallCont{j} = [time;pos];
            else
                bigCont{j} = [time;pos];
            end
            if (norm(pos - goal) < 1e-4)
                posOut{j} = pos;
                tsqOut{j} = tsq;
                break;
            end
        end
    end
end

end
