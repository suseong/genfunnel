function generateFunnelCSV(data,filename)

if nargin < 3
    if size(traj,2) == 8
        type = 'state';
        filename = [type,'_',datestr(now,'yyyymmdd'),'_',datestr(now,'HHMMSS'),'.csv'];
    elseif size(traj,2) == 2
        type = 'input';
        filename = [type,'_',datestr(now,'yyyymmdd'),'_',datestr(now,'HHMMSS'),'.csv'];
    elseif size(traj,2) == 10
        type = 'gain';
        filename = [type,'_',datestr(now,'yyyymmdd'),'_',datestr(now,'HHMMSS'),'.csv'];
    else
        error('trajectory not identified');
    end
else
    if size(traj,2) == 8
        type = 'state';
    elseif size(traj,2) == 2
        type = 'input';
    elseif size(traj,2) == 10
        type = 'gain';
    else
        error('trajectory not identified');
    end    
    typecheck(filename,'char');
end

N = length(ts);
traj_matrix = [];
[sR,sC] = size(traj);

if strcmp(type,'state')
       num = 9; order = 7;
       if and(sR == (N-1)*num, sC == order+1)
          for k = 1:N-1
              traj_matrix(k,:) = reshape(traj(num*(k-1)+1:num*k,:)',1,num*(order+1));
          end
       else
          error('dimension is not correct');
       end    
elseif strcmp(type,'input')
       [sR,sC] = size(traj);
       num = 4; order = 1;
       if and(sR == (N-1)*num, sC == order+1)
          for k = 1:N-1
              traj_matrix(k,:) = reshape(traj(num*(k-1)+1:num*k,:)',1,num*(order+1));
          end
       else
          error('dimension is not correct');
       end    
elseif strcmp(type,'gain')
       num = 36; order = 9;
       if and(sR == (N-1)*num, sC == order+1)
          for k = 1:N-1
              traj_matrix(k,:) = reshape(traj(num*(k-1)+1:num*k,:)',1,num*(order+1));
          end
       else
          error('dimension is not correct');
       end
end

traj_matrix = [ts' [traj_matrix;zeros(1,num*(order+1))]];
csvwrite(filename,traj_matrix);

end