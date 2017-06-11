function [cost,T,c_x,c_y] = calc_min_cost4(T,p,opts)

k = 1/3;
dt = 0.03;
del_t(1,:) = dt*[0  1  1-k 1-2*k 0];
del_t(2,:) = dt*[0 -k  1-k 1-2*k 0];
del_t(3,:) = dt*[0 -k -2*k 1-2*k 0];
del_t(4,:) = dt*[0 -k -2*k -3*k  0];

[cost,c_x,c_y] = calc_cost4(T,p,opts);

cnt = 0;
sw = 0;
convergence = zeros(1,4);
for i=1:200
    cnt = cnt+1;
    if cnt > 4
       cnt = 1;
       if (sum(convergence) > 0)
           convergence = zeros(1,4);
       else
           sw = 1;
       end
    end
    
    T_ = T + del_t(cnt,:);
    
    check = T_(2:5) - T_(1:4);
    
    for j = 1:4
       if check(j) < 0
          disp('time constraint');
          sw = 1;
       end
    end
    
    if sw == 1
        break;
    end
    [cost_,c_x_,c_y_] = calc_cost4(T_,p,opts);
    if(cost_ < cost)
       convergence(cnt) = 1;
       T = T_;
       cost = cost_;
       cnt = cnt - 1;
       c_x = c_x_;
       c_y = c_y_;
    end
end
i
end