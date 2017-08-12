function [pos,acc,tsq] = calc_1_2(input,initState,finalState,tf)

cond_pos = and(and(input(1) > 0, initState(3) >= input(2)), finalState(3) <=-input(2));
cond_neg = and(and(input(1) < 0, initState(3) <= input(2)), finalState(3) >=-input(2));

if and(or(cond_pos, cond_neg),tf == abs(finalState(3)-initState(3))/input(1))
    u = input(1);
%     am = input(2);
    
    x0 = initState(1);
    v0 = initState(2);
    a0 = initState(3);
    
%     xf = finalState(1);
%     vf = finalState(2);
    af = finalState(3);
    
%     t4 = (a0 - af + tf*u)/(2*u);
%     a4 = a0-t4*u;
    
    pos = - (u*tf^3)/6 + (a0*tf^2)/2 + v0*tf + x0;
%     vel = v0 + a0*tf - (tf^2*u)/2
    acc = [a0 a0 af af af];
    tsq = [0 0 tf tf tf];
    
else
    pos = -1e10;
    tsq = [-100 -100 -100 -100 -100];
    acc = [-100 -100 -100 -100 -100];
end

% test
% [a,b,c] = calc_1_2([20 5],[0 0 10],-[3 5 14],1.2)

end