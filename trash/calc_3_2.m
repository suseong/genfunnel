function [pos,acc,tsq] = calc_3_2(input,initState,finalState,tf)

cond_pos = and(and(input(1) > 0, initState(3) >= input(2)), finalState(3) <=-input(2));
cond_neg = and(and(input(1) < 0, initState(3) <= input(2)), finalState(3) >=-input(2));

if or(cond_pos, cond_neg)
    u = input(1);
    am = input(2);
    
    x0 = initState(1);
    v0 = initState(2);
    a0 = initState(3);
    
    xf = finalState(1);
    vf = finalState(2);
    af = finalState(3);

    t1 = 0;
    t2 = 0;
    t3 = (a0 + am)/u;
    t4 = (af + am + tf*u)/u;
    
    acc = [a0 a0 -am -am af];
    tsq = [t1 t2 t3 t4 tf];
    pos = x0 + tf*v0 - (tf^2*(a0 + am))/2 + (a0*tf^2)/2 + (tf^2*(af + am + tf*u))/2 - (tf^3*u)/6 - (a0 + am)^3/(6*u^2) + (af + am + tf*u)^3/(6*u^2) + (tf*(a0 + am)^2)/(2*u) - (tf*(af + am + tf*u)^2)/(2*u);
    
else
    pos = -1e10;
    tsq = [-100 -100 -100 -100 -100];
    acc = [-100 -100 -100 -100 -100];
end

% test
% [a,b,c] = calc_3_2(-[20 5],-[0 0 10],-[1.667 -1.25 -10],1.25)

end