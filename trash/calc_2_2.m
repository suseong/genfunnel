function [pos,acc,tsq] = calc_2_2(input,initState,finalState,tf)

cond_pos = and(and(input(1) > 0, initState(3) >= input(2)), finalState(3) <=-input(2));
cond_neg = and(and(input(1) < 0, initState(3) <= input(2)), finalState(3) >=-input(2));

if or(cond_pos, cond_neg)
    u = input(1);
    am = input(2);
    
    x0 = initState(1);
    v0 = initState(2);
    a0 = initState(3);
    
%     xf = finalState(1);
%     vf = finalState(2);
    af = finalState(3);

    t1 = (a0 - am)/u;
    t2 = (af - am + tf*u)/u;
    
    x5 = (- a0^3 + 3*a0^2*am + 3*a0^2*tf*u - 3*a0*am^2 - 6*a0*am*tf*u + am^3 + 3*am^2*tf*u + 3*am*tf^2*u^2 + t2^3*u^3 - 3*t2^2*tf*u^3 + 3*t2*tf^2*u^3 - tf^3*u^3 + 6*v0*tf*u^2 + 6*x0*u^2)/(6*u^2);

    a1 = a0 - t1*u;
    a2 = a1;
    a3 = a2;
    a4 = a3;
    a5 = a4 - (tf-t2)*u;
    
    tsq = [t1 t2 t2 t2 tf];
    acc = [a1 a2 a3 a4 a5];
    pos = x5;
else
    pos = -1e10;
    tsq = [-100 -100 -100 -100 -100];
    acc = [-100 -100 -100 -100 -100];
end

% test
% [a,b,c] = calc_2_2(-[20 5],-[0 0 10],-[3 5 -15],2)

end