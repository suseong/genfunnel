function [pos,acc,tsq] = calc_3_4(input,initState,finalState,tf)

cond_pos = and(and(input(1)>0, initState(3)<=-input(2)), finalState(3)<=-input(2));
cond_neg = and(and(input(1)<0, initState(3)>=-input(2)), finalState(3)>=-input(2));

if or(cond_pos, cond_neg)
    u = input(1);
    am = input(2);
    
    x0 = initState(1);
    v0 = initState(2);
    a0 = initState(3);
    
    xf = finalState(1);
    vf = finalState(2);
    af = finalState(3);

    t4 = tf + (af + am)/u;
    
    A = u/4;
    B = a0/2 + am/2;
    C = - vf - (a0^2 + 2*a0*am + 2*af^2 + 4*af*am + 3*am^2 + 4*tf*u*am - 4*u*v0)/(4*u);
    
    t3_1 = (-B + sqrt(B^2 -4*A*C)) / (2*A);
    t1_1 = -(a0 + am - t3_1*u)/(2*u);
    t1 = t1_1; t2 = t1_1; t3 = t3_1;
    x5_1 = (3*af*am^2 - t3^3*u^3 + 3*af^2*am + 6*u^2*x0 + af^3 + am^3 - (a0 + am - t3*u)^3/4 - (3*tf*u*(a0 + am - t3*u)^2)/2 + 6*tf*u^2*v0 + 3*a0*tf^2*u^2 - 3*tf^2*u^2*(a0 + am - t3*u) - 3*t3*tf^2*u^3 + 3*t3^2*tf*u^3)/(6*u^2);
    a1 = a0 + t1*u;
    a2 = a1;
    a3 = a2 - (t3-t2)*u;
    a4 = a3;
    a5 = a4 - (tf-t4)*u;
    acc1 = [a1 a2 a3 a4 a5];

    t3_2 = (-B - sqrt(B^2 -4*A*C)) / (2*A);
    t1_2 = -(a0 + am - t3_2*u)/(2*u);
    t1 = t1_2; t2 = t1_2; t3 = t3_2;
    x5_2 = (3*af*am^2 - t3^3*u^3 + 3*af^2*am + 6*u^2*x0 + af^3 + am^3 - (a0 + am - t3*u)^3/4 - (3*tf*u*(a0 + am - t3*u)^2)/2 + 6*tf*u^2*v0 + 3*a0*tf^2*u^2 - 3*tf^2*u^2*(a0 + am - t3*u) - 3*t3*tf^2*u^3 + 3*t3^2*tf*u^3)/(6*u^2);
    a1 = a0 + t1*u;
    a2 = a1;
    a3 = a2 - (t3-t2)*u;
    a4 = a3;
    a5 = a4 - (tf-t4)*u;
    acc2 = [a1 a2 a3 a4 a5];
    
    tsq1 = [t1_1 t1_1 t3_1 t4 tf];
    tsq2 = [t1_2 t1_2 t3_2 t4 tf];
    
    pos = real([x5_1;x5_2]);
    tsq = real([tsq1;tsq2]);
    acc = real([acc1;acc2]);
else
    pos = -1e10;
    tsq = [-100 -100 -100 -100 -100];
    acc = [-100 -100 -100 -100 -100];
end

% test
% [a,b,c] = calc_3_4( [20 5], [0 0 -10],[-7.188 -8.75 -10],1.75)
% [a,b,c] = calc_3_4( -[20 5], -[0 0 -10],-[-7.188 -8.75 -10],1.75)

end