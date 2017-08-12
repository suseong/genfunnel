function [pos,acc,tsq] = calc_2_3(input,initState,finalState,tf)

cond_pos = and(and(input(1)>0, initState(3)<= input(2)), finalState(3)>=-input(2));
cond_neg = and(and(input(1)<0, initState(3)>= input(2)), finalState(3)<=-input(2));

if or(cond_pos, cond_neg)
    u = input(1);
    am = input(2);
    
    x0 = initState(1);
    v0 = initState(2);
    a0 = initState(3);
    
    xf = finalState(1);
    vf = finalState(2);
    af = finalState(3);

    t1 = -(a0 - am)/u;
    
    A = -u;
    B = 2*am - 2*af + 2*tf*u;
    C = - vf - (a0^2 - 2*a0*am + af^2 - 2*af*am - 4*af*tf*u + 2*am^2 + 2*am*tf*u + 2*tf^2*u^2 - 2*v0*u)/(2*u);
    
    t3_1 = (-B + sqrt(B^2 -4*A*C)) / (2*A);
    t2_1 = -(a0 - af + t1*u - 2*t3_1*u + tf*u)/u;
    t2 = t2_1;
    t3 = t3_1;
    x5_1 = x0 - (a0 - am)^3/(6*u^2) + (tf^2*(af - am + 2*t3*u - tf*u))/2 + tf*v0 + (a0*tf^2)/2 + (af - am + 2*t3*u - tf*u)^3/(6*u^2) - (t3^3*u)/3 + (tf^3*u)/6 - (tf^2*(a0 - am))/2 - (tf*(af - am + 2*t3*u - tf*u)^2)/(2*u) - t3*tf^2*u + t3^2*tf*u - (tf*(a0 - am)^2)/(2*u);
    a1 = a0 + t1*u;
    a2 = a1;
    a3 = a2 - (t3-t2)*u;
    a4 = a3;
    a5 = a4 + (tf-t3)*u;    
    acc1 = [a1 a2 a3 a4 a5];
    
    t3_2 = (-B - sqrt(B^2 -4*A*C)) / (2*A);
    t2_2 = -(a0 - af + t1*u - 2*t3_2*u + tf*u)/u;
    t2 = t2_2;
    t3 = t3_2;
    t3_2 = (-B - sqrt(B^2 -4*A*C)) / (2*A);
    t2_2 = -(a0 - af + t1*u - 2*t3_2*u + tf*u)/u;
    x5_2 = x0 - (a0 - am)^3/(6*u^2) + (tf^2*(af - am + 2*t3*u - tf*u))/2 + tf*v0 + (a0*tf^2)/2 + (af - am + 2*t3*u - tf*u)^3/(6*u^2) - (t3^3*u)/3 + (tf^3*u)/6 - (tf^2*(a0 - am))/2 - (tf*(af - am + 2*t3*u - tf*u)^2)/(2*u) - t3*tf^2*u + t3^2*tf*u - (tf*(a0 - am)^2)/(2*u);
    a1 = a0 + t1*u;
    a2 = a1;
    a3 = a2 - (t3-t2)*u;
    a4 = a3;
    a5 = a4 + (tf-t3)*u;
    acc2 = [a1 a2 a3 a4 a5];

    tsq1 = [t1 t2_1 t3_1 t3_1 tf];
    tsq2 = [t1 t2_2 t3_2 t3_2 tf];
    
    pos = real([x5_1;x5_2]);
    tsq = real([tsq1;tsq2]);
    acc = real([acc1;acc2]);
else
    pos = -1e10;
    tsq = [-100 -100 -100 -100 -100];
    acc = [-100 -100 -100 -100 -100];
end

% test
% [a,b,c] = calc_2_3([20 5],[0 0 -10],[-1.354 2.5 10],1.75)
% [a,b,c] = calc_2_3(-[20 5],-[0 0 -10],-[-1.354 2.5 10],1.75)

end