function [pos,acc,tsq] = calc_1_3(input,initState,finalState,tf)

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

    t3 = (3*tf^2*u^2 - 2*a0*af - 4*u*v0 + 4*u*vf + a0^2 + af^2 + 2*a0*tf*u - 6*af*tf*u)/(4*u*(a0 - af + tf*u));
    t1 = -(a0 - af - 2*t3*u + tf*u)/(2*u);
    
    pos = x0 - (tf^2*(a0 - af - 2*t3*u + tf*u))/2 + tf*v0 + (a0*tf^2)/2 - (a0 - af - 2*t3*u + tf*u)^3/(24*u^2) + (tf^3*u)/6 - (tf*(a0 - af - 2*t3*u + tf*u)^2)/(4*u) - (tf^2*(3*tf^2*u^2 - 2*a0*af - 4*u*v0 + 4*u*vf + a0^2 + af^2 + 2*a0*tf*u - 6*af*tf*u))/(4*(a0 - af + tf*u)) - (3*tf^2*u^2 - 2*a0*af - 4*u*v0 + 4*u*vf + a0^2 + af^2 + 2*a0*tf*u - 6*af*tf*u)^3/(192*u^2*(a0 - af + tf*u)^3) + (tf*(3*tf^2*u^2 - 2*a0*af - 4*u*v0 + 4*u*vf + a0^2 + af^2 + 2*a0*tf*u - 6*af*tf*u)^2)/(16*u*(a0 - af + tf*u)^2);
    vel = v0 + a0*tf - (a0 - af - 2*t3*u + tf*u)^2/(4*u) + (tf^2*u)/2 - tf*(a0 - af - 2*t3*u + tf*u) + (3*tf^2*u^2 - 2*a0*af - 4*u*v0 + 4*u*vf + a0^2 + af^2 + 2*a0*tf*u - 6*af*tf*u)^2/(16*u*(a0 - af + tf*u)^2) - (tf*(3*tf^2*u^2 - 2*a0*af - 4*u*v0 + 4*u*vf + a0^2 + af^2 + 2*a0*tf*u - 6*af*tf*u))/(2*(a0 - af + tf*u));
    tsq = [t1 t1 t3 t3 tf];
    
    a1 = a0 + t1*u;
    a2 = a1;
    a3 = a2 - (t3-t1)*u;
    a4 = a3;
    a5 = a4 + (tf-t3)*u;
    
    acc = [a1 a2 a3 a4 a5];
else
    pos = -1e10;
    tsq = [-100 -100 -100 -100 -100];
    acc = [-100 -100 -100 -100 -100];
end

% test
% [a,b,c] = calc_1_3(-[20 5],-[0 0 -10],-[-2.627 0 10],1.4)
% [a,b,c] = calc_1_3([20 5],[0 0 -3],[-0.0625 -0.25 -3],0.5)

end