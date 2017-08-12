function [pos,acc,tsq] = calc_4_1(input,initState,finalState,tf)

cond_pos = and(and(input(1)>0, initState(3)>= input(2)), finalState(3)>= input(2));
cond_neg = and(and(input(1)<0, initState(3)<= input(2)), finalState(3)<= input(2));

if or(cond_pos, cond_neg)
    u = input(1);
    am = input(2);
    
    x0 = initState(1);
    v0 = initState(2);
    a0 = initState(3);
    
    xf = finalState(1);
    vf = finalState(2);
    af = finalState(3);

    t1 = (a0 - am)/u;
    t4 = tf - (af + am)/u;
    t3 = (- a0^2 + 2*a0*am - af^2 - 2*af*am + 2*am^2 + 2*tf*u*am - 2*u*v0 + 2*u*vf)/(4*am*u);
    t2 = -(2*am - t3*u)/u;
    
    a1 = a0 - t1*u;
    a2 = a1;
    a3 = a2 - (t3-t2)*u;
    a4 = a3;
    a5 = a4 + (tf-t4)*u;
    
    acc = [a1 a2 a3 a4 a5];
    tsq = [t1 t2 t3 t4 tf];
    pos = -(3*a0^4 - 4*a0^3*am + 6*a0^2*af^2 + 12*a0^2*af*am - 12*a0^2*am*tf*u + 12*a0^2*u*v0 - 12*a0^2*u*vf - 12*a0*af^2*am - 24*a0*af*am^2 + 24*a0*am^2*tf*u - 24*a0*am*u*v0 + 24*a0*am*u*vf + 3*af^4 + 4*af^3*am + 12*af^2*am*tf*u + 12*af^2*u*v0 - 12*af^2*u*vf + 24*af*am^2*tf*u + 24*af*am*u*v0 - 24*af*am*u*vf + 12*am^4 - 12*am^2*tf^2*u^2 + 24*am^2*u*v0 - 24*am^2*u*vf - 24*am*tf*u^2*v0 - 24*am*tf*u^2*vf - 48*x0*am*u^2 + 12*u^2*v0^2 - 24*u^2*v0*vf + 12*u^2*vf^2)/(48*am*u^2);
    
else
    pos = -1e10;
    tsq = [-100 -100 -100 -100 -100];
    acc = [-100 -100 -100 -100 -100];
end

% test
% [a,b,c] = calc_4_1([20 5],[0 0 10],[7.188 3.75 10],2.5)
% [a,b,c] = calc_4_1(-[20 5],-[0 0 10],-[7.188 3.75 10],2.5)

end