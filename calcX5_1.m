function [x5_,tsq] = calcX5_1(input,initState,finalState,tf)

u = input(1);
% am = input(2);

x0 = initState(1);
v0 = initState(2);
a0 = initState(3);

% xf = finalState(1);
vf = finalState(2);
af = finalState(3);

t3 = (4*u*vf + (a0^2 - 2*a0*af + 2*a0*tf*u + af^2 - 6*af*tf*u + 3*tf^2*u^2 - 4*v0*u)) / (4*(a0-af)*u + 4*tf*u^2);
t1 = (af - a0 + 2*u*t3 - u*tf)/(2*u);

t4 = t3-tf;
t5 = t1-t3;
t6 = t1.^2;
t7 = a0.*t1;
t8 = t6.*u.*(1.0./2.0);
t9 = a0.*(1.0./2.0);
t10 = t1.*u.*(1.0./2.0);
t11 = t5.^2;
t12 = t4.^2;
x5_ = x0+t11.*(t9+t10)+a0.*t6.*(1.0./2.0)+t1.*v0+t12.*(t9+t10+t5.*u.*(1.0./2.0))-t5.*(t7+t8+v0)-t4.*(t7+t8+v0-t11.*u.*(1.0./2.0)-t5.*(a0+t1.*u))+t1.*t6.*u.*(1.0./6.0)-t4.*t12.*u.*(1.0./6.0)+t5.*t11.*u.*(1.0./6.0);

tsq = [t1 t1 t3 t3 tf];

end