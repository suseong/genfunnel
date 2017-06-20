function [pos,tsq] = calcX5_4(input,initState,finalState,tf)

u = input(1);
am = input(2);

x0 = initState(1);
v0 = initState(2);
a0 = initState(3);

xf = finalState(1);
vf = finalState(2);
af = finalState(3);

t4 = -(af + am - tf*u)/u;
t2 = -(- a0^2 + 2*a0*am + af^2 + 2*af*am + 4*am^2 - 2*tf*u*am + 2*u*v0 - 2*u*vf)/(4*am*u);
t3 = t2 + 2*am/u;
t1 = (am-a0)/u;

x5 = x0-1.0./u.^2.*(a0-am).^3.*(1.0./6.0)-(t4-tf).*(v0+((a0-am).^2.*(1.0./2.0))./u+am.*(t2+(a0-am)./u)+am.*(t2-t4+(am.*2.0)./u)-(a0.*(a0-am))./u)-am.*(t2-t4+(am.*2.0)./u).^2.*(1.0./2.0)+am.*(t2+(a0-am)./u).^2.*(1.0./2.0)-am.*(t4-tf).^2.*(1.0./2.0)-u.*(t4-tf).^3.*(1.0./6.0)-(t2-t4+(am.*2.0)./u).*(v0+((a0-am).^2.*(1.0./2.0))./u+am.*(t2+(a0-am)./u)-(a0.*(a0-am))./u)+(t2+(a0-am)./u).*(v0+((a0-am).^2.*(1.0./2.0))./u-(a0.*(a0-am))./u)+am.^3.*1.0./u.^2.*(2.0./3.0)-(v0.*(a0-am))./u+a0.*1.0./u.^2.*(a0-am).^2.*(1.0./2.0)+(am.*(v0+((a0-am).^2.*(1.0./2.0))./u+am.*(t2+(a0-am)./u)-(a0.*(a0-am))./u).*2.0)./u;

tsq = [t1 t2 t3 t4 tf];

pos = x5;

end