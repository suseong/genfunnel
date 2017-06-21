function [p,v] = show_pv(tsq,acc,initState,actLimit,act,ptNum)

t = linspace(0,tsq(end),ptNum);
p = [];
v = [];

a(1) = initState(3);
for k=2:7
    a(k) = acc(k-1);
end

u = actLimit(1)*act;

t1 = tsq(1);
t2 = tsq(2);
t3 = tsq(3);
t4 = tsq(4);
tf = tsq(5);

a0 = acc(1); 
a1 = acc(2); 
a2 = acc(3); 
a3 = acc(4); 
a4 = acc(5); 
a5 = acc(6);

v0 = initState(2);
v1 = v0 + a0*t1      + 1/2*u*t1^2;
v2 = v1 + a1*(t2-t1);
v3 = v2 + a2*(t3-t2) - 1/2*u*(t3-t2)^2;
v4 = v3 + a3*(t4-t3);
v5 = v4 + a4*(tf-t4) + 1/2*u*(tf-t4)^2;

x0 = initState(1);
x1 = x0 + v0*t1      + 1/2*a0*t1^2      + 1/6*u*t1^3;
x2 = x1 + v1*(t2-t1) + 1/2*a1*(t2-t1)^2;
x3 = x2 + v2*(t3-t2) + 1/2*a2*(t3-t2)^2 - 1/6*u*(t3-t2)^3;
x4 = x3 + v3*(t4-t3) + 1/2*a3*(t4-t3)^2;
x5 = x4 + v4*(tf-t4) + 1/2*a4*(tf-t4)^2 + 1/6*u*(tf-t4)^3;

for k=1:length(t)
    
    if t(k) < tsq(1)
        p(k) = x0 + v0*t(k)      + 1/2*a0*t(k)^2      + 1/6*u*t(k)^3;
        v(k) = v0 + a0*t(k)      + 1/2*u*t(k)^2;
    elseif and(t(k) >= tsq(1),t(k) < tsq(2))
        p(k) = x1 + v1*(t(k)-t1) + 1/2*a1*(t(k)-t1)^2;
        v(k) = v1 + a1*(t(k)-t1);
    elseif and(t(k) >= tsq(2),t(k) < tsq(3))
        p(k) = x2 + v2*(t(k)-t2) + 1/2*a2*(t(k)-t2)^2 - 1/6*u*(t(k)-t2)^3;
        v(k) = v2 + a2*(t(k)-t2) - 1/2*u*(t(k)-t2)^2;
    elseif and(t(k) >= tsq(3),t(k) < tsq(4))
        p(k) = x3 + v3*(t(k)-t3) + 1/2*a3*(t(k)-t3)^2;
        v(k) = v3 + a3*(t(k)-t3);
    elseif and(t(k) >= tsq(4),t(k) <= tsq(5))
        p(k) = x4 + v4*(t(k)-t4) + 1/2*a4*(t(k)-t4)^2 + 1/6*u*(t(k)-t4)^3;
        v(k) = v4 + a4*(t(k)-t4) + 1/2*u*(t(k)-t4)^2;
    end
   
end

end