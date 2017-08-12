clear all
close all
clc

syms t1 t2 t3 t4 tf real
syms u am x0 v0 a0 xf vf af real

a1 = a0 + t1*u;
a2 = a1;
a3 = a2 - (t3-t2)*u;
a4 = a3;
a5 = a4 + (tf-t4)*u;

v1 = v0 + a0*t1      + 1/2*u*t1^2;
v2 = v1 + a1*(t2-t1);
v3 = v2 + a2*(t3-t2) - 1/2*u*(t3-t2)^2;
v4 = v3 + a3*(t4-t3);
v5 = v4 + a4*(tf-t4) + 1/2*u*(tf-t4)^2;

x1 = x0 + v0*t1      + 1/2*a0*t1^2      + 1/6*u*t1^3;
x2 = x1 + v1*(t2-t1) + 1/2*a1*(t2-t1)^2;
x3 = x2 + v2*(t3-t2) + 1/2*a2*(t3-t2)^2 - 1/6*u*(t3-t2)^3;
x4 = x3 + v3*(t4-t3) + 1/2*a3*(t4-t3)^2;
x5 = x4 + v4*(tf-t4) + 1/2*a4*(tf-t4)^2 + 1/6*u*(tf-t4)^3;

%% case 1
a5_ = subs(a5,[t2 t4],[t1 t3]);
v5_ = subs(v5,[t2 t4],[t1 t3]);
x5_ = subs(x5,[t2 t4],[t1 t3]);

chad = simplify(af - a5_);
A = simplify(diff(chad,t1));
B = simplify(chad - A*t1);
t1_ = -B / A;

v5_ = simplify(subs(v5_,t1,t1_));
x5_ = simplify(subs(x5_,t1,t1_));

chad = simplify(vf - v5_);
A = simplify(diff(chad,t3));
B = simplify(chad - A*t3);
t3_ = -B / A;

x5_ = simplify(subs(x5_,t3,t3_));


-(a0^4 - 4*a0^3*af + 4*a0^3*tf*u + 6*a0^2*af^2 - 12*a0^2*af*tf*u + 6*a0^2*tf^2*u^2 - 4*a0*af^3 + 12*a0*af^2*tf*u + 36*a0*af*tf^2*u^2 - 12*a0*tf^3*u^3 
- 96*a0*tf*u^2*vf - 96*x0*a0*u^2 + af^4 - 4*af^3*tf*u + 6*af^2*tf^2*u^2 + 12*af*tf^3*u^3 + 96*af*tf*u^2*v0 + 96*x0*af*u^2 
- 3*tf^4*u^4 - 48*tf^2*u^3*v0 - 48*tf^2*u^3*vf - 96*x0*tf*u^3 + 48*u^2*v0^2 - 96*u^2*v0*vf + 48*u^2*vf^2)
/(96*u^2*(a0 - af + tf*u))


% t1_ = (af - a0 + 2*u*t3 - u*tf)/(2*u);
% t3_ = (4*u*vf + (a0^2 - 2*a0*af + 2*a0*tf*u + af^2 - 6*af*tf*u + 3*tf^2*u^2 - 4*v0*u)) / (4*(a0-af)*u + 4*tf*u^2);
% 
% t3__ = subs(t3_,[u am x0 v0 a0 xf vf af tf],      [10 5 0 0 0 5 2 0 tfCand]);
% t1__ = subs(t1_,[u am x0 v0 a0 xf vf af tf t3],   [10 5 0 0 0 5 2 0 tfCand double(t3__)]);
% x5__ = subs(x5_,[u am x0 v0 a0 xf vf af tf t3 t1],[10 5 0 0 0 5 2 0 tfCand double(t3__) double(t1__)]);

%% case 2
clear all
close all
clc

syms t1 t2 t3 t4 tf real
syms u am x0 v0 a0 xf vf af real
syms dt1 dt2 dt3 dt4 dt5 real

a1 = a0 + dt1*u;
a2 = a1;
a3 = a2 - dt3*u;
a4 = a3;
a5 = a4 + dt5*u;

v1 = v0 + a0*dt1 + 1/2*u*dt1^2;
v2 = v1 + a1*dt2;
v3 = v2 + a2*dt3 - 1/2*u*dt3^2;
v4 = v3 + a3*dt4;
v5 = v4 + a4*dt5 + 1/2*u*dt5^2;

x1 = x0 + v0*dt1 + 1/2*a0*dt1^2 + 1/6*u*dt1^3;
x2 = x1 + v1*dt2 + 1/2*a1*dt2^2;
x3 = x2 + v2*dt3 + 1/2*a2*dt3^2 - 1/6*u*dt3^3;
x4 = x3 + v3*dt4 + 1/2*a3*dt4^2;
x5 = x4 + v4*dt5 + 1/2*a4*dt5^2 + 1/6*u*dt5^3;

a5_ = simplify(subs(a5,[dt1 dt4],[(am-a0)/u 0]));
v5_ = simplify(subs(v5,[dt1 dt4],[(am-a0)/u 0]));
x5_ = simplify(subs(x5,[dt1 dt4],[(am-a0)/u 0]));

chad = simplify(af - a5_);
A = simplify(diff(chad,dt5));
B = simplify(chad - A*dt5);
dt5_ = -B / A;

v5_ = simplify(subs(v5_,dt5,dt5_));
x5_ = simplify(subs(x5_,dt5,dt5_));

chad = simplify(vf - v5_);
A = simplify(diff(chad,dt2));
B = simplify(chad - A*dt2);
dt2_ = -B / A;

x5_ = simplify(subs(x5_,dt2,dt2_));
(- 3*a0^4 + 8*a0^3*am - 6*a0^2*am^2 + 12*a0^2*u*v0 - 24*a0*am*u*v0 + 3*af^4 - 8*af^3*am + 6*af^2*am^2 
- 12*af^2*dt3^2*u^2 - 12*af^2*u*vf + 24*af*am*u*vf + 12*am^2*dt3^2*u^2 + 12*am^2*u*v0 - 12*am^2*u*vf 
- 24*am*dt3^3*u^3 + 24*x0*am*u^2 + 12*dt3^4*u^4 + 24*dt3^2*u^3*vf - 12*u^2*v0^2 + 12*u^2*vf^2)
/(24*am*u^2)


a5_ = simplify(subs(a5,[t1 t4],[(am-a0)/u t3]));
v5_ = simplify(subs(v5,[t1 t4],[(am-a0)/u t3]));
x5_ = simplify(subs(x5,[t1 t4],[(am-a0)/u t3]));

chad = simplify(af - a5_);
A = simplify(diff(chad,tf));
B = simplify(chad - A*tf);
tf_ = -B / A;

v5_ = simplify(subs(v5_,tf,tf_));
x5_ = simplify(subs(x5_,tf,tf_));

chad = simplify(vf - v5_);
A = simplify(diff(diff(chad,t3),t3)/2);
B = simplify(diff(chad,t3) - 2*A*t3);
C = simplify(chad - A*t3^2 - B*t3);

t31_ = (-B + sqrt(B^2-4*A*C))/(2*A);
t32_ = (-B - sqrt(B^2-4*A*C))/(2*A);

chad1 = simplify(subs(x5_,t2,t2_));
chad1 = simplify(subs(chad1,t3,t31_));



% % tfCand = 2.14;
% % A_ = double(subs(A,[u am x0 v0 a0 xf vf af tf],[10 5 0 0 0 5 2 0 tfCand]));
% % B_ = double(subs(B,[u am x0 v0 a0 xf vf af tf],[10 5 0 0 0 5 2 0 tfCand]));
% % C_ = double(subs(C,[u am x0 v0 a0 xf vf af tf],[10 5 0 0 0 5 2 0 tfCand]));
% % 
% % t1__ = double(subs((am-a0)/u,[u am a0],[10 5 0]));
% % % t3__ = (-B_ + sqrt(B_^2-4*A_*C_))/(2*A_)
% % % t2__ = double(subs(t2_,[u am x0 v0 a0 xf vf af tf t3],[10 5 0 0 0 5 2 0 tfCand t3__]))
% % 
% % t3__ = (-B_ - sqrt(B_^2-4*A_*C_))/(2*A_);
% % t2__ = double(subs(t2_,[u am x0 v0 a0 xf vf af tf t3],[10 5 0 0 0 5 2 0 tfCand t3__]));
% % 
% % x5__ = double(subs(x5_,[u am x0 v0 a0 xf vf af tf t3 t2 t1],[10 5 0 0 0 5 2 0 tfCand double(t3__) double(t2__) double(t1__)]))
% % 
% % a1_ = double(subs(a1,[u am x0 v0 a0 xf vf af tf t4 t3 t2 t1],[10 5 0 0 0 5 2 0 tfCand double(t3__) double(t3__) double(t2__) double(t1__)]));
% % a2_ = double(subs(a2,[u am x0 v0 a0 xf vf af tf t4 t3 t2 t1],[10 5 0 0 0 5 2 0 tfCand double(t3__) double(t3__) double(t2__) double(t1__)]));
% % a3_ = double(subs(a3,[u am x0 v0 a0 xf vf af tf t4 t3 t2 t1],[10 5 0 0 0 5 2 0 tfCand double(t3__) double(t3__) double(t2__) double(t1__)]));
% % a4_ = double(subs(a4,[u am x0 v0 a0 xf vf af tf t4 t3 t2 t1],[10 5 0 0 0 5 2 0 tfCand double(t3__) double(t3__) double(t2__) double(t1__)]));
% % a5_ = double(subs(a5,[u am x0 v0 a0 xf vf af tf t4 t3 t2 t1],[10 5 0 0 0 5 2 0 tfCand double(t3__) double(t3__) double(t2__) double(t1__)]));

%% u > 0, am > 0 or u < 0, am < 0
a5_ = subs(a5,[t2 t4],[t1,tf-(af+am)/u]);
v5_ = subs(v5,[t2 t4],[t1,tf-(af+am)/u]);
x5_ = subs(x5,[t2 t4],[t1,tf-(af+am)/u]);

chad = simplify(af-a5_);
A = simplify(diff(chad,t1));
B = simplify(chad - A*t1);
t1_ = -B / A;

v5_ = subs(v5_,t1,t1_);
x5_ = simplify(subs(x5_,t1,t1_));

chad = simplify(vf - v5_);
A = simplify(diff(chad,tf));
B = simplify(chad - A*tf);
tf_ = -B / A;

x5_ = simplify(subs(x5_,tf,tf_));

% (3*a0^4 + 8*a0^3*am - 12*a0^3*t3*u + 6*a0^2*am^2 - 24*a0^2*am*t3*u + 6*a0^2*t3^2*u^2 - 24*a0^2*u*v0 
% - 12*a0*am^2*t3*u - 48*a0*am*u*v0 + 12*a0*t3^3*u^3 + 48*a0*t3*u^2*v0 - 12*af^4 - 32*af^3*am - 24*af^2*am^2 
% + 48*af^2*u*vf + 96*af*am*u*vf + 3*am^4 - 6*am^2*t3^2*u^2 - 24*am^2*u*v0 + 48*am^2*u*vf + 48*am*t3*u^2*v0 
% + 96*x0*am*u^2 + 3*t3^4*u^4 + 24*t3^2*u^3*v0 + 48*u^2*v0^2 - 48*u^2*vf^2)
% /(96*am*u^2)

chad = simplify(vf - v5_);
A = simplify(diff(diff(chad,t3),t3)/2);
B = simplify(diff(chad,t3) - 2*A*t3);
C = simplify(chad - A*t3^2 - B*t3);

% % tfCand = 2.1;
% % A_ = double(subs(A,[u am x0 v0 a0 xf vf af tf],[-10 5 0 0 0 5 2 0 tfCand]));
% % B_ = double(subs(B,[u am x0 v0 a0 xf vf af tf],[-10 5 0 0 0 5 2 0 tfCand]));
% % C_ = double(subs(C,[u am x0 v0 a0 xf vf af tf],[-10 5 0 0 0 5 2 0 tfCand]));
% % 
% % t3__ = (-B_ - sqrt(B_^2-4*A_*C_))/(2*A_)
% % t3__ = (-B_ + sqrt(B_^2-4*A_*C_))/(2*A_)
% % 
% % 
% % t1__ = double(subs(t1_,[u am x0 v0 a0 xf vf af tf t4 t3 t2 t1],[-10 5 0 0 0 5 2 0 tfCand 0 double(t3__) 0 0]));
% % t2__ = t1__;
% % t4_ = tf - (af - am) / u;
% % t4__ = double(subs(t4_,[u am x0 v0 a0 xf vf af tf t4 t3 t2 t1],[-10 5 0 0 0 5 2 0 tfCand 0 double(t3__) double(t1__) double(t1__)]));
% % 
% % x5__ = double(subs(x5_,[u am x0 v0 a0 xf vf af tf t4 t3 t2 t1],[-10 5 0 0 0 5 2 0 tfCand 0 double(t3__) double(t1__) double(t1__)]));
% % 
% % a1_ = double(subs(a1,[u am x0 v0 a0 xf vf af tf t4 t3 t2 t1],[10 5 0 0 0 5 2 0 tfCand double(t4__) double(t3__) double(t1__) double(t1__)]));
% % a2_ = double(subs(a2,[u am x0 v0 a0 xf vf af tf t4 t3 t2 t1],[10 5 0 0 0 5 2 0 tfCand double(t4__) double(t3__) double(t1__) double(t1__)]));
% % a3_ = double(subs(a3,[u am x0 v0 a0 xf vf af tf t4 t3 t2 t1],[10 5 0 0 0 5 2 0 tfCand double(t4__) double(t3__) double(t1__) double(t1__)]));
% % a4_ = double(subs(a4,[u am x0 v0 a0 xf vf af tf t4 t3 t2 t1],[10 5 0 0 0 5 2 0 tfCand double(t4__) double(t3__) double(t1__) double(t1__)]));
% % a5_ = double(subs(a5,[u am x0 v0 a0 xf vf af tf t4 t3 t2 t1],[10 5 0 0 0 5 2 0 tfCand double(t4__) double(t3__) double(t1__) double(t1__)]));

%%
a5_ = subs(a5,[t1 t3],[(am-a0)/u, t2 + 2*am/u]);
v5_ = subs(v5,[t1 t3],[(am-a0)/u, t2 + 2*am/u]);
x5_ = subs(x5,[t1 t3],[(am-a0)/u, t2 + 2*am/u]);

chad = simplify(af-a5_);
A = simplify(diff(chad,t4));
B = simplify(chad - A*t4);
t4_ = -B / A;

v5_ = subs(v5_,t4,t4_);
x5_ = subs(x5_,t4,t4_);

chad = simplify(vf - v5_);
A = simplify(diff(chad,t2));
B = simplify(chad - A*t2);
t2_ = -B / A;

chad = simplify(subs(x5_,t2,t2_));

A = simplify(diff(diff(chad,t2),t2)/2);
B = simplify(diff(chad,t2) - 2*A*t2);
C = simplify(chad - A*t2^2 - B*t2);

chad = simplify(diff(chad,tf));
A = simplify(diff(chad,tf));
B = simplify(chad - A*tf);

tf_star = -B / A;

chad = simplify(subs(x5_ - xf,[t2 t4],[t2_ t4_]));
chad = simplify(subs(chad,tf,tf_star));

A = simplify(diff(chad,tf));
B = simplify(chad -xf - A*tf);

tf_star = -B / A;









































