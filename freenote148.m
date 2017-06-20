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
% a5_ = subs(a5,[t2 t4],[t1 t3]);
% v5_ = subs(v5,[t2 t4],[t1 t3]);
% x5_ = subs(x5,[t2 t4],[t1 t3]);
% 
% t1_ = (af - a0 + 2*u*t3 - u*tf)/(2*u);
% t3_ = (4*u*vf + (a0^2 - 2*a0*af + 2*a0*tf*u + af^2 - 6*af*tf*u + 3*tf^2*u^2 - 4*v0*u)) / (4*(a0-af)*u + 4*tf*u^2);
% 
% t3__ = subs(t3_,[u am x0 v0 a0 xf vf af tf],      [10 5 0 0 0 5 2 0 tfCand]);
% t1__ = subs(t1_,[u am x0 v0 a0 xf vf af tf t3],   [10 5 0 0 0 5 2 0 tfCand double(t3__)]);
% x5__ = subs(x5_,[u am x0 v0 a0 xf vf af tf t3 t1],[10 5 0 0 0 5 2 0 tfCand double(t3__) double(t1__)]);

%% case 2
% a5_ = subs(a5,[t1 t4],[(am-a0)/u t3]);
% v5_ = subs(v5,[t1 t4],[(am-a0)/u t3]);
% x5_ = subs(x5,[t1 t4],[(am-a0)/u t3]);
% 
% chad = simplify(af - a5_);
% A = simplify(diff(chad,t2));
% B = simplify(chad - A*t2);
% t2_ = -B / A;
% 
% v5_ = subs(v5_,t2,t2_);
% 
% chad = simplify(vf - v5_);
% A = simplify(diff(diff(chad,t3),t3)/2);
% B = simplify(diff(chad,t3) - 2*A*t3);
% C = simplify(chad - A*t3^2 - B*t3);

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
% a5_ = subs(a5,[t2 t4],[t1,tf-(af+am)/u]);
% v5_ = subs(v5,[t2 t4],[t1,tf-(af+am)/u]);
% x5_ = subs(x5,[t2 t4],[t1,tf-(af+am)/u]);
% 
% chad = simplify(af-a5_);
% A = simplify(diff(chad,t1));
% B = simplify(chad - A*t1);
% t1_ = -B / A;
% 
% v5_ = subs(v5_,t1,t1_);
% 
% chad = simplify(vf - v5_);
% A = simplify(diff(diff(chad,t3),t3)/2);
% B = simplify(diff(chad,t3) - 2*A*t3);
% C = simplify(chad - A*t3^2 - B*t3);

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

chad = simplify(vf - v5_);
A = simplify(diff(chad,t2));
B = simplify(chad - A*t2);
t2_ = -B / A;



A = simplify(diff(diff(chad,t2),t2)/2);
B = simplify(diff(chad,t2) - 2*A*t2);
C = simplify(chad - A*t2^2 - B*t2);








































