clear all
% close all
% clc

% syms t1 t2 t3 t4 tf real
syms u am x0 v0 a0 xf vf af real
syms dt real;

%% acc case 2 :  2 -> 2
clear
syms u am x0 v0 a0 xf vf af real
syms dt real;

t_in = -(am + a0) / u;
t_out =-(am + af) / u;

a1 = -am; % a0 - u_in*t_in = a0 - u_in*(a0 + am)/u_in
a2 = -am;
a3 = af; % am + u_out*t_out = -am + u_out*(af + am)/u_out

v1 = v0 + a0*t_in + 1/2*u*t_in^2;
v2 = v1 + a1*dt;
v3 = simplify(v2 + a2*t_out - 1/2*u*t_out^2); % vf

x1 = x0 + v0*t_in + 1/2*a0*t_in^2 + 1/6*u*t_in^3;
x2 = x1 + v1*dt + 1/2*a1*dt^2;
x3 = simplify(x2 + v2*t_out + 1/2*a2*t_out^2 - 1/6*u*t_out^3); % xf

% solve(vf - v3 == 0, dt)
dt_ = -(vf - (- a0^2 - af^2 + 2*am^2 + 2*u*v0)/(2*u))/am;
% x3_ = simplify(subs(x3-xf,dt,dt_));
err = -(12*u^2*vf^2 - 12*u^2*v0^2 - 8*a0^3*am + 8*af^3*am - 3*a0^4 + 3*af^4 - 6*a0^2*am^2 + 6*af^2*am^2 + 12*a0^2*u*v0 + 12*am^2*u*v0 + 12*af^2*u*vf + 12*am^2*u*vf - 24*am*u^2*x0 + 24*am*u^2*xf + 24*a0*am*u*v0 + 24*af*am*u*vf)/(24*am*u^2)

% dt__ = []; vf__ = []; xf = [];
% am__ = 0.1:0.1:7.9;
% vf_ = -10.5;
% a0_ = -10;
% af_ = -10;
% 
% for kk = 1:length(am__)
%     dt__(kk) = subs(dt_,[u x0 v0 a0 am vf af],[20 3 2 a0_ am__(kk) vf_ af_]);
%     xf(kk) = double(subs(x3,[u x0 v0 a0 am vf af dt],[20 3 2 a0_ am__(kk) vf_ af_ dt__(kk)]));
% end
% 
% figure(110)
% plot(am__,xf)
% 
%% acc case 2 :  2 -> 1
clear
syms u am x0 v0 a0 xf vf af real
syms dt real;

t_in = -(am + a0) / u;
t_out = (am + af) / u;

a1 = -am; % a0 - u_in*t_in = a0 - u_in*(a0 + am)/u_in
a2 = -am;
a3 = af; % am + u_out*t_out = -am + u_out*(af + am)/u_out

v1 = v0 + a0*t_in + 1/2*u*t_in^2;
v2 = v1 + a1*dt;
v3 = simplify(v2 + a2*t_out + 1/2*u*t_out^2); % vf

x1 = x0 + v0*t_in + 1/2*a0*t_in^2 + 1/6*u*t_in^3;
x2 = x1 + v1*dt + 1/2*a1*dt^2;
x3 = simplify(x2 + v2*t_out + 1/2*a2*t_out^2 + 1/6*u*t_out^3); % xf

% solve(vf - v3 == 0, dt)
dt_ = -(vf - (- a0^2 + af^2 + 2*u*v0)/(2*u))/am;
% x3_ = simplify(subs(x3-xf,dt,dt_));
err = (12*u^2*v0^2 - 12*u^2*vf^2 + 8*a0^3*am - 8*af^3*am + 3*a0^4 - 3*af^4 + 6*a0^2*am^2 - 6*af^2*am^2 - 12*a0^2*u*v0 - 12*am^2*u*v0 + 12*af^2*u*vf + 12*am^2*u*vf + 24*am*u^2*x0 - 24*am*u^2*xf - 24*a0*am*u*v0 + 24*af*am*u*vf)/(24*am*u^2)

% dt__ = []; vf__ = []; xf = [];
% am__ = 0.1:0.1:7.9;
% v0_ = -0.2467;
% vf_ = -1.2436;
% a0_ = -2.446;
% af_ = 3.363;
% 
% for kk = 1:length(am__)
%     dt__(kk) = subs(dt_,[u x0 v0 a0 am vf af],[20 0 v0_ a0_ am__(kk) vf_ af_]);
%     xf(kk) = double(subs(x3,[u x0 v0 a0 am vf af dt],[20 0 v0_ a0_ am__(kk) vf_ af_ dt__(kk)]));
% end
% 
% figure(110)
% plot(am__,xf)
% 
%% acc case 2 :  1 -> 2
clear
syms u am x0 v0 a0 xf vf af real
syms dt real;

t_in = (am + a0) / u;
t_out = -(am + af) / u;

a1 = -am; % a0 - u_in*t_in = a0 - u_in*(a0 + am)/u_in
a2 = -am;
a3 = af; % am + u_out*t_out = -am + u_out*(af + am)/u_out

v1 = v0 + a0*t_in - 1/2*u*t_in^2;
v2 = v1 + a1*dt;
v3 = simplify(v2 + a2*t_out - 1/2*u*t_out^2); % vf

x1 = x0 + v0*t_in + 1/2*a0*t_in^2 - 1/6*u*t_in^3;
x2 = x1 + v1*dt + 1/2*a1*dt^2;
x3 = simplify(x2 + v2*t_out + 1/2*a2*t_out^2 - 1/6*u*t_out^3); % xf

% solve(vf - v3 == 0, dt)
dt_ = -(vf - (a0^2 - af^2 + 2*u*v0)/(2*u))/am;
% x3_ = simplify(subs(x3-xf,dt,dt_));
err = (12*u^2*v0^2 - 12*u^2*vf^2 + 8*a0^3*am - 8*af^3*am + 3*a0^4 - 3*af^4 + 6*a0^2*am^2 - 6*af^2*am^2 + 12*a0^2*u*v0 + 12*am^2*u*v0 - 12*af^2*u*vf - 12*am^2*u*vf + 24*am*u^2*x0 - 24*am*u^2*xf + 24*a0*am*u*v0 - 24*af*am*u*vf)/(24*am*u^2)

% % dt__ = []; vf__ = []; xf = [];
% % am__ = 0.1:0.1:7.9;
% % vf_ = -4.25;
% % a0_ = 10;
% % af_ = -10;
% % 
% % for kk = 1:length(am__)
% %     dt__(kk) = subs(dt_,[u x0 v0 a0 am vf af],[20 3 2 a0_ am__(kk) vf_ af_]);
% %     xf(kk) = double(subs(x3,[u x0 v0 a0 am vf af dt],[20 3 2 a0_ am__(kk) vf_ af_ dt__(kk)]));
% % end
% % 
% % figure(110)
% % plot(am__,xf)
% 
%% acc case 2 :  1 -> 1
clear
syms u am x0 v0 a0 xf vf af real
syms dt real;

t_in = (am + a0) / u;
t_out = (am + af) / u;

a1 = -am; % a0 - u_in*t_in = a0 - u_in*(a0 + am)/u_in
a2 = -am;
a3 = af; % am + u_out*t_out = -am + u_out*(af + am)/u_out

v1 = v0 + a0*t_in - 1/2*u*t_in^2;
v2 = v1 + a1*dt;
v3 = simplify(v2 + a2*t_out + 1/2*u*t_out^2); % vf

x1 = x0 + v0*t_in + 1/2*a0*t_in^2 - 1/6*u*t_in^3;
x2 = x1 + v1*dt + 1/2*a1*dt^2;
x3 = simplify(x2 + v2*t_out + 1/2*a2*t_out^2 + 1/6*u*t_out^3); % xf

solve(vf - v3 == 0, dt)
dt_ = -(vf - (a0^2 + af^2 - 2*am^2 + 2*u*v0)/(2*u))/am;
x3_ = simplify(subs(x3-xf,dt,dt_));
% err = (12*u^2*v0^2 - 12*u^2*vf^2 + 8*a0^3*am - 8*af^3*am + 3*a0^4 - 3*af^4 + 6*a0^2*am^2 - 6*af^2*am^2 + 12*a0^2*u*v0 + 12*am^2*u*v0 + 12*af^2*u*vf + 12*am^2*u*vf + 24*am*u^2*x0 - 24*am*u^2*xf + 24*a0*am*u*v0 + 24*af*am*u*vf)/(24*am*u^2)


% dt__ = []; vf__ = []; xf = [];
% am__ = 0.1:0.1:7.9;
% v0_ = 1.4684;
% vf_ = 1.1015;
% a0_ = 1.49;
% af_ = 1.702;
% 
% for kk = 1:length(am__)
%     dt__(kk) = subs(dt_,[u x0 v0 a0 am vf af],[20 0 v0_ a0_ am__(kk) vf_ af_]);
%     xf(kk) = double(subs(x3,[u x0 v0 a0 am vf af dt],[20 0 v0_ a0_ am__(kk) vf_ af_ dt__(kk)]));
% end
% 
% figure(110)
% plot(am__,xf)

%% acc case 1 :  2 -> 2
syms u am x0 v0 a0 xf vf af real
syms dt real;

t_in = (am - a0) / u;
t_out = (am - af) / u;

a1 = am; % a0 - u_in*t_in = a0 - u_in*(a0 - am)/u_in
a2 = am;
a3 = af; % am + u_out*t_out = am + u_out*(af - am)/u_out

v1 = v0 + a0*t_in + 1/2*u*t_in^2;
v2 = v1 + a1*dt;
v3 = simplify(v2 + a2*t_out - 1/2*u*t_out^2); % vf

x1 = x0 + v0*t_in + 1/2*a0*t_in^2 + 1/6*u*t_in^3;
x2 = x1 + v1*dt + 1/2*a1*dt^2;
x3 = simplify(x2 + v2*t_out + 1/2*a2*t_out^2 - 1/6*u*t_out^3); % xf

% solve(vf - v3 == 0, dt)
dt_ = (vf - (- a0^2 - af^2 + 2*am^2 + 2*u*v0)/(2*u))/am;
x3_ = simplify(subs(x3-xf,dt,dt_));
err = (12*u^2*vf^2 - 12*u^2*v0^2 + 8*a0^3*am - 8*af^3*am - 3*a0^4 + 3*af^4 - 6*a0^2*am^2 + 6*af^2*am^2 + 12*a0^2*u*v0 + 12*am^2*u*v0 + 12*af^2*u*vf + 12*am^2*u*vf + 24*am*u^2*x0 - 24*am*u^2*xf - 24*a0*am*u*v0 - 24*af*am*u*vf)/(24*am*u^2)

% dt__ = []; vf__ = []; xf = [];
% am__ = 0.1:0.1:7.9;
% v0_ = 1.6083;
% vf_ = 1.9095;
% a0_ = -0.2785;
% af_ = -2.6832;
% 
% for kk = 1:length(am__)
%     dt__(kk) = subs(dt_,[u x0 v0 a0 am vf af],[20 0 v0_ a0_ am__(kk) vf_ af_]);
%     xf(kk) = double(subs(x3,[u x0 v0 a0 am vf af dt],[20 0 v0_ a0_ am__(kk) vf_ af_ dt__(kk)]));
% end
% 
% figure(110)
% plot(am__,xf)

%% acc case 1 :  2 -> 1
syms u am x0 v0 a0 xf vf af real
syms dt real;

t_in = (am - a0) / u;
t_out = (af - am) / u;

a1 = am; % a0 - u_in*t_in = a0 - u_in*(a0 - am)/u_in
a2 = am;
a3 = af; % am + u_out*t_out = am + u_out*(af - am)/u_out

v1 = v0 + a0*t_in + 1/2*u*t_in^2;
v2 = v1 + a1*dt;
v3 = simplify(v2 + a2*t_out + 1/2*u*t_out^2); % vf

x1 = x0 + v0*t_in + 1/2*a0*t_in^2 + 1/6*u*t_in^3;
x2 = x1 + v1*dt + 1/2*a1*dt^2;
x3 = simplify(x2 + v2*t_out + 1/2*a2*t_out^2 + 1/6*u*t_out^3); % xf

% dt_ = (vf - (- a0^2 + af^2 + 2*u*v0)/(2*u))/am;
% x3_ = simplify(subs(x3-xf,dt,dt_));
% err = -(12*u^2*v0^2 - 12*u^2*vf^2 - 8*a0^3*am + 8*af^3*am + 3*a0^4 - 3*af^4 + 6*a0^2*am^2 - 6*af^2*am^2 - 12*a0^2*u*v0 - 12*am^2*u*v0 + 12*af^2*u*vf + 12*am^2*u*vf - 24*am*u^2*x0 + 24*am*u^2*xf + 24*a0*am*u*v0 - 24*af*am*u*vf)/(24*am*u^2)

% 
% dt__ = []; vf__ = []; xf = [];
% am__ = 0.1:0.1:7.9;
% vf_ = 11.375;
% a0_ = -10;
% af_ = 15;
% 
% for kk = 1:length(am__)
%     dt__(kk) = subs(dt_,[u x0 v0 a0 am vf af],[20 3 2 a0_ am__(kk) vf_ af_]);
%     xf(kk) = double(subs(x3,[u x0 v0 a0 am vf af dt],[20 3 2 a0_ am__(kk) vf_ af_ dt__(kk)]));
% end
% 
% figure(110)
% plot(am__,xf)

%% acc case 1 :  1 -> 2
clear
syms u am x0 v0 a0 xf vf af real
syms dt real;

t_in = (a0 - am) / u;
t_out = (-af + am) / u;

a1 =  am; % a0 - u*t_in = a0 - u*(a0 - am)/u
a2 =  am;
a3 =  af; % am - u*t_out = am - u*(-af + am)/u

v1 = v0 + a0*t_in - 1/2*u*t_in^2;
v2 = v1 + a1*dt;
v3 = simplify(v2 + a2*t_out - 1/2*u*t_out^2); % vf

x1 = x0 + v0*t_in + 1/2*a0*t_in^2 - 1/6*u*t_in^3;
x2 = x1 + v1*dt + 1/2*a1*dt^2;
x3 = simplify(x2 + v2*t_out + 1/2*a2*t_out^2 - 1/6*u*t_out^3); % xf

% solve(v3 - vf == 0, dt)
dt_ = (vf - (a0^2 - af^2 + 2*u*v0)/(2*u))/am;
% x3_ = simplify(subs(x3-xf,dt,dt_));
% x3 - xf = -(12*u^2*v0^2 - 12*u^2*vf^2 - 8*a0^3*am + 8*af^3*am + 3*a0^4 - 3*af^4 + 6*a0^2*am^2 - 6*af^2*am^2 + 12*a0^2*u*v0 + 12*am^2*u*v0 - 12*af^2*u*vf - 12*am^2*u*vf - 24*am*u^2*x0 + 24*am*u^2*xf - 24*a0*am*u*v0 + 24*af*am*u*vf)/(24*am*u^2)

% dt__ = []; vf__ = []; xf = [];
% am__ = 4.10:.1:5.3;
% vf_ = 12.025;
% af_ = -7;
% 
% for kk = 1:length(am__)
%     dt__(kk) = subs(dt_,[u x0 v0 a0 am vf af],[20 3 2 10 am__(kk) vf_ af_]);
% %     vf__(kk) = subs(v3,[u x0 v0 a0 am vf af dt],[20 3 2 10 am__(kk) vf_ af_ dt__(kk)]);
%     xf(kk) = double(subs(x3,[u x0 v0 a0 am vf af dt],[20 3 2 10 am__(kk) vf_ af_ dt__(kk)]));
% end
% 
% figure(110)
% plot(am__,xf)

%% acc case 1 :  1 -> 1
syms u am x0 v0 a0 xf vf af real
syms dt real;

t_in = (a0 - am) / u;
t_out = (af - am) / u;

a1 = am; % a0 - u_in*t_in = a0 - u_in*(a0 - am)/u_in
a2 = am;
a3 = af; % am + u_out*t_out = am + u_out*(af - am)/u_out

v1 = v0 + a0*t_in - 1/2*u*t_in^2;
v2 = v1 + a1*dt;
v3 = simplify(v2 + a2*t_out + 1/2*u*t_out^2); % vf

x1 = x0 + v0*t_in + 1/2*a0*t_in^2 - 1/6*u*t_in^3;
x2 = x1 + v1*dt + 1/2*a1*dt^2;
x3 = simplify(x2 + v2*t_out + 1/2*a2*t_out^2 + 1/6*u*t_out^3); % xf

dt_ = (vf - (a0^2 + af^2 - 2*am^2 + 2*u*v0)/(2*u))/am;
% x3_ = simplify(subs(x3-xf,dt,dt_));
%x3 - xf = -(12*u^2*v0^2 - 12*u^2*vf^2 - 8*a0^3*am + 8*af^3*am + 3*a0^4 - 3*af^4 + 6*a0^2*am^2 - 6*af^2*am^2 + 12*a0^2*u*v0 + 12*am^2*u*v0 + 12*af^2*u*vf + 12*am^2*u*vf - 24*am*u^2*x0 + 24*am*u^2*xf - 24*a0*am*u*v0 - 24*af*am*u*vf)/(24*am*u^2)


% solve(v3 - vf == 0, dt)
% x3_ = simplify(subs(x3,dt,dt_));
% solve(x3_ - xf == 0, am)

% dt__ = []; vf__ = []; xf = [];
% am__ = 0.1:0.1:7.9;
% vf_ = 17.625;
% af_ = 15;

% for kk = 1:length(am__)
%     dt__(kk) = subs(dt_,[u x0 v0 a0 am vf af],[20 3 2 10 am__(kk) vf_ af_]);
% %     vf__(kk) = subs(v3,[u x0 v0 a0 am vf af dt],[20 3 2 10 am__(kk) vf_ vf_ dt__(kk)]);
%     xf(kk) = double(subs(x3,[u x0 v0 a0 am vf af dt],[20 3 2 10 am__(kk) vf_ af_ dt__(kk)]));
% end
% 
% figure(110)
% plot(am__,xf)

%%

% %% case 4-4
% a1 = a0 + t1*u;
% a2 = a1;
% a3 = a2 - (t3-t2)*u;
% a4 = a3;
% a5 = a4 - (tf-t4)*u;
% 
% v1 = v0 + a0*t1      + 1/2*u*t1^2;
% v2 = v1 + a1*(t2-t1);
% v3 = v2 + a2*(t3-t2) - 1/2*u*(t3-t2)^2;
% v4 = v3 + a3*(t4-t3);
% v5 = v4 + a4*(tf-t4) - 1/2*u*(tf-t4)^2;
% 
% x1 = x0 + v0*t1      + 1/2*a0*t1^2      + 1/6*u*t1^3;
% x2 = x1 + v1*(t2-t1) + 1/2*a1*(t2-t1)^2;
% x3 = x2 + v2*(t3-t2) + 1/2*a2*(t3-t2)^2 - 1/6*u*(t3-t2)^3;
% x4 = x3 + v3*(t4-t3) + 1/2*a3*(t4-t3)^2;
% x5 = x4 + v4*(tf-t4) + 1/2*a4*(tf-t4)^2 - 1/6*u*(tf-t4)^3;
% 
% t1_ = -(a0 - am)/u;
% t4_ = tf + (af + am)/u;
% 
% a5 = simplify(subs(a5,t4,t4_));
% v5 = simplify(subs(v5,t4,t4_));
% x5 = simplify(subs(x5,t4,t4_));
% a5 = simplify(subs(a5,t1,t1_));
% v5 = simplify(subs(v5,t1,t1_));
% x5 = simplify(subs(x5,t1,t1_));
% 
% A = simplify(diff(a5,t2));
% B = simplify(a5 - A*t2);
% simplify(a5 - A*t2 - B)
% t2_ = simplify((af - B) / A);
% 
% a5 = simplify(subs(a5,t2,t2_)); 
% v5 = simplify(subs(v5,t2,t2_));
% x5 = simplify(subs(x5,t2,t2_));
% 
% A = simplify(diff(v5,t3));
% B = simplify(v5 - A*t3);
% simplify(v5 - A*t3 - B)
% t3_ = simplify((vf - B) / A);
% x5 = simplify(subs(x5,t3,t3_));


%% case 4-3
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

t1_ = -(a0 - am)/u;
t4_ = tf - (af + am)/u;

a5 = simplify(subs(a5,t4,t4_));
v5 = simplify(subs(v5,t4,t4_));
x5 = simplify(subs(x5,t4,t4_));
a5 = simplify(subs(a5,t1,t1_));
v5 = simplify(subs(v5,t1,t1_));
x5 = simplify(subs(x5,t1,t1_));

A = simplify(diff(a5,t2));
B = simplify(a5 - A*t2);
simplify(a5 - A*t2 - B)
t2_ = simplify((af - B) / A);

a5 = simplify(subs(a5,t2,t2_)); 
v5 = simplify(subs(v5,t2,t2_));
x5 = simplify(subs(x5,t2,t2_));

A = simplify(diff(v5,t3));
B = simplify(v5 - A*t3);
simplify(v5 - A*t3 - B)
t3_ = simplify((vf - B) / A);
x5 = simplify(subs(x5,t3,t3_));

%%

% %% case 4-2
% a1 = a0 - t1*u;
% a2 = a1;
% a3 = a2 - (t3-t2)*u;
% a4 = a3;
% a5 = a4 - (tf-t4)*u;
% 
% v1 = v0 + a0*t1      - 1/2*u*t1^2;
% v2 = v1 + a1*(t2-t1);
% v3 = v2 + a2*(t3-t2) - 1/2*u*(t3-t2)^2;
% v4 = v3 + a3*(t4-t3);
% v5 = v4 + a4*(tf-t4) - 1/2*u*(tf-t4)^2;
% 
% x1 = x0 + v0*t1      + 1/2*a0*t1^2      - 1/6*u*t1^3;
% x2 = x1 + v1*(t2-t1) + 1/2*a1*(t2-t1)^2;
% x3 = x2 + v2*(t3-t2) + 1/2*a2*(t3-t2)^2 - 1/6*u*(t3-t2)^3;
% x4 = x3 + v3*(t4-t3) + 1/2*a3*(t4-t3)^2;
% x5 = x4 + v4*(tf-t4) + 1/2*a4*(tf-t4)^2 - 1/6*u*(tf-t4)^3;
% 
% t1_ = (a0 - am)/u;
% t4_ = tf + (af + am)/u;
% 
% a5 = simplify(subs(a5,t1,t1_)); 
% v5 = simplify(subs(v5,t1,t1_));
% x5 = simplify(subs(x5,t1,t1_));
% 
% a5 = simplify(subs(a5,t4,t4_)); 
% v5 = simplify(subs(v5,t4,t4_));
% x5 = simplify(subs(x5,t4,t4_));
% 
% A = simplify(diff(a5,t2));
% B = simplify(a5 - A*t2);
% simplify(a5 - A*t2 - B)
% t2_ = simplify((af - B) / A);
% 
% a5 = simplify(subs(a5,t2,t2_)); 
% v5 = simplify(subs(v5,t2,t2_));
% x5 = simplify(subs(x5,t2,t2_));
% 
% A = simplify(diff(v5,t3));
% B = simplify(v5 - A*t3);
% simplify(v5 - A*t3 - B)
% t3_ = simplify((vf - B) / A);
% 
% v5 = simplify(subs(v5,t3,t3_));
% x5 = simplify(subs(x5,t3,t3_));

% %% case 4-1
% a1 = a0 - t1*u;
% a2 = a1;
% a3 = a2 - (t3-t2)*u;
% a4 = a3;
% a5 = a4 + (tf-t4)*u;
% 
% v1 = v0 + a0*t1      - 1/2*u*t1^2;
% v2 = v1 + a1*(t2-t1);
% v3 = v2 + a2*(t3-t2) - 1/2*u*(t3-t2)^2;
% v4 = v3 + a3*(t4-t3);
% v5 = v4 + a4*(tf-t4) + 1/2*u*(tf-t4)^2;
% 
% x1 = x0 + v0*t1      + 1/2*a0*t1^2      - 1/6*u*t1^3;
% x2 = x1 + v1*(t2-t1) + 1/2*a1*(t2-t1)^2;
% x3 = x2 + v2*(t3-t2) + 1/2*a2*(t3-t2)^2 - 1/6*u*(t3-t2)^3;
% x4 = x3 + v3*(t4-t3) + 1/2*a3*(t4-t3)^2;
% x5 = x4 + v4*(tf-t4) + 1/2*a4*(tf-t4)^2 + 1/6*u*(tf-t4)^3;
% 
% t1_ = (a0 - am)/u;
% t4_ = tf - (af + am)/u;
% 
% a5 = simplify(subs(a5,t1,t1_)); 
% v5 = simplify(subs(v5,t1,t1_));
% x5 = simplify(subs(x5,t1,t1_));
% 
% a5 = simplify(subs(a5,t4,t4_)); 
% v5 = simplify(subs(v5,t4,t4_));
% x5 = simplify(subs(x5,t4,t4_));
% 
% A = simplify(diff(a5,t2));
% B = simplify(a5 - A*t2);
% simplify(a5 - A*t2 - B)
% t2_ = simplify((af - B) / A);
% 
% a5 = simplify(subs(a5,t2,t2_)); 
% v5 = simplify(subs(v5,t2,t2_));
% x5 = simplify(subs(x5,t2,t2_));
% 
% A = simplify(diff(v5,t3));
% B = simplify(v5 - A*t3);
% simplify(v5 - A*t3 - B)
% t3_ = simplify((vf - B) / A);
% 
% v5 = simplify(subs(v5,t3,t3_));
% x5 = simplify(subs(x5,t3,t3_));


% %% case 3-4
% a1 = a0 + t1*u;
% a2 = a1;
% a3 = a2 - (t3-t2)*u;
% a4 = a3;
% a5 = a4 - (tf-t4)*u;
% 
% v1 = v0 + a0*t1      + 1/2*u*t1^2;
% v2 = v1 + a1*(t2-t1);
% v3 = v2 + a2*(t3-t2) - 1/2*u*(t3-t2)^2;
% v4 = v3 + a3*(t4-t3);
% v5 = v4 + a4*(tf-t4) - 1/2*u*(tf-t4)^2;
% 
% x1 = x0 + v0*t1      + 1/2*a0*t1^2      + 1/6*u*t1^3;
% x2 = x1 + v1*(t2-t1) + 1/2*a1*(t2-t1)^2;
% x3 = x2 + v2*(t3-t2) + 1/2*a2*(t3-t2)^2 - 1/6*u*(t3-t2)^3;
% x4 = x3 + v3*(t4-t3) + 1/2*a3*(t4-t3)^2;
% x5 = x4 + v4*(tf-t4) + 1/2*a4*(tf-t4)^2 - 1/6*u*(tf-t4)^3;
% 
% a5 = simplify(subs(a5,t2,t1)); 
% v5 = simplify(subs(v5,t2,t1));
% x5 = simplify(subs(x5,t2,t1));
% 
% t4_ = tf + (af + am)/u;
% 
% a5 = simplify(subs(a5,t4,t4_)); 
% v5 = simplify(subs(v5,t4,t4_));
% x5 = simplify(subs(x5,t4,t4_));
% 
% A = simplify(diff(a5,t1));
% B = simplify(a5 - A*t1);
% simplify(a5 - A*t1 - B)
% t1_ = simplify((af - B) / A);
% 
% a5 = simplify(subs(a5,t1,t1_)); 
% v5 = simplify(subs(v5,t1,t1_));
% x5 = simplify(subs(x5,t1,t1_));
% 
% % A = simplify(diff(a5,t2));
% % B = simplify(a5 - A*t2);
% % simplify(a5 - A*t2 - B)
% % t2_ = simplify((af - B) / A);
% % 
% % a5 = simplify(subs(a5,t2,t2_)); 
% % v5 = simplify(subs(v5,t2,t2_));
% % x5 = simplify(subs(x5,t2,t2_));
% 
% dv5 = simplify(diff(v5,t3));
% ddv5 = simplify(diff(dv5,t3));
% 
% A = simplify(ddv5/2);
% B = simplify(dv5 - 2*A*t3);
% C = simplify(v5 - A*t3^2 - B*t3);
% simplify(v5 - A*t3^2 - B*t3 - C)
% % A*t3^2 + B*t3 + C - vf = 0
% C = C - vf;


% %% case 3-3
% a1 = a0 + t1*u;
% a2 = a1;
% a3 = a2 - (t3-t2)*u;
% a4 = a3;
% a5 = a4 + (tf-t4)*u;
% 
% v1 = v0 + a0*t1      + 1/2*u*t1^2;
% v2 = v1 + a1*(t2-t1);
% v3 = v2 + a2*(t3-t2) - 1/2*u*(t3-t2)^2;
% v4 = v3 + a3*(t4-t3);
% v5 = v4 + a4*(tf-t4) + 1/2*u*(tf-t4)^2;
% 
% x1 = x0 + v0*t1      + 1/2*a0*t1^2      + 1/6*u*t1^3;
% x2 = x1 + v1*(t2-t1) + 1/2*a1*(t2-t1)^2;
% x3 = x2 + v2*(t3-t2) + 1/2*a2*(t3-t2)^2 - 1/6*u*(t3-t2)^3;
% x4 = x3 + v3*(t4-t3) + 1/2*a3*(t4-t3)^2;
% x5 = x4 + v4*(tf-t4) + 1/2*a4*(tf-t4)^2 + 1/6*u*(tf-t4)^3;
% 
% a5 = simplify(subs(a5,t2,t1)); 
% v5 = simplify(subs(v5,t2,t1));
% x5 = simplify(subs(x5,t2,t1));
% 
% t4_ = (-am - af)/u + tf;
% 
% a5 = simplify(subs(a5,t4,t4_)); 
% v5 = simplify(subs(v5,t4,t4_));
% x5 = simplify(subs(x5,t4,t4_));
% 
% A = simplify(diff(a5,t1));
% B = simplify(a5 - A*t1);
% simplify(a5 - A*t1 - B)
% t1_ = simplify((af - B) / A);
% 
% a5 = simplify(subs(a5,t1,t1_)); 
% v5 = simplify(subs(v5,t1,t1_));
% x5 = simplify(subs(x5,t1,t1_));
% 
% A = simplify(diff(a5,t2));
% B = simplify(a5 - A*t2);
% simplify(a5 - A*t2 - B)
% t2_ = simplify((af - B) / A);
% 
% a5 = simplify(subs(a5,t2,t2_)); 
% v5 = simplify(subs(v5,t2,t2_));
% x5 = simplify(subs(x5,t2,t2_));
% 
% dv5 = simplify(diff(v5,t3));
% ddv5 = simplify(diff(dv5,t3));
% 
% A = simplify(ddv5/2);
% B = simplify(dv5 - 2*A*t3);
% C = simplify(v5 - A*t3^2 - B*t3);
% simplify(v5 - A*t3^2 - B*t3 - C)
% % A*t3^2 + B*t3 + C - vf = 0
% C = C - vf;

% %% case 3-2
% a1 = a0;
% a2 = a1;
% a3 = a2 - (t3-t2)*u;
% a4 = a3;
% a5 = a4 - (tf-t4)*u;
% 
% v1 = v0 + a0*t1      ;
% v2 = v1 + a1*(t2-t1);
% v3 = v2 + a2*(t3-t2) - 1/2*u*(t3-t2)^2;
% v4 = v3 + a3*(t4-t3);
% v5 = v4 + a4*(tf-t4) - 1/2*u*(tf-t4)^2;
% 
% x1 = x0 + v0*t1      + 1/2*a0*t1^2      ;
% x2 = x1 + v1*(t2-t1) + 1/2*a1*(t2-t1)^2;
% x3 = x2 + v2*(t3-t2) + 1/2*a2*(t3-t2)^2 - 1/6*u*(t3-t2)^3;
% x4 = x3 + v3*(t4-t3) + 1/2*a3*(t4-t3)^2;
% x5 = x4 + v4*(tf-t4) + 1/2*a4*(tf-t4)^2 - 1/6*u*(tf-t4)^3;
% 
% a5 = simplify(subs(a5,t1,0)); 
% a5 = simplify(subs(a5,t2,0)); 
% v5 = simplify(subs(v5,t1,0));
% v5 = simplify(subs(v5,t2,0));
% x5 = simplify(subs(x5,t1,0));
% x5 = simplify(subs(x5,t2,0));
% 
% A = simplify(diff(a3,t3));
% B = simplify(a3 - A*t3);
% simplify(a3 - A*t3 - B)
% t3_ = simplify((-am - B) / A);
% t3_ = subs(t3_,t2,0);
% 
% a5_ = subs(a5,t3,t3_);
% 
% A = simplify(diff(a5_,t4));
% B = simplify(a5_ - A*t4);
% simplify(a5_ - A*t4 - B)
% t4_ = simplify((af - B) / A);
% 
% x5 = subs(x5,[t3 t4],[t3_ t4_]);


% %% case 3-1
% a1 = a0;
% a2 = a1;
% a3 = a2 - (t3-t2)*u;
% a4 = a3;
% a5 = a4 + (tf-t4)*u;
% 
% v1 = v0 + a0*t1      ;
% v2 = v1 + a1*(t2-t1);
% v3 = v2 + a2*(t3-t2) - 1/2*u*(t3-t2)^2;
% v4 = v3 + a3*(t4-t3);
% v5 = v4 + a4*(tf-t4) + 1/2*u*(tf-t4)^2;
% 
% x1 = x0 + v0*t1      + 1/2*a0*t1^2      ;
% x2 = x1 + v1*(t2-t1) + 1/2*a1*(t2-t1)^2;
% x3 = x2 + v2*(t3-t2) + 1/2*a2*(t3-t2)^2 - 1/6*u*(t3-t2)^3;
% x4 = x3 + v3*(t4-t3) + 1/2*a3*(t4-t3)^2;
% x5 = x4 + v4*(tf-t4) + 1/2*a4*(tf-t4)^2 + 1/6*u*(tf-t4)^3;
% 
% a5 = simplify(subs(a5,t1,0)); 
% a5 = simplify(subs(a5,t2,0)); 
% v5 = simplify(subs(v5,t1,0));
% v5 = simplify(subs(v5,t2,0));
% x5 = simplify(subs(x5,t1,0));
% x5 = simplify(subs(x5,t2,0));
% 
% A = simplify(diff(a3,t3));
% B = simplify(a3 - A*t3);
% simplify(a3 - A*t3 - B)
% t3_ = simplify((-am - B) / A);
% t3_ = subs(t3_,t2,0);
% 
% a5_ = subs(a5,t3,t3_);
% 
% A = simplify(diff(a5_,t4));
% B = simplify(a5_ - A*t4);
% simplify(a5_ - A*t4 - B)
% t4_ = simplify((af - B) / A);
% 
% x5 = subs(x5,[t3 t4],[t3_ t4_]);

% %% case 2-3
% a1 = a0 + t1*u;
% a2 = a1;
% a3 = a2 - (t3-t2)*u;
% a4 = a3;
% a5 = a4 + (tf-t4)*u;
% 
% v1 = v0 + a0*t1      + 1/2*u*t1^2;
% v2 = v1 + a1*(t2-t1);
% v3 = v2 + a2*(t3-t2) - 1/2*u*(t3-t2)^2;
% v4 = v3 + a3*(t4-t3);
% v5 = v4 + a4*(tf-t4) + 1/2*u*(tf-t4)^2;
% 
% x1 = x0 + v0*t1      + 1/2*a0*t1^2      + 1/6*u*t1^3;
% x2 = x1 + v1*(t2-t1) + 1/2*a1*(t2-t1)^2;
% x3 = x2 + v2*(t3-t2) + 1/2*a2*(t3-t2)^2 - 1/6*u*(t3-t2)^3;
% x4 = x3 + v3*(t4-t3) + 1/2*a3*(t4-t3)^2;
% x5 = x4 + v4*(tf-t4) + 1/2*a4*(tf-t4)^2 + 1/6*u*(tf-t4)^3;
% 
% a5 = simplify(subs(a5,t4,t3)); 
% v5 = simplify(subs(v5,t4,t3));
% x5 = simplify(subs(x5,t4,t3));
% 
% A = simplify(diff(a1,t1));
% B = simplify(a1 - A*t1);
% simplify(a1 - A*t1 - B)
% t1_ = simplify((am - B) / A);
% 
% a5 = simplify(subs(a5,t1,t1_)); 
% v5 = simplify(subs(v5,t1,t1_));
% x5 = simplify(subs(x5,t1,t1_));
% 
% A = simplify(diff(a5,t2));
% B = simplify(a5 - A*t2);
% simplify(a5 - A*t2 - B)
% t2_ = simplify((af - B) / A);
% 
% a5 = simplify(subs(a5,t2,t2_)); 
% v5 = simplify(subs(v5,t2,t2_));
% x5 = simplify(subs(x5,t2,t2_));
% 
% dv5 = simplify(diff(v5,t3));
% ddv5 = simplify(diff(dv5,t3));
% 
% A = simplify(ddv5/2);
% B = simplify(dv5 - 2*A*t3);
% C = simplify(v5 - A*t3^2 - B*t3);
% simplify(v5 - A*t3^2 - B*t3 - C)
% % A*t3^2 + B*t3 + C - vf = 0
% C = C - vf;

% %% case 2-2
% a1 = a0 - t1*u;
% a2 = a1;
% a3 = a2;
% a4 = a3;
% a5 = a4 - (tf-t2)*u;
% 
% v1 = v0 + a0*t1      - 1/2*u*t1^2;
% v2 = v1 + a1*(t2-t1);
% v3 = v2 + a2*(t3-t2);
% v4 = v3 + a3*(t4-t3);
% v5 = v4 + a4*(tf-t4) - 1/2*u*(tf-t2)^2;
% 
% x1 = x0 + v0*t1      + 1/2*a0*t1^2      - 1/6*u*t1^3;
% x2 = x1 + v1*(t2-t1) + 1/2*a1*(t2-t1)^2;
% x3 = x2 + v2*(t3-t2) + 1/2*a2*(t3-t2)^2;
% x4 = x3 + v3*(t4-t3) + 1/2*a3*(t4-t3)^2;
% x5 = x4 + v4*(tf-t4) + 1/2*a4*(tf-t4)^2 - 1/6*u*(tf-t2)^3;
% 
% A = simplify(diff(a1,t1));
% B = simplify(a1 - A*t1);
% simplify(a1 - A*t1 - B)
% t1_ = simplify((am - B) / A);
% 
% 
% a5 = simplify(subs(a5,t1,t1_));
% % a5 = simplify(subs(a5,t2,0));
% a5 = simplify(subs(a5,t3,t2));
% a5 = simplify(subs(a5,t4,t2));
% 
% v5 = simplify(subs(v5,t1,t1_));
% % v5 = simplify(subs(v5,t2,0));
% v5 = simplify(subs(v5,t3,t2));
% v5 = simplify(subs(v5,t4,t2));
% 
% x5 = simplify(subs(x5,t1,t1_));
% % x5 = simplify(subs(x5,t2,0));
% x5 = simplify(subs(x5,t3,t2));
% x5 = simplify(subs(x5,t4,t2));
% 
% A = simplify(diff(a5,t2));
% B = simplify(a5 - A*t2);
% simplify(a5 - A*t2 - B)
% t2 = simplify((af - B) / A);

% %% case 2-1
% a1 = a0 - t1*u;
% a2 = a1;
% a3 = a2 - (t3-t2)*u;
% a4 = a3;
% a5 = a4 + (tf-t4)*u;
% 
% v1 = v0 + a0*t1      - 1/2*u*t1^2;
% v2 = v1 + a1*(t2-t1);
% v3 = v2 + a2*(t3-t2) - 1/2*u*(t3-t2)^2;
% v4 = v3 + a3*(t4-t3);
% v5 = v4 + a4*(tf-t4) + 1/2*u*(tf-t4)^2;
% 
% x1 = x0 + v0*t1      + 1/2*a0*t1^2      - 1/6*u*t1^3;
% x2 = x1 + v1*(t2-t1) + 1/2*a1*(t2-t1)^2;
% x3 = x2 + v2*(t3-t2) + 1/2*a2*(t3-t2)^2 - 1/6*u*(t3-t2)^3;
% x4 = x3 + v3*(t4-t3) + 1/2*a3*(t4-t3)^2;
% x5 = x4 + v4*(tf-t4) + 1/2*a4*(tf-t4)^2 + 1/6*u*(tf-t4)^3;
% 
% A = simplify(diff(a1,t1));
% B = simplify(a1 - A*t1);
% simplify(a1 - A*t1 - B)
% t1_ = simplify((am - B) / A);
% 
% a5 = simplify(subs(a5,t1,t1_));
% a5 = simplify(subs(a5,t4,t3));
% 
% v5 = simplify(subs(v5,t1,t1_));
% v5 = simplify(subs(v5,t4,t3));
% 
% x5 = simplify(subs(x5,t1,t1_));
% x5 = simplify(subs(x5,t4,t3));
% 
% A = simplify(diff(a5,t2));
% B = simplify(a5 - A*t2);
% simplify(a5 - A*t2 - B)
% t2_ = simplify((af - B) / A);
% 
% v5 = simplify(subs(v5,t2,t2_));
% 
% dv5 = simplify(diff(v5,t3));
% ddv5 = simplify(diff(dv5,t3));
% 
% A = simplify(ddv5/2);
% B = simplify(dv5 - 2*A*t3);
% C = simplify(v5 - A*t3^2 - B*t3);
% simplify(v5 - A*t3^2 - B*t3 - C)
% % A*t3^2 + B*t3 + C - vf = 0
% C = C - vf;
% 
% %% case 1-4 = 1-1

% %% case 1-3
% 
% a1 = a0 + t1*u;
% a2 = a1;
% a3 = a2 - (t3-t2)*u;
% a4 = a3;
% a5 = a4 + (tf-t4)*u;
% 
% v1 = v0 + a0*t1      + 1/2*u*t1^2;
% v2 = v1 + a1*(t2-t1);
% v3 = v2 + a2*(t3-t2) - 1/2*u*(t3-t2)^2;
% v4 = v3 + a3*(t4-t3);
% v5 = v4 + a4*(tf-t4) + 1/2*u*(tf-t4)^2;
% 
% x1 = x0 + v0*t1      + 1/2*a0*t1^2      + 1/6*u*t1^3;
% x2 = x1 + v1*(t2-t1) + 1/2*a1*(t2-t1)^2;
% x3 = x2 + v2*(t3-t2) + 1/2*a2*(t3-t2)^2 - 1/6*u*(t3-t2)^3;
% x4 = x3 + v3*(t4-t3) + 1/2*a3*(t4-t3)^2;
% x5 = x4 + v4*(tf-t4) + 1/2*a4*(tf-t4)^2 + 1/6*u*(tf-t4)^3;
% 
% % a5 = simplify(subs(a5,t1,t2));
% a5 = simplify(subs(a5,t2,t1));
% % a5 = simplify(subs(a5,t3,tf));
% a5 = simplify(subs(a5,t4,t3));
% 
% % v5 = simplify(subs(v5,t1,0));
% v5 = simplify(subs(v5,t2,t1));
% % v5 = simplify(subs(v5,t3,tf));
% v5 = simplify(subs(v5,t4,t3));
% 
% % x5 = simplify(subs(x5,t1,0));
% x5 = simplify(subs(x5,t2,t1));
% % x5 = simplify(subs(x5,t3,tf));
% x5 = simplify(subs(x5,t4,t3));
% 
% A = simplify(diff(a5,t1));
% B = simplify(a5 - A*t1);
% simplify(a5 - A*t1 - B)
% t1_ = simplify((af - B) / A);
% 
% v5_ = simplify(subs(v5,t1,t1_));
% A = simplify(diff(v5_,t3));
% B = simplify(v5_ - A*t3);
% simplify(v5_ - A*t3 - B)
% t3_ = simplify((vf - B) / A);
% x5_ = simplify(subs(x5,[t1 t3],[t1_ t3_]));
% v5_ = simplify(subs(v5,[t1 t3],[t1_ t3_]));

%% case 1-2

% a1 = a0 + t1*u;
% a2 = a1;
% a3 = a2 - (t3-t2)*u;
% a4 = a3;
% a5 = a4 + (tf-t4)*u;
% 
% v1 = v0 + a0*t1      + 1/2*u*t1^2;
% v2 = v1 + a1*(t2-t1);
% v3 = v2 + a2*(t3-t2) - 1/2*u*(t3-t2)^2;
% v4 = v3 + a3*(t4-t3);
% v5 = v4 + a4*(tf-t4) + 1/2*u*(tf-t4)^2;
% 
% x1 = x0 + v0*t1      + 1/2*a0*t1^2      + 1/6*u*t1^3;
% x2 = x1 + v1*(t2-t1) + 1/2*a1*(t2-t1)^2;
% x3 = x2 + v2*(t3-t2) + 1/2*a2*(t3-t2)^2 - 1/6*u*(t3-t2)^3;
% x4 = x3 + v3*(t4-t3) + 1/2*a3*(t4-t3)^2;
% x5 = x4 + v4*(tf-t4) + 1/2*a4*(tf-t4)^2 + 1/6*u*(tf-t4)^3;
% 
% a5 = simplify(subs(a5,t1,0));
% a5 = simplify(subs(a5,t2,0));
% a5 = simplify(subs(a5,t3,tf));
% a5 = simplify(subs(a5,t4,tf));
% 
% v5 = simplify(subs(v5,t1,0));
% v5 = simplify(subs(v5,t2,0));
% v5 = simplify(subs(v5,t3,tf));
% v5 = simplify(subs(v5,t4,tf));
% 
% x5 = simplify(subs(x5,t1,0));
% x5 = simplify(subs(x5,t2,0));
% x5 = simplify(subs(x5,t3,tf));
% x5 = simplify(subs(x5,t4,tf));

%% case 1-1

% a1 = a0 - t4*u;
% a2 = a1 + 0;
% a3 = a2 - 0;
% a4 = a3 + 0;
% a5 = a4 + (tf-t4)*u;
% 
% v1 = v0 + a0*t1      - 1/2*u*t1^2;
% v2 = v1 + a1*(t2-t1);
% v3 = v2 + a2*(t3-t2) - 1/2*u*(t3-t2)^2;
% v4 = v3 + a3*(t4-t3);
% v5 = v4 + a4*(tf-t4) + 1/2*u*(tf-t4)^2;
% 
% x1 = x0 + v0*t1      + 1/2*a0*t1^2      - 1/6*u*t1^3;
% x2 = x1 + v1*(t2-t1) + 1/2*a1*(t2-t1)^2;
% x3 = x2 + v2*(t3-t2) + 1/2*a2*(t3-t2)^2 + 1/6*u*(t3-t2)^3;
% x4 = x3 + v3*(t4-t3) + 1/2*a3*(t4-t3)^2;
% x5 = x4 + v4*(tf-t4) + 1/2*a4*(tf-t4)^2 + 1/6*u*(tf-t4)^3;
% 
% a5 = simplify(subs(a5,t1,t4));
% a5 = simplify(subs(a5,t2,t4));
% a5 = simplify(subs(a5,t3,t4));
% 
% v5 = simplify(subs(v5,t1,t4));
% v5 = simplify(subs(v5,t2,t4));
% v5 = simplify(subs(v5,t3,t4));
% 
% x5 = simplify(subs(x5,t1,t4));
% x5 = simplify(subs(x5,t2,t4));
% x5 = simplify(subs(x5,t3,t4));
% 
% A = diff(a5,t4);
% B = a5 - A*t4;
% a5 - A*t4 - B % to check
% t4_ = simplify((af - B) / A);




