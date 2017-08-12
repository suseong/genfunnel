clear
clc

xinit =  [10*(rand(1,2)-0.5) 15*(rand(1)-0.5)];
xfinal = [10*(rand(1,2)-0.5) 15*(rand(1)-0.5)];

yinit =  [10*(rand(1,2)-0.5) 15*(rand(1)-0.5)];
yfinal = [10*(rand(1,2)-0.5) 15*(rand(1)-0.5)];

zinit =  [10*(rand(1,2)-0.5) 15*(rand(1)-0.5)];
zfinal = [10*(rand(1,2)-0.5) 15*(rand(1)-0.5)];

init{1} =  xinit;
init{2} =  yinit;
init{3} =  zinit;
final{1} = xfinal;
final{2} = yfinal;
final{3} = zfinal;

% xinit  = [  0];     init{1} =  xinit;
% yinit  = [  0];     init{2} =  yinit;
% zinit  = [  0];     init{3} =  zinit;
% xfinal = [  0];    final{1} = xfinal;
% yfinal = [  0];    final{2} = yfinal;
% zfinal = [  0];    final{3} = zfinal;

[tsq,pos,acc,act] = mintime_traj(init,final);

% xacc = acc{1}; yacc = acc{2}; zacc = acc{3};
% figure(151);clf;
% totalAcc = sqrt(xacc.^2 + yacc.^2 + (zacc + 9.8).^2);

%%
clear all
close all
clc

u = 20;
amax = 3;

x0 = 0;
v0 = 0;
a0 = -17;

t1 = 1.0;
t2 = t1 + 0;
t3 = t2 + 0.2;
t4 = t3 + 0;
tf = t4 + 1;

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

x5
disp([num2str(t1),' ',num2str(t2),' ',num2str(t3),' ',num2str(t4),' ',num2str(tf)])
[a,b] = calcX5_4([20 3],[0 0 -17],[x5 v5 19],tf)





