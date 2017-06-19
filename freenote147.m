clear all
close all
clc

syms t1 t2 t3 t4 tf real

u = -3;
amax = 2;

tt = [];

x0 = 3;
v0 = 5;
a0 = 0;

xf = 0;
vf = 0;
af = 0;

a1 = a0 + t1*u;
a2 = a1;
a3 = a2 - (t3 - t2)*u;
a4 = a3;
a5 = a4 + (tf - t4)*u;

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

%%
a5_ = subs(a5,[t2 t4],[t1 t3]);
v5_ = subs(v5,[t2 t4],[t1 t3]);
x5_ = subs(x5,[t2 t4],[t1 t3]);

[solt1,solt3,soltf] = solve([a5_ == af,v5_ == vf,x5_ == xf],[t1,t3,tf]);

for k=1:length(solt1)
   tt = [tt;double([solt1(k) solt1(k) solt3(k) solt3(k) soltf(k)])]; 
end

%%
a5_ = subs(a5,[t1 t4],[(-amax-a0)/u,t3]);
v5_ = subs(v5,[t1 t4],[(-amax-a0)/u,t3]);
x5_ = subs(x5,[t1 t4],[(-amax-a0)/u,t3]);

[solt2,solt3,soltf] = solve([a5_ == af,v5_ == vf,x5_ == xf],[t2,t3,tf]);

for k=1:length(solt2)
   tt = [tt;double([(-amax-a0)/u solt2(k) solt3(k) solt3(k) soltf(k)])]; 
end

%%
a5_ = subs(a5,[t2 tf],[t1,t4 + (af - amax)/u]);
v5_ = subs(v5,[t2 tf],[t1,t4 + (af - amax)/u]);
x5_ = subs(x5,[t2 tf],[t1,t4 + (af - amax)/u]);

[solt1,solt3,solt4] = solve([a5_ == af,v5_ == vf,x5_ == xf],[t1,t3,t4]);

for k=1:length(solt1)
   tt = [tt;double([solt1(k) solt1(k) solt3(k) solt4(k) solt4(k)+(af - amax)/u])]; 
end

%%
a5_ = subs(a5,[t1 t3],[(-amax-a0)/u,t2 + (-amax-amax)/u]);
v5_ = subs(v5,[t1 t3],[(-amax-a0)/u,t2 + (-amax-amax)/u]);
x5_ = subs(x5,[t1 t3],[(-amax-a0)/u,t2 + (-amax-amax)/u]);

[solt2,solt4,soltf] = solve([a5_ == af,v5_ == vf,x5_ == xf],[t2,t4,tf]);

for k=1:length(solt2)
   tt = [tt;double([(-amax-a0)/u solt2(k) solt2(k)+(-amax-amax)/u solt4(k) soltf(k)])]; 
end

%%
for k=1:size(tt,1)
   a5 = double(subs(a5,[t1 t2 t3 t4 tf],real(tt(k,:))));
   v5 = double(subs(v5,[t1 t2 t3 t4 tf],real(tt(k,:))));
   x5 = double(subs(x5,[t1 t2 t3 t4 tf],real(tt(k,:))));   
   disp([num2str(x5),' ',num2str(v5),' ',num2str(a5)]);
end

