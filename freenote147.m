clear all
close all
clc

syms t1 t2 t3 t4 tf real

u = 8;
amax = -5;

tt = [];

x0 = -3;
v0 = -5;
a0 = 0;

xf = -4;
vf = -10;
af = 0;

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
chad = diff(real(tt)')';

idx = [];
for k=1:size(chad,1)
    temp = sum(find(chad(k,:) >= 0));
    if temp == 10
        idx = [idx;k];
    end
end

if idx ~= 0
    for aa = 1:length(idx)
        ttt = real(tt(idx(aa),:));
        
        x0_ = x0;
        v0_ = v0;
        a0_ = a0;
        a1_ = subs(a1,[t1 t2 t3 t4 tf],ttt);
        a2_ = subs(a2,[t1 t2 t3 t4 tf],ttt);
        a3_ = subs(a3,[t1 t2 t3 t4 tf],ttt);
        a4_ = subs(a4,[t1 t2 t3 t4 tf],ttt);
        a5_ = subs(a5,[t1 t2 t3 t4 tf],ttt);
        
        %%
        tf_ = ttt(end);
        tpiece = double([0 ttt]);
        apiece = double([a0_ a1_ a2_ a3_ a4_ a5_]);
        temp = find(abs(apiece > abs(amax)+1e-5));
        
        if(isempty(temp))        
            sim('simMinTimeTraj');
            xTraj = traj;
            figure(2)
            for k=1:3
                subplot(3,1,k)
                plot(xTraj.Time,xTraj.Data(:,k))
                grid on
                box on
            end
        end
    end
end













