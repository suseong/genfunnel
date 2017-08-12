clear all
close all
% clc

syms u t1 t2 t3 t4 tf real

u_ = 20;
amax = 3;

tt = [];

init = [ 3.0033   -0.4620   -0.6761];
final = [ 3.2531   -4.1653   -3.6683];
x0 = init(1); v0 = init(2); a0 = init(3);
xf = final(1);vf = final(2);af = final(3);

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
a5_ = subs(a5,[t2 t4 u],[t1 t3 u_]);
v5_ = subs(v5,[t2 t4 u],[t1 t3 u_]);
x5_ = subs(x5,[t2 t4 u],[t1 t3 u_]);

[solt1,solt3,soltf] = solve([a5_ == af,v5_ == vf,x5_ == xf],[t1,t3,tf]);

for k=1:length(solt1)
   tt = [tt;double([solt1(k) solt1(k) solt3(k) solt3(k) soltf(k) u_ 1.1])]; 
end

a5_ = subs(a5,[t2 t4 u],[t1 t3 -u_]);
v5_ = subs(v5,[t2 t4 u],[t1 t3 -u_]);
x5_ = subs(x5,[t2 t4 u],[t1 t3 -u_]);

[solt1,solt3,soltf] = solve([a5_ == af,v5_ == vf,x5_ == xf],[t1,t3,tf]);

for k=1:length(solt1)
   tt = [tt;double([solt1(k) solt1(k) solt3(k) solt3(k) soltf(k) -u_ 1.2])]; 
end

%% 
% u > 0, 
a5_ = subs(a5,[t1 t4 u],[(amax-a0)/u_ t3 u_]);
v5_ = subs(v5,[t1 t4 u],[(amax-a0)/u_ t3 u_]);
x5_ = subs(x5,[t1 t4 u],[(amax-a0)/u_ t3 u_]);

[solt2,solt3,soltf] = solve([a5_ == af,v5_ == vf,x5_ == xf],[t2,t3,tf]);

for k=1:length(solt2)
   tt = [tt;double([(amax-a0)/u_ solt2(k) solt3(k) solt3(k) soltf(k) u_ 2.1])]; 
end

% u < 0
a5_ = subs(a5,[t1 t4 u],[(-amax-a0)/-u_,t3,-u_]);
v5_ = subs(v5,[t1 t4 u],[(-amax-a0)/-u_,t3,-u_]);
x5_ = subs(x5,[t1 t4 u],[(-amax-a0)/-u_,t3,-u_]);

[solt2,solt3,soltf] = solve([a5_ == af,v5_ == vf,x5_ == xf],[t2,t3,tf]);

for k=1:length(solt2)
   tt = [tt;double([(-amax-a0)/-u_ solt2(k) solt3(k) solt3(k) soltf(k) -u_ 2.2])]; 
end

%%
% u > 0
a5_ = subs(a5,[t2 tf u],[t1 t4 + (af+amax)/u_ u_]);
v5_ = subs(v5,[t2 tf u],[t1 t4 + (af+amax)/u_ u_]);
x5_ = subs(x5,[t2 tf u],[t1 t4 + (af+amax)/u_ u_]);

[solt1,solt3,solt4] = solve([a5_ == af,v5_ == vf,x5_ == xf],[t1,t3,t4]);

for k=1:length(solt1)
   tt = [tt;double([solt1(k) solt1(k) solt3(k) solt4(k) solt4(k)+(af+amax)/u_ u_ 3.1])]; 
end

% u < 0
a5_ = subs(a5,[t2 tf u],[t1,t4 + (af-amax)/-u_ -u_]);
v5_ = subs(v5,[t2 tf u],[t1,t4 + (af-amax)/-u_ -u_]);
x5_ = subs(x5,[t2 tf u],[t1,t4 + (af-amax)/-u_ -u_]);

[solt1,solt3,solt4] = solve([a5_ == af,v5_ == vf,x5_ == xf],[t1,t3,t4]);

for k=1:length(solt1)
   tt = [tt;double([solt1(k) solt1(k) solt3(k) solt4(k) solt4(k)+(af-amax)/-u_ -u_ 3.2])]; 
end

%%
% u > 0
a5_ = subs(a5,[t1 t3 u],[(amax-a0)/u_ t2 + 2*amax/u_ u_]);
v5_ = subs(v5,[t1 t3 u],[(amax-a0)/u_ t2 + 2*amax/u_ u_]);
x5_ = subs(x5,[t1 t3 u],[(amax-a0)/u_ t2 + 2*amax/u_ u_]);

[solt2,solt4,soltf] = solve([a5_ == af,v5_ == vf,x5_ == xf],[t2,t4,tf]);

for k=1:length(solt2)
   tt = [tt;double([(amax-a0)/u_ solt2(k) solt2(k)+2*amax/u_ solt4(k) soltf(k) u_ 4.1])]; 
end

% u < 0
a5_ = subs(a5,[t1 t3 u],[(-amax-a0)/-u_ t2 + 2*amax/u_ -u_]);
v5_ = subs(v5,[t1 t3 u],[(-amax-a0)/-u_ t2 + 2*amax/u_ -u_]);
x5_ = subs(x5,[t1 t3 u],[(-amax-a0)/-u_ t2 + 2*amax/u_ -u_]);

[solt2,solt4,soltf] = solve([a5_ == af,v5_ == vf,x5_ == xf],[t2,t4,tf]);

for k=1:length(solt2)
   tt = [tt;double([(-amax-a0)/-u_ solt2(k) solt2(k)+2*amax/u_ solt4(k) soltf(k) -u_ 4.2])]; 
end

%%
a1 = a0 - t1*u;
a2 = a1;
a3 = a2 - (t3-t2)*u;
a4 = a3;
a5 = a4 + (tf-t4)*u;

v1 = v0 + a0*t1      - 1/2*u*t1^2;
v2 = v1 + a1*(t2-t1);
v3 = v2 + a2*(t3-t2) - 1/2*u*(t3-t2)^2;
v4 = v3 + a3*(t4-t3);
v5 = v4 + a4*(tf-t4) + 1/2*u*(tf-t4)^2;

x1 = x0 + v0*t1      + 1/2*a0*t1^2      - 1/6*u*t1^3;
x2 = x1 + v1*(t2-t1) + 1/2*a1*(t2-t1)^2;
x3 = x2 + v2*(t3-t2) + 1/2*a2*(t3-t2)^2 - 1/6*u*(t3-t2)^3;
x4 = x3 + v3*(t4-t3) + 1/2*a3*(t4-t3)^2;
x5 = x4 + v4*(tf-t4) + 1/2*a4*(tf-t4)^2 + 1/6*u*(tf-t4)^3;

%%
% u > 0
a5_ = subs(a5,[t1 t3 u],[(a0-amax)/u_ t4 u_]);
v5_ = subs(v5,[t1 t3 u],[(a0-amax)/u_ t4 u_]);
x5_ = subs(x5,[t1 t3 u],[(a0-amax)/u_ t4 u_]);

[solt2,solt4,soltf] = solve([a5_ == af,v5_ == vf,x5_ == xf],[t2,t4,tf]);

for k=1:length(solt2)
   tt = [tt;double([(a0-amax)/u_ solt2(k) solt4(k) solt4(k) soltf(k) u_ 5.1])]; 
end

% u < 0
a5_ = subs(a5,[t1 t3 u],[(a0+amax)/-u_ t4 -u_]);
v5_ = subs(v5,[t1 t3 u],[(a0+amax)/-u_ t4 -u_]);
x5_ = subs(x5,[t1 t3 u],[(a0+amax)/-u_ t4 -u_]);

[solt2,solt4,soltf] = solve([a5_ == af,v5_ == vf,x5_ == xf],[t2,t4,tf]);

for k=1:length(solt2)
   tt = [tt;double([(a0+amax)/-u_ solt2(k) solt4(k) solt4(k) soltf(k) -u_ 5.2])]; 
end

%%
% u > 0
a5_ = subs(a5,[t1 t4 u],[(a0-amax)/u_ (af+amax)/u_ u_]);
v5_ = subs(v5,[t1 t4 u],[(a0-amax)/u_ (af+amax)/u_ u_]);
x5_ = subs(x5,[t1 t4 u],[(a0-amax)/u_ (af+amax)/u_ u_]);

[solt2,solt3,soltf] = solve([a5_ == af,v5_ == vf,x5_ == xf],[t2,t3,tf]);

for k=1:length(solt2)
   tt = [tt;double([(a0-amax)/u_ solt2(k) solt3(k) (af+amax)/u_ soltf(k) u_ 6.1])]; 
end

% u < 0
a5_ = subs(a5,[t1 t4 u],[(a0+amax)/-u_ (af-amax)/-u_ -u_]);
v5_ = subs(v5,[t1 t4 u],[(a0+amax)/-u_ (af-amax)/-u_ -u_]);
x5_ = subs(x5,[t1 t4 u],[(a0+amax)/-u_ (af-amax)/-u_ -u_]);

[solt2,solt3,soltf] = solve([a5_ == af,v5_ == vf,x5_ == xf],[t2,t3,tf]);

for k=1:length(solt2)
   tt = [tt;double([(a0+amax)/-u_ solt2(k) solt3(k) (af-amax)/-u_ soltf(k) -u_ 6.2])]; 
end

%%
a1 = a0 - t1*u;
a2 = a1;
a3 = a2 - (t3-t2)*u;
a4 = a3;
a5 = a4 - (tf-t4)*u;

v1 = v0 + a0*t1      - 1/2*u*t1^2;
v2 = v1 + a1*(t2-t1);
v3 = v2 + a2*(t3-t2) - 1/2*u*(t3-t2)^2;
v4 = v3 + a3*(t4-t3);
v5 = v4 + a4*(tf-t4) - 1/2*u*(tf-t4)^2;

x1 = x0 + v0*t1      + 1/2*a0*t1^2      - 1/6*u*t1^3;
x2 = x1 + v1*(t2-t1) + 1/2*a1*(t2-t1)^2;
x3 = x2 + v2*(t3-t2) + 1/2*a2*(t3-t2)^2 - 1/6*u*(t3-t2)^3;
x4 = x3 + v3*(t4-t3) + 1/2*a3*(t4-t3)^2;
x5 = x4 + v4*(tf-t4) + 1/2*a4*(tf-t4)^2 - 1/6*u*(tf-t4)^3;

%%
% u > 0
a5_ = subs(a5,[t1 t4 u],[(a0-amax)/u_ (-af-amax)/u_ u_]);
v5_ = subs(v5,[t1 t4 u],[(a0-amax)/u_ (-af-amax)/u_ u_]);
x5_ = subs(x5,[t1 t4 u],[(a0-amax)/u_ (-af-amax)/u_ u_]);

[solt2,solt3,soltf] = solve([a5_ == af,v5_ == vf,x5_ == xf],[t2,t3,tf]);

for k=1:length(solt2)
   tt = [tt;double([(a0-amax)/u_ solt2(k) solt3(k) (-af-amax)/u_ soltf(k) u_ 7.1])]; 
end

% u < 0
a5_ = subs(a5,[t1 t4 u],[(a0+amax)/-u_ (-af+amax)/-u_ -u_]);
v5_ = subs(v5,[t1 t4 u],[(a0+amax)/-u_ (-af+amax)/-u_ -u_]);
x5_ = subs(x5,[t1 t4 u],[(a0+amax)/-u_ (-af+amax)/-u_ -u_]);

[solt2,solt3,soltf] = solve([a5_ == af,v5_ == vf,x5_ == xf],[t2,t3,tf]);

for k=1:length(solt2)
   tt = [tt;double([(a0+amax)/-u_ solt2(k) solt3(k) (-af+amax)/-u_ soltf(k) -u_ 7.2])]; 
end

%%
a1 = a0 + t1*u;
a2 = a1;
a3 = a2 - (t3-t2)*u;
a4 = a3;
a5 = a4 - (tf-t4)*u;

v1 = v0 + a0*t1      + 1/2*u*t1^2;
v2 = v1 + a1*(t2-t1);
v3 = v2 + a2*(t3-t2) - 1/2*u*(t3-t2)^2;
v4 = v3 + a3*(t4-t3);
v5 = v4 + a4*(tf-t4) - 1/2*u*(tf-t4)^2;

x1 = x0 + v0*t1      + 1/2*a0*t1^2      + 1/6*u*t1^3;
x2 = x1 + v1*(t2-t1) + 1/2*a1*(t2-t1)^2;
x3 = x2 + v2*(t3-t2) + 1/2*a2*(t3-t2)^2 - 1/6*u*(t3-t2)^3;
x4 = x3 + v3*(t4-t3) + 1/2*a3*(t4-t3)^2;
x5 = x4 + v4*(tf-t4) + 1/2*a4*(tf-t4)^2 - 1/6*u*(tf-t4)^3;

%%
% u > 0
a5_ = subs(a5,[t1 t4 u],[t2 (-af-amax)/u_ u_]);
v5_ = subs(v5,[t1 t4 u],[t2 (-af-amax)/u_ u_]);
x5_ = subs(x5,[t1 t4 u],[t2 (-af-amax)/u_ u_]);

[solt2,solt3,soltf] = solve([a5_ == af,v5_ == vf,x5_ == xf],[t2,t3,tf]);

for k=1:length(solt2)
   tt = [tt;double([solt2(k) solt2(k) solt3(k) (-af-amax)/u_ soltf(k) u_ 8.1])]; 
end

% u < 0
a5_ = subs(a5,[t1 t4 u],[t2 (-af+amax)/-u_ -u_]);
v5_ = subs(v5,[t1 t4 u],[t2 (-af+amax)/-u_ -u_]);
x5_ = subs(x5,[t1 t4 u],[t2 (-af+amax)/-u_ -u_]);

[solt2,solt3,soltf] = solve([a5_ == af,v5_ == vf,x5_ == xf],[t2,t3,tf]);

for k=1:length(solt2)
   tt = [tt;double([solt2(k) solt2(k) solt3(k) (-af+amax)/-u_ soltf(k) -u_ 8.2])]; 
end

%%
% u > 0
a5_ = subs(a5,[t1 t4 u],[(amax-a0)/u_ (-af-amax)/u_ u_]);
v5_ = subs(v5,[t1 t4 u],[(amax-a0)/u_ (-af-amax)/u_ u_]);
x5_ = subs(x5,[t1 t4 u],[(amax-a0)/u_ (-af-amax)/u_ u_]);

[solt2,solt3,soltf] = solve([a5_ == af,v5_ == vf,x5_ == xf],[t2,t3,tf]);

for k=1:length(solt2)    
    tt = [tt;double([(amax-a0)/u_ solt2(k) solt3(k) (-af-amax)/u_ soltf(k) u_ 9.1])];
end

% u < 0
a5_ = subs(a5,[t1 t4 u],[(-amax-a0)/-u_ (-af+amax)/-u_ -u_]);
v5_ = subs(v5,[t1 t4 u],[(-amax-a0)/-u_ (-af+amax)/-u_ -u_]);
x5_ = subs(x5,[t1 t4 u],[(-amax-a0)/-u_ (-af+amax)/-u_ -u_]);

[solt2,solt3,soltf] = solve([a5_ == af,v5_ == vf,x5_ == xf],[t2,t3,tf]);

for k=1:length(solt2)
   tt = [tt;double([(-amax-a0)/-u_ solt2(k) solt3(k) (-af+amax)/-u_ soltf(k) -u_ 9.2])]; 
end

%%
tt = real(tt);
chad = diff(tt')';

idx = [];
for k=1:size(chad,1)
    temp = sum(find(chad(k,1:end-2) >= 0));
    if temp == 10
        idx = [idx;k];
    end
end

if idx ~= 0
    for aa = 1:length(idx)
        ttt = real(tt(idx(aa),1:end-1));
        
        x0_ = x0;
        v0_ = v0;
        a0_ = a0;
        a1_ = subs(a1,[t1 t2 t3 t4 tf u],ttt);
        a2_ = subs(a2,[t1 t2 t3 t4 tf u],ttt);
        a3_ = subs(a3,[t1 t2 t3 t4 tf u],ttt);
        a4_ = subs(a4,[t1 t2 t3 t4 tf u],ttt);
        a5_ = subs(a5,[t1 t2 t3 t4 tf u],ttt);
        
        %%
        tf_ = ttt(end-1);
        tpiece = double([0 ttt(1:end-1)]);
        apiece = double([a1_ a2_ a3_ a4_]);
        temp = find(abs(apiece) > abs(amax)+1e-5);
        temp1 = find(tpiece < 0);
        
        if and(isempty(temp),isempty(temp1))
            tpiece
            apiece
%             sim('simMinTimeTraj');
%             xTraj = traj;
%             figure(2)
%             for k=1:3
%                 subplot(3,1,k)
%                 plot(xTraj.Time,xTraj.Data(:,k))
%                 grid on
%                 box on
%             end
        end
    end
end













