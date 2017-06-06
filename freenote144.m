clear all
close all
clc

checkDependency('yalmip');

e = sdpvar(6,1);
dt = 0.005;

initRegion = diag(1./[0.2 0.2 0.2 0.2 0.2 0.2].^2);
max_ar = 0;

Kp = diag([10 10 15]);
Kd = diag([4 4 6]);
A = [zeros(3,3) eye(3); -Kp -Kd];    

P = lyap(A',-0.1*eye(6));
P = P/P(1,1);

%%
maxPpv = max(max(P(1:3,4:6)));
maxPp = max(max(P(1:3,1:3)));
maxPv = max(max(P(4:6,4:6)));
maxKp = max(max(Kp));
maxKd = max(max(Kd));

Er = 0.1;

rho = sdpvar(1,1);
rhodot = sdpvar(1,1);
epbar = sdpvar(1,1);
edbar = sdpvar(1,1);

epep = e(1:3)'*e(1:3);
eded = e(4:6)'*e(4:6);
V = e'*P*e;

Vdot = e'*(P*A+A'*P)*e ...
       + Er*(maxPpv*epbar+maxPv*edbar)*(maxKp*epbar+maxKd*edbar + 9.8 + max_ar);

monomialOrder = 2;
[L_init,coeff_init] = polynomial(e,monomialOrder);

c1 = sos(L_init);
c2 = sos(rho - V - L_init*(1-e'*initRegion*e));

sol = solvesos([c1 c2],rho,[],coeff_init);
rhoInit = value(rho);

N = 20;
rhoTemp = rhoInit;
rhoCont = rhoInit;

for i=1:N
    [L1,coeff1] = polynomial([e;epbar;edbar],monomialOrder);
    [L2,coeff2] = polynomial([e;epbar;edbar],monomialOrder);
    [L3,coeff3] = polynomial([e;epbar;edbar],monomialOrder);
    [L4,coeff4] = polynomial([e;epbar;edbar],monomialOrder);
    [L5,coeff5] = polynomial([e;epbar;edbar],monomialOrder);
    
    vars = [coeff1;coeff2;coeff3;coeff4;coeff5];
    
    c1 = sos(rhodot - Vdot ...
             - L1*(rhoTemp - V) ...
             - L2*(epbar^2 - epep) ...
             - L3*(edbar^2 - eded) ...
             - L4*(edbar) ...
             - L5*(epbar));
    c2 = sos(L4);
    c3 = sos(L5);
    
    sol = solvesos([c1 c2 c3],rhodot,[],vars);
    
    rhoTemp = rhoTemp + value(rhodot)*dt;
    rhoCont(i+1) = rhoTemp;
end

%%
ts = linspace(0,dt*(N-1),N);
rhodot = diff(rhoCont)/dt;

[coeffL1,coeffL3,S_] = findL(P,P*A +A'*P,rhoCont,rhodot,initRegion,Kp,Kd);

%%
[rho,sVars,p,solProblem] = findRho(dt,rhoInit,A,coeffL1,coeffL3,initRegion,Kp,Kd);

%%
% Kp = diag([10 10 15]);
% Kd = diag([4 4 6]);
% A = [zeros(3,3) eye(3); -Kp -Kd];    

ang = -pi:0.2:pi;
for jj = 1:N-1
    figure(101);clf;
    hold on
    P = reshape(double(sVars(:,jj)),6,6);
    kk = 3;
    p1 = [P(kk,kk) P(kk,kk+3);P(kk+3,kk) P(kk+3,kk+3)];
    invp1 = inv(sqrtm(p1));
    p2 = [initRegion(kk,kk) initRegion(kk,kk+3);initRegion(kk+3,kk) initRegion(kk+3,kk+3)];
    invp2 = inv(sqrtm(p2));
for k=1:length(ang)
   xx = invp1*[cos(ang(k));sin(ang(k))];
   plot(xx(1),xx(2),'.') 
   yy = invp2*[cos(ang(k));sin(ang(k))];
   my_phase(yy,jj);
end
   axis([-0.5 0.5 -0.5 0.5]);
   axis equal
   pause(0.1);
end









