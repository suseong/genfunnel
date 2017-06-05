clear all
close all
clc

checkDependency('yalmip');

e = sdpvar(6,1);
dt = 0.03;

initRegion = diag(1./[0.5 0.5 0.5 1.6 1.6 1.6].^2);

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

Er = 0.08;

rho = sdpvar(1,1);
rhodot = sdpvar(1,1);
epbar = sdpvar(1,1);
edbar = sdpvar(1,1);

epep = e(1:3)'*e(1:3);
eded = e(4:6)'*e(4:6);
V = e'*P*e;

Vdot = e'*(P*A+A'*P)*e ...
       + Er*(maxPpv*epbar+maxPv*edbar)*(maxKp*epbar+maxKd*edbar + 18);

monomialOrder = 2;
[L_init,coeff_init] = polynomial(e,monomialOrder);

c1 = sos(L_init);
c2 = sos(rho - V - L_init*(1-e'*initRegion*e));

sol = solvesos([c1 c2],rho,[],coeff_init);
rhoInit = value(rho);

N = 200;
rhoTemp = rhoInit;
rhoCont = rhoInit;

for i=1:N
%     [L1,coeff1] = polynomial(e,monomialOrder);
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
for kk = 1:3
[rho,sVars,p,solProblem] = findRho(dt,rhoInit,P,A,coeffL1,coeffL3,initRegion,Kp,Kd);

%%
ts = linspace(0,dt*(N-1),N);
rhodot = diff(rho)/dt;

% P = [   1    0    0  p(2)   0    0;
%         0    1    0    0  p(2)   0;
%         0    0  p(1)   0    0  p(3);
%       p(2)   0    0  p(4)   0    0;
%         0  p(2)   0    0  p(4)   0;
%         0    0  p(3)   0    0  p(5)];
    
[coeffL1,coeffL3,S_] = findL(P,P*A +A'*P,rho,rhodot,initRegion,Kp,Kd);

end
%%
figure(100);
P = reshape(double(sVars(:,end)),6,6);
% P = SS;
for kk=1:3
subplot(1,3,kk)
% p1 = [P(kk,kk) P(kk,kk+3);P(kk+3,kk) P(kk+3,kk+3)]/(value(rho));
p1 = [P(kk,kk) P(kk,kk+3);P(kk+3,kk) P(kk+3,kk+3)];
p2 = [initRegion(kk,kk) initRegion(kk,kk+3);initRegion(kk+3,kk) initRegion(kk+3,kk+3)];
invp1 = inv(sqrtm(p1));
invp2 = inv(sqrtm(p2));
ang = -pi:0.02:pi;
hold on
for k=1:length(ang)
   xx = invp1*[cos(ang(k));sin(ang(k))];
   plot(xx(1),xx(2),'*') 
   yy = invp2*[cos(ang(k));sin(ang(k))];
   plot(yy(1),yy(2),'.')    
end
axis equal
end
