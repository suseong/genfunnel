clear all
close all
clc

checkDependency('yalmip');

e = sdpvar(6,1);
dt = 0.01;

initRegion = diag(1./[0.1 0.1 0.2 0.1 0.1 0.2].^2);
Er = 0.0;
max_ar = 0;

Kp = diag([10 10 15]);
Kd = diag([4 4 6]);
A = [zeros(3,3) eye(3); -Kp -Kd];    

P = lyap(A',-0.1*eye(6));
P = P/P(1,1);

N = 50;

%%
rho = sdpvar(1,1);
rhodot = sdpvar(1,1);

V = e'*P*e;
Vdot = e'*(P*A+A'*P)*e ;

monomialOrder = 2;
[L_init,coeff_init] = polynomial(e,monomialOrder);

c1 = sos(L_init);
c2 = sos(rho - V ...
         - L_init*(1-e'*initRegion*e));

sol = solvesos([c1 c2],rho,[],coeff_init);
rhoInit = value(rho);

rhoTemp = rhoInit;
rhoCont = rhoInit;

for i=1:N
    [L1,coeff1] = polynomial(e,monomialOrder);
    
    vars = [coeff1];
    
    c1 = sos(rhodot - Vdot ...
             - L1*(rhoTemp - V));
    
    sol = solvesos(c1,rhodot,[],vars);
    
    rhoTemp = rhoTemp + value(rhodot)*dt;
    rhoCont(i+1) = rhoTemp;
end

%%
ts = linspace(0,dt*(N-1),N);
rhodot = diff(rhoCont)/dt;

P_temp = [];
Q_temp = [];

for k = 1:N+1
    P_temp{k} = P;
    Q_temp{k} = P*A + A'*P;
end

[coeffL1,coeffL3,S_] = findL(dt,P_temp,Q_temp,rhoCont,rhodot,Kp,Kd);

rho__ = [];

%%
[rho,sVars,p,solProblem] = findRho(dt,A,coeffL1,coeffL3,initRegion,Kp,Kd,P_temp);

% P_temp = []; Q_temp = [];
% for k = 1:N-2
%    p_temp = p(6*(k-1)+1:6*k);
%    P_temp{k} = [p_temp(1)        0         0  p_temp(3)        0         0;
%                        0  p_temp(1)        0         0  p_temp(3)        0;
%                        0         0  p_temp(2)        0         0  p_temp(4);
%                 p_temp(3)        0         0  p_temp(5)        0         0;
%                        0  p_temp(3)        0         0  p_temp(5)        0;
%                        0         0  p_temp(4)        0         0  p_temp(6)];
%                    
%    Q_temp{k} = P_temp{k}*A +A'*P_temp{k};
% end
% 
% rhodot = diff(rho(1:end-2))/dt;
% 
% [coeffL1,coeffL3,S_] = findL(dt,P_temp,Q_temp,rho(1:end-2),rhodot,Kp,Kd);

%%
ang = -pi:0.2:pi;
for jj = 1:N-2
    figure(101);clf;
    hold on
    P = reshape(double(sVars(:,jj)),6,6);
    kk = 1;
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

