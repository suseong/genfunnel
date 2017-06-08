clear all
close all
clc

checkDependency('yalmip');

e = sdpvar(6,1);
dt = 0.01;

initRegion = diag(1./[0.1 0.1 0.1 0.2 0.2 0.2].^2);
Er = 0.0;
ar = 15;

Kp = diag([10 10 15]);
Kd = diag([4 4 6]);
A = [zeros(3,3) eye(3); -Kp -Kd];    

P = lyap(A',-0.1*eye(6));
P = P/P(1,1);

N = 10;

%%
rho = sdpvar(1,1);
rhodot = sdpvar(1,1);
epbar = sdpvar(1,1);
edbar = sdpvar(1,1);
Ppv = sdpvar(1,1);
Pv = sdpvar(1,1);
x = sdpvar(1,1);

maxKp = max(max(Kp));
maxKd = max(max(Kd));

V = e'*P*e;
Vdot = e'*(P*A+A'*P)*e ...
       + Er*(maxKp*epbar + maxKd*edbar + ar)*(Ppv*epbar + Pv*edbar);

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
    [Lrho,Crho] =       polynomial([e;epbar;edbar ],monomialOrder);
    [Lep,Cep] =         polynomial([e;epbar;edbar ],monomialOrder);
    [Led,Ced] =         polynomial([e;epbar;edbar ],monomialOrder);
    [Lepsign,Cepsign] = polynomial([e;epbar;edbar ],monomialOrder);
    [Ledsign,Cedsign] = polynomial([e;epbar;edbar ],monomialOrder);
    [LPpv1,CPpv1] =     polynomial(x,monomialOrder);
    [LPpv2,CPpv2] =     polynomial(x,monomialOrder);
    [LPv1,CPv1] =       polynomial(x,monomialOrder);
    [LPv2,CPv2] =       polynomial(x,monomialOrder);
    
    vars = [Crho;Cep;Ced;Cepsign;Cedsign];
    vars = [vars;CPpv1;CPpv2;CPv1;CPv2];
    vars = [vars;rhodot;Ppv;Pv];
    
    c1 = sos(rhodot - Vdot ...
             - Lrho*(rhoTemp - V) ...
             - Lep*(epbar^2 - e(1:3)'*e(1:3)) ...
             - Led*(edbar^2 - e(4:6)'*e(4:6)) ...
             - Lepsign*(epbar) ...
             - Ledsign*(edbar) ...
             );
    c2 = sos(Lepsign);
    c3 = sos(Ledsign);

    c4 = sos(Ppv - P(1,4) - LPpv1);
    c5 = sos(Ppv + P(1,4) - LPpv1);
    c6 = sos(Ppv - P(3,6) - LPpv2);
    c7 = sos(Ppv + P(3,6) - LPpv2);
    c8 = sos(LPpv1);
    c9 = sos(LPpv2);
    
    c10 = sos(Pv - P(4,4) - LPv1);
    c11 = sos(Pv + P(4,4) - LPv1);
    c12 = sos(Pv - P(6,6) - LPv2);
    c13 = sos(Pv + P(6,6) - LPv2);
    c14 = sos(LPv1);
    c15 = sos(LPv2);
    
    constraints = [c1 c2 c3 c4 c5 c6 c7 c8 c9 c10 c11 c12 c13 c14 c15];
    sol = solvesos(constraints,rhodot,[],vars);
    
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

%%
[coeffL1,coeffL3,S_] = findL(dt,P_temp,Q_temp,rhoCont,rhodot,Kp,Kd);

%%
[rho,sVars,p,solProblem] = findRho(dt,A,coeffL1,coeffL3,initRegion,Kp,Kd);

%%
ang = -pi:0.2:pi;
for jj = 1:N-2
    figure(101);clf;
    hold on
    P = reshape(double(sVars(:,jj)),6,6);
%     P = reshape(double(S_(:,jj)),6,6);
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

