clear all
close all
clc

checkDependency('yalmip');

e = sdpvar(6,1);
dt = 0.03;

initRegion = diag(1./[0.5 0.5 0.5 1.6 1.6 1.6].^2);
% initRegion = diag(1./[1 1 1 1 1 1].^2);
Kp = diag([10 10 15]);
Kd = diag([4 4 6]);
A = [zeros(3,3) eye(3); -Kp -Kd];    

P = lyap(A',-0.1*eye(6));
P = P/P(1,1);

%%
maxPpv = max(max(P(1:3,4:6)));
maxPp = max(max(P(1:3,1:3)));
maxPv = max(max(P(4:6,4:6)));

rho = sdpvar(1,1);
y1 = sdpvar(1,1);

V = e'*P*e;
Vdot = e'*(P*A+A'*P)*e + y1;

monimialOrder = 3;
[L_init,coeff_init] = polynomial(e,monimialOrder);
[L__,coeff__] = polynomial(e,monimialOrder);

c1 = sos(L_init);
c2 = sos(rho - V - L_init*(1 - e'*initRegion*e));

sol = solvesos([c1 c2],rho,[],[coeff_init]);
rhoInit = value(rho);

ll = sdpvar(1,1);
N = 100;
rhoTemp = rhoInit;
rhoCont = [rhoInit];
KK(1:3,1:3) = Kp^2; KK(4:6,4:6) = Kd^2;

for i=1:N-1
    P__ = P(4:6,1:6);
    eKKe = e'*KK*e;
    eP2e = e'*P__'*P__*e;
    c_ = sos(ll-eP2e + L__*(rhoTemp - e'*P*e));
    sol = solvesos(c_,ll,[],coeff__);
    Pe_norm = sqrt(value(ll));
    c_ = sos(ll-eKKe + L__*(rhoTemp - e'*P*e));    
    sol = solvesos(c_,ll,[],coeff__);
    Ke_norm = sqrt(value(ll));

    epbar = sqrt(rhoTemp/P(3,3));
    evbar = sqrt(rhoTemp/P(6,6));
    Er = 0.01;
%     maxD = Er*(maxPpv*epbar+maxPv*evbar)*(15*epbar+6*evbar+15);
    maxD = Er*Pe_norm*(Ke_norm+15);
    
    rhodot = sdpvar(1,1);
    [L1,coeff1] = polynomial(e,monimialOrder);
    [L_chad,coeff_chad] = polynomial(e,monimialOrder);
    c3 = sos(rhodot - Vdot - L1*(rhoTemp - V) - L_chad*(maxD - y1));
    c4 = sos(L_chad);
    sol = solvesos([c3 c4],rhodot,[],[coeff1;coeff_chad]);
    rhoTemp = rhoTemp + value(rhodot)*dt;
    rhoCont(i+1) = rhoTemp;
end

%%
ts = linspace(0,dt*(N-1),N);
rhodot = diff(rhoCont)/dt;
% rho = rho_first*exp(-0.0*(ts-ts(1))/(ts(end)-ts(1))); 
% rho = rhoInit*ones(1,N);
% rho = rhoInit*exp(-0.1*(ts-ts(1))/ts(end)-ts(1));
% rhodot = diff(rho)/dt;

% for k=1:N
%     P_{k} = P;
%     Q_{k} = P*A +A'*P;
% end

[coeffL1,coeffL3,S_] = findL(P,P*A +A'*P,rhoCont,rhodot,initRegion);

%%
[rho,sVars,p,solProblem] = findRho(dt,A,coeffL1,coeffL3,initRegion);

%%
figure(100);
P = reshape(double(sVars(:,end)),6,6);
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
