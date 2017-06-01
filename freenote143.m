clear all
close all

checkDependency('yalmip');

e = sdpvar(6,1);

dt = 0.01;

initRegion = diag(1./[0.5 0.5 0.5 2 2 2].^2);
% initRegion = diag(1./[1 1 1 1 1 1].^2);
A = [zeros(3,3) eye(3); -diag([10 10 15]) -diag([4 4 6])];    

P = lyap(A',-0.1*eye(6));
P = P/P(1,1);

%%
V = e'*P*e;
Vdot = e'*(P*A+A'*P)*e;

monimialOrder = 2;
[L_init,coeff_init] = polynomial(e,monimialOrder);

rho = sdpvar(1,1);

c1 = sos(L_init);
c2 = sos(rho - V - L_init*(1 - e'*initRegion*e));
sol = solvesos([c1 c2],rho,[],coeff_init);
rhoInit = value(rho);

N = 300;
rhoTemp = rhoInit;
rhoCont = [rhoInit];
for i=1:N-1
    rhodot = sdpvar(1,1);
    [L1,coeff1] = polynomial(e,monimialOrder);
    c3 = sos(rhodot - Vdot - L1*(rhoTemp - V));
    sol = solvesos(c3,rhodot,[],[coeff_init;coeff1]);
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

for k=1:N
    P_{k} = P;
    Q_{k} = P*A +A'*P;
end

[coeffL1,coeffL2,coeffL3,S_] = findL(P_,Q_,rhoCont,rhodot,initRegion);

%%
[rho,sVars,solProblem] = findRho(dt,A,coeffL1,coeffL3,initRegion);

%%
figure(100);
P = reshape(double(sVars(:,1)),6,6);
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
