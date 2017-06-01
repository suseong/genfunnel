clear all
close all

checkDependency('yalmip');

e = sdpvar(6,1);

dt = 0.05;

initRegion = diag(1./[0.5 0.5 0.5 2 2 2].^2);
A = [zeros(3,3) eye(3); -diag([10 10 15]) -diag([4 4 6])];    

P = lyap(A',-0.1*eye(6));
P = P/P(1,1);

%%
V = e'*P*e;
Vdot = e'*(P*A+A'*P)*e;

[L_init,coeff_init] = polynomial(e,2);

rho = sdpvar(1,1);

c1 = sos(L_init);
c2 = sos(rho - V - L_init*(1 - e'*initRegion*e));
sol = solvesos([c1 c2],rho,[],coeff_init);

rhoInit = value(rho);

N = 2;

rho = rhoInit*ones(1,N);
rhodot = diff(rho)/dt;

for k=1:N
    P_{k} = P;
    Q_{k} = P*A +A'*P;
end

[coeffL1,coeffL2,coeffL3] = findL(P_,Q_,rho,rhodot,initRegion);

%%
[rho,sVars,solProblem] = findRho(dt,A,coeffL1,coeffL3,initRegion);

