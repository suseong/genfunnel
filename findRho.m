function [rho,sOut,p,solProblem] = findRho(dt,rhoInit,P,A,Lc1,Lc3,initRegion,Kp,Kd)

checkDependency('yalmip');
monomialOrder = 2;
Er = 0.08;

N = size(Lc1,2);
rho = sdpvar(N+1,1);
e = sdpvar(6,1);
epbar = sdpvar(1,1);
edbar = sdpvar(1,1);
epep = e(1:3)'*e(1:3);
eded = e(4:6)'*e(4:6);
assign(rho(1),rhoInit);

% P = sdpvar(6,6);
p = sdpvar(5,1);
P = [   1    0    0  p(2)   0    0;
        0    1    0    0  p(2)   0;
        0    0  p(1)   0    0  p(3);
      p(2)   0    0  p(4)   0    0;
        0  p(2)   0    0  p(4)   0;
        0    0  p(3)   0    0  p(5)];

maxKp = max(max(Kp));
maxKd = max(max(Kd));
% maxPpv = max(max(P(1:3,4:6)));
% maxPv = max(max(P(4:6,4:6)));
maxPpv = p(2); maxPv = p(4); % need to fix
% maxPpv = 0.04; maxPv = 0.1;
    
V = e'*P*e;
Vdot = e'*(A'*P + P*A)*e ...
       + Er*(maxPpv*epbar+maxPv*edbar)*(maxKp*epbar+maxKd*edbar + 18);
   
constraints = [];
vars = rho;
cost = 0;
pVars = p;
% pVars = [];
sVars = [];
sOut = [];

k=1;

S = sdpvar(6,6);

[L1,coeff1] = polynomial([e;epbar;edbar],monomialOrder);
[L2,coeff2] = polynomial(e,monomialOrder);
[L3,coeff3] = polynomial(e,monomialOrder);

[L4,coeff4] = polynomial([e;epbar;edbar],monomialOrder);
[L5,coeff5] = polynomial([e;epbar;edbar],monomialOrder);
[L6,coeff6] = polynomial([e;epbar;edbar],monomialOrder);
[L7,coeff7] = polynomial([e;epbar;edbar],monomialOrder);

rhodot(k) = (rho(k+1) - rho(k))/dt;

L1 = replace(L1,coeff1,value(Lc1(:,k)));
L3 = replace(L3,coeff3,value(Lc3(:,k)));

c1 = sos(rhodot(k) - Vdot ...
         - L1*(rho(k) - V) ...
         - L4*(epbar^2 - epep) ...
         - L5*(edbar^2 - eded) ...
         - L6*(edbar) ...
         - L7*(epbar));
c2 = sos(rho(k) - V - L2*(1-e'*initRegion*e));
c3 = sos(L2);
c4 = sos(1 - e'*S*e - L3*(rho(k) - V));
c5 = sos(L6);
c6 = sos(L7);
c7 = S >= 0;
constraints = [constraints c1 c2 c3 c4 c5 c6 c7];
cost = cost + geomean(S);

vars = [vars;coeff2;coeff4;coeff5;coeff6;coeff7];
sVars = [sVars;S(:)];
sOut = [sOut S(:)];

for k=2:N
    [L1,coeff1] = polynomial([e;epbar;edbar],monomialOrder);
    [L3,coeff3] = polynomial(e,monomialOrder);

    [L4,coeff4] = polynomial([e;epbar;edbar],monomialOrder);
    [L5,coeff5] = polynomial([e;epbar;edbar],monomialOrder);
    [L6,coeff6] = polynomial([e;epbar;edbar],monomialOrder);
    [L7,coeff7] = polynomial([e;epbar;edbar],monomialOrder);

    rhodot(k) = (rho(k+1) - rho(k))/dt;

    L1 = replace(L1,coeff1,value(Lc1(:,k)));
    L3 = replace(L3,coeff3,value(Lc3(:,k)));
        
    S = sdpvar(6,6);
    sVars = [sVars;coeff4;coeff5;coeff6;coeff7;S(:)];
    sOut = [sOut S(:)];
    
    V = e'*P*e;
    Vdot = e'*(A'*P + P*A)*e ...
           + Er*(maxPpv*epbar+maxPv*edbar)*(maxKp*epbar+maxKd*edbar + 18);
    
    c1 = sos(rhodot(k) - Vdot ...
             - L1*(rho(k) - V) ...
             - L4*(epbar^2 - epep) ...
             - L5*(edbar^2 - eded) ...
             - L6*(edbar) ...
             - L7*(epbar));    
    c4 = sos(1 - e'*S*e - L3*(rho(k) - V)); 
    c5 = sos(L6);
    c6 = sos(L7);
    c7 = S >= 0;

    constraints = [constraints c1 c4 c5 c6 c7];
    cost = cost + geomean(S);
end

sol = solvesos(constraints,-cost,[],[vars;pVars;sVars]);
solProblem = sol.problem;

rho = value(rho);
sOut = value(sOut);
p = value(p);
% p = 0;

end

