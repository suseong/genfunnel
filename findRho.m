function [rho,sOut,p_,solProblem] = findRho(dt,A,Lc1,Lc3,LcPe,initRegion,Kp,Kd)

checkDependency('yalmip');
monomialOrder = 2;
Er = 0.1;
max_ar = 0;

N = size(Lc1,2)-1;
rho = sdpvar(N+1,1);
e = sdpvar(6,1);
% Pe = sdpvar(1,1);
Ke = sdpvar(1,1);

p_ = sdpvar(6*(N+1),1);
p = p_(1:6);
P = [p(1)   0    0  p(3)   0    0;
       0  p(1)   0    0  p(3)   0;
       0    0  p(2)   0    0  p(4);
     p(3)   0    0  p(5)   0    0;
       0  p(3)   0    0  p(5)   0;
       0    0  p(4)   0    0  p(6)];

k=1;

p = p_(7:12);
Pnext = [p(1)   0    0  p(3)   0    0;
           0  p(1)   0    0  p(3)   0;
           0    0  p(2)   0    0  p(4);
         p(3)   0    0  p(5)   0    0;
           0  p(3)   0    0  p(5)   0;
           0    0  p(4)   0    0  p(6)];

V = e'*P*e;
Vdot = e'*(A'*P + P*A)*e ...
       + Er*1*(Ke + max_ar) ...
       + e'*(Pnext - P)*e/dt;
   
constraints = [];
vars = rho;
cost = 0;
sVars = [];
sOut = [];

S = sdpvar(6,6);

[L1,coeff1] = polynomial([e;Ke],monomialOrder);
[L2,coeff2] = polynomial([e;Ke],monomialOrder);
[L3,coeff3] = polynomial([e;Ke],monomialOrder);
% [LPe,coeffPe] = polynomial([e;Pe;Ke],monomialOrder);
[LKe,coeffKe] = polynomial([e;Ke],monomialOrder);

rhodot(k) = (rho(k+1) - rho(k))/dt;

L1 = replace(L1,coeff1,value(Lc1(:,k)));
L3 = replace(L3,coeff3,value(Lc3(:,k)));
% LPe = replace(LPe,coeffPe,value(LcPe(:,k)));

c1 = sos(rhodot(k) - Vdot ...
         - L1*(rho(k) - V) ...
         - LKe*(Ke - Kp*e(1:3) - Kd*e(4:6)));
c2 = sos(rho(k) - V ...
         - L2*(1-e'*initRegion*e));
c3 = sos(L2);
c4 = sos(1 - e'*S*e ...
         - L3*(rho(k) - V));
c5 = S >= 0;

constraints = [constraints c1 c2 c3 c4 c5];
cost = cost + geomean(S);

vars = [vars;coeff2;coeffKe]; 
sVars = [sVars;S(:)];
sOut = [sOut S(:)];

for k=2:N-1
    [L1,coeff1] = polynomial([e;Ke],monomialOrder);
    [L3,coeff3] = polynomial([e;Ke],monomialOrder);
%     [LPe,coeffPe] = polynomial([e;Pe;Ke],monomialOrder);
    [LKe,coeffKe] = polynomial([e;Ke],monomialOrder);

    rhodot(k) = (rho(k+1) - rho(k))/dt;

    L1 = replace(L1,coeff1,value(Lc1(:,k)));
    L3 = replace(L3,coeff3,value(Lc3(:,k)));
%     LPe = replace(LPe,coeffPe,value(LcPe(:,k)));

    p = p_(6*(k-1)+1:6*k);
    P = [p(1)   0    0  p(3)   0    0;
           0  p(1)   0    0  p(3)   0;
           0    0  p(2)   0    0  p(4);
         p(3)   0    0  p(5)   0    0;
           0  p(3)   0    0  p(5)   0;
           0    0  p(4)   0    0  p(6)];

    p = p_(6*k+1:6*(k+1));
    Pnext = [p(1)   0    0  p(3)   0    0;
               0  p(1)   0    0  p(3)   0;
               0    0  p(2)   0    0  p(4);
             p(3)   0    0  p(5)   0    0;
               0  p(3)   0    0  p(5)   0;
               0    0  p(4)   0    0  p(6)];
    S = sdpvar(6,6);
    V = e'*P*e;
    Vdot = e'*(A'*P + P*A)*e ...
           + Er*1*(Ke + max_ar) ...
           + e'*(Pnext - P)*e / dt;
    
    c1 = sos(rhodot(k) - Vdot ...
             - L1*(rho(k) - V) ...
             - LKe*(Ke - Kp*e(1:3) - Kd*e(4:6)));
    c4 = sos(1 - e'*S*e ...
             - L3*(rho(k) - V));
    c7 = S >= 0;

    constraints = [constraints c1 c4 c7];
    vars = [vars;coeffKe];
    sVars = [sVars;S(:)];
    sOut = [sOut S(:)];
    
    cost = cost + geomean(S);
end

sol = solvesos(constraints,-cost,[],[vars;p_;sVars]);
solProblem = sol.problem;

rho = value(rho);
sOut = value(sOut);
p_ = value(p_);

end

%          - LPe*(Pe - P(1:3,4:6)*e(1:3) - P(4:6,4:6)*e(4:6)) ...
%              - LPe*(Pe - P(1:3,4:6)*e(1:3) - P(4:6,4:6)*e(4:6)) ...
