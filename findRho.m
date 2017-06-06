function [rho,sOut,p_,solProblem] = findRho(dt,A,Lc1,Lc3,initRegion,Kp,Kd)

checkDependency('yalmip');
monomialOrder = 2;
Er = 0.1;
max_ar = 5;

N = size(Lc1,2)+1;
rho = sdpvar(N,1);
e = sdpvar(6,1);
epbar = sdpvar(1,1);
edbar = sdpvar(1,1);
epep = e(1:3)'*e(1:3);
eded = e(4:6)'*e(4:6);
maxKp = max(max(Kp));
maxKd = max(max(Kd));

p_ = sdpvar(6*(N+1),1);
p = p_(1:6);
P = [p(1)   0    0  p(3)   0    0;
       0  p(1)   0    0  p(3)   0;
       0    0  p(2)   0    0  p(4);
     p(3)   0    0  p(5)   0    0;
       0  p(3)   0    0  p(5)   0;
       0    0  p(4)   0    0  p(6)];

maxPpv = 0; maxPv = p(5); % need to fix

p = p_(7:12);
Pnext = [p(1)   0    0  p(3)   0    0;
           0  p(1)   0    0  p(3)   0;
           0    0  p(2)   0    0  p(4);
         p(3)   0    0  p(5)   0    0;
           0  p(3)   0    0  p(5)   0;
           0    0  p(4)   0    0  p(6)];

V = e'*P*e;
Vdot = e'*(A'*P + P*A)*e ...
       + Er*(maxPpv*epbar+maxPv*edbar)*(maxKp*epbar+maxKd*edbar +  9.8 + max_ar) ...
       + e'*(Pnext - P)*e/dt;
   
constraints = [];
vars = rho;
cost = 0;
% pVars = p;
% pVars = [];
sVars = [];
sOut = [];

k=1;

S = sdpvar(6,6);

[L1,coeff1] = polynomial([e;epbar;edbar],monomialOrder);
[L2,coeff2] = polynomial([e;epbar;edbar],monomialOrder);
[L3,coeff3] = polynomial([e;epbar;edbar],monomialOrder);

[L4,coeff4] = polynomial([e;epbar;edbar],monomialOrder);
[L5,coeff5] = polynomial([e;epbar;edbar],monomialOrder);
[L6,coeff6] = polynomial([e;epbar;edbar],monomialOrder);
[L7,coeff7] = polynomial([e;epbar;edbar],monomialOrder);

[L8,coeff8] = polynomial([e;epbar;edbar],monomialOrder);
[L9,coeff9] = polynomial([e;epbar;edbar],monomialOrder);
[L10,coeff10] = polynomial([e;epbar;edbar],monomialOrder);
[L11,coeff11] = polynomial([e;epbar;edbar],monomialOrder);

[L12,coeff12] = polynomial([e;epbar;edbar],monomialOrder);
[L13,coeff13] = polynomial([e;epbar;edbar],monomialOrder);
[L14,coeff14] = polynomial([e;epbar;edbar],monomialOrder);
[L15,coeff15] = polynomial([e;epbar;edbar],monomialOrder);

rhodot(k) = (rho(k+1) - rho(k))/dt;

L1 = replace(L1,coeff1,value(Lc1(:,k)));
L3 = replace(L3,coeff3,value(Lc3(:,k)));

c1 = sos(rhodot(k) - Vdot ...
         - L1*(rho(k) - V) ...
         - L4*(epbar^2 - epep) ...
         - L5*(edbar^2 - eded) ...
         - L6*(edbar) ...
         - L7*(epbar));
c2 = sos(rho(k) - V ...
         - L2*(1-e'*initRegion*e) ...
         - L8*(epbar^2 - epep) ...
         - L9*(edbar^2 - eded) ...
         - L10*(edbar) ...
         - L11*(epbar));
c3 = sos(L2);
c4 = sos(1 - e'*S*e ...
         - L3*(rho(k) - V) ...
         - L12*(epbar^2 - epep) ...
         - L13*(edbar^2 - eded) ...
         - L14*(edbar) ...
         - L15*(epbar));
c5 = sos(L6);
c6 = sos(L7);
c8 = sos(L10);
c9 = sos(L11);
c10 = sos(L14);
c11 = sos(L15);
c7 = S >= 0;

constraints = [constraints c1 c2 c3 c4 c5 c6 c7 c8 c9 c10 c11];
cost = cost + geomean(S);

vars = [vars;coeff2;coeff4;coeff5;coeff6;coeff7;coeff8;coeff9;coeff10;coeff11;coeff12;coeff13;coeff14;coeff15];
sVars = [sVars;S(:)];
sOut = [sOut S(:)];

for k=2:N-1
    [L1,coeff1] = polynomial([e;epbar;edbar],monomialOrder);
    [L3,coeff3] = polynomial([e;epbar;edbar],monomialOrder);

    [L4,coeff4] = polynomial([e;epbar;edbar],monomialOrder);
    [L5,coeff5] = polynomial([e;epbar;edbar],monomialOrder);
    [L6,coeff6] = polynomial([e;epbar;edbar],monomialOrder);
    [L7,coeff7] = polynomial([e;epbar;edbar],monomialOrder);
    
    [L8,coeff8] = polynomial([e;epbar;edbar],monomialOrder);
    [L9,coeff9] = polynomial([e;epbar;edbar],monomialOrder);
    [L10,coeff10] = polynomial([e;epbar;edbar],monomialOrder);
    [L11,coeff11] = polynomial([e;epbar;edbar],monomialOrder);
    
    rhodot(k) = (rho(k+1) - rho(k))/dt;

    L1 = replace(L1,coeff1,value(Lc1(:,k)));
    L3 = replace(L3,coeff3,value(Lc3(:,k)));

%     p = sdpvar(5,1);
    p = p_(6*(k-1)+1:6*k);
    P = [p(1)   0    0  p(3)   0    0;
           0  p(1)   0    0  p(3)   0;
           0    0  p(2)   0    0  p(4);
         p(3)   0    0  p(5)   0    0;
           0  p(3)   0    0  p(5)   0;
           0    0  p(4)   0    0  p(6)];
    maxPpv = p(3); maxPv = p(5); % need to fix    

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
           + Er*(maxPpv*epbar+maxPv*edbar)*(maxKp*epbar+maxKd*edbar + 9.8 + max_ar) ...
           + e'*(Pnext - P)*e / dt;
    
    c1 = sos(rhodot(k) - Vdot ...
             - L1*(rho(k) - V) ...
             - L4*(epbar^2 - epep) ...
             - L5*(edbar^2 - eded) ...
             - L6*(edbar) ...
             - L7*(epbar));    
    c4 = sos(1 - e'*S*e ...
             - L3*(rho(k) - V) ... 
             - L8*(epbar^2 - epep) ...
             - L9*(edbar^2 - eded) ...
             - L10*(edbar) ...
             - L11*(epbar));
    c5 = sos(L6);
    c6 = sos(L7);
    c8 = sos(L10);
    c9 = sos(L11);
    c7 = S >= 0;

    constraints = [constraints c1 c4 c5 c6 c7 c8 c9];
    sVars = [sVars;coeff4;coeff5;coeff6;coeff7;coeff8;coeff9;coeff10;coeff11;S(:)];
    sOut = [sOut S(:)];
    
    cost = cost + geomean(S);
end

% sol = solvesos(constraints,-cost,[],[vars;pVars;sVars]);
sol = solvesos(constraints,-cost,[],[vars;p_;sVars]);
solProblem = sol.problem;

rho = value(rho);
sOut = value(sOut);
% p = value(p);
p_ = value(p_);
% p = 0;

end

