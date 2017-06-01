function [rho,sVars,solProblem] = findRho(dt,A,Lc1,Lc3,initRegion)

checkDependency('yalmip');
monomialOrder = 3;
N = size(Lc1,2);
rho = sdpvar(N+1,1);
e = sdpvar(6,1);

constraints = [];
vars = rho;
cost = 0;
pVars = [];
sVars = [];

k=1;

[L1,coeff1] = polynomial(e,monomialOrder);
[L2,coeff2] = polynomial(e,monomialOrder);
[L3,coeff3] = polynomial(e,monomialOrder);

rhodot(k) = (rho(k+1) - rho(k))/dt;

L1 = replace(L1,coeff1,value(Lc1(:,k)));
% L2 = replace(L2,coeff2,value(Lc2));
L3 = replace(L3,coeff3,value(Lc3(:,k)));

% p = sdpvar(6,1);
% P = [ p(1)   0    0  p(3)   0    0;
%         0  p(1)   0    0  p(3)   0;
%         0    0  p(2)   0    0  p(4);
%       p(3)   0    0  p(5)   0    0;
%         0  p(3)   0    0  p(5)   0;
%         0    0  p(4)   0    0  p(6)];

p = sdpvar(5,1);
P = [   1    0    0  p(2)   0    0;
        0    1    0    0  p(2)   0;
        0    0  p(1)   0    0  p(3);
      p(2)   0    0  p(4)   0    0;
        0  p(2)   0    0  p(4)   0;
        0    0  p(3)   0    0  p(5)];

V = e'*P*e;
Vdot = e'*(P*A+A'*P)*e;

S = sdpvar(6,6);
    
vars = [vars;coeff2];
pVars = [pVars;p];
sVars = [sVars;S(:)];

c1 = sos(rhodot(k) - Vdot - L1*(rho(k) - V));
c2 = sos(rho(k) - V - L2*(1-e'*initRegion*e));
c3 = sos(L2);
c4 = sos(1 - e'*S*e - L3*(rho(k) - V));
c5 = S > 0;
constraints = [constraints c1 c2 c3 c4 c5];
cost = cost + geomean(S);

for k=2:N
    [L1,coeff1] = polynomial(e,monomialOrder);
    [L3,coeff3] = polynomial(e,monomialOrder);

    rhodot(k) = (rho(k+1) - rho(k))/dt;

    L1 = replace(L1,coeff1,Lc1(:,k));
    L3 = replace(L3,coeff3,value(Lc3(:,k)));
        
    p = sdpvar(5,1);
    P = [   1    0    0  p(2)   0    0;
            0    1    0    0  p(2)   0;
            0    0  p(1)   0    0  p(3);
          p(2)   0    0  p(4)   0    0;
            0  p(2)   0    0  p(4)   0;
            0    0  p(3)   0    0  p(5)];
    
%     p = sdpvar(6,1);
%     P = [ p(1)   0    0  p(3)   0    0;
%             0  p(1)   0    0  p(3)   0;
%             0    0  p(2)   0    0  p(4);
%           p(3)   0    0  p(5)   0    0;
%             0  p(3)   0    0  p(5)   0;
%             0    0  p(4)   0    0  p(6)];

    V = e'*P*e;
    Vdot = e'*(P*A+A'*P)*e;

    S = sdpvar(6,6);
        
    pVars = [pVars;p];
    sVars = [sVars;S(:)];

    c1 = sos(rhodot(k) - Vdot - L1*(rho(k) - V));
%     c3 = sos(L3);
    c4 = sos(1 - e'*S*e - L3*(rho(k) - V));
    c5 = S > 0;
    constraints = [constraints c1 c4 c5];
    cost = cost + geomean(S);
end

sol = solvesos(constraints,-cost,[],[vars;pVars;sVars]);
solProblem = sol.problem;

end


% S = [ s(1)   0    0  s(3)   0    0;
%         0  s(1)   0    0  s(3)   0;
%         0    0  s(2)   0    0  s(4);
%       s(3)   0    0  s(5)   0    0;
%         0  s(3)   0    0  s(5)   0;
%         0    0  s(4)   0    0  s(6)];

% p = sdpvar(5,1);
% P = [   1    0    0  p(2)   0    0;
%         0    1    0    0  p(2)   0;
%         0    0  p(1)   0    0  p(3);
%       p(2)   0    0  p(4)   0    0;
%         0  p(2)   0    0  p(4)   0;
%         0    0  p(3)   0    0  p(5)];


