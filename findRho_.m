function [rho,rhodot,solProblem] = findRho_(x,V,Vdot,L1_c,L2_c,dts,initRegion)

checkDependency('yalmip');
N = length(V)-1;
rho = sdpvar(N+1,1);
% rho(1) = rho1;

constraints = [];
cost = 0;
vars = [];
xs = sdpvar(9,1);

i=1;
[L1,coeff1_] = polynomial(xs,2);
[L2,coeff2_] = polynomial(xs,2);
rhodot(i) = (rho(i+1) - rho(i))/dts(i);

L1 = replace(L1,coeff1_,L1_c(:,i));
L2 = replace(L2,coeff2_,L2_c);

Vdots = msspoly2sdpvar(x,xs,Vdot{i});
Vs = msspoly2sdpvar(x,xs,V{i});

c1 = sos(rhodot(1) - Vdots - L1*(rho(i) - Vs));
c2 = sos(rho(1) - Vs - L2*(1-xs'*initRegion*xs));
constraints = [constraints c1 c2];

D{i} = [];
cost = cost + rho(1);

for i=2:N
    [L1,coeff1_] = polynomial(xs,2);
    rhodot(i) = (rho(i+1) - rho(i))/dts(i);

    L1 = replace(L1,coeff1_,L1_c(:,i));

    Vdots = msspoly2sdpvar(x,xs,Vdot{i});
    Vs = msspoly2sdpvar(x,xs,V{i});
    
    c1 = sos(rhodot(i) - Vdots - L1*(rho(i) - Vs));
    constraints = [constraints c1];
    
    D{i} = [];
    cost = cost + rho(i);
end

cost = cost + rho(i+1);
vars = rho;
sol = solvesos(constraints,cost,[],vars);
solProblem = sol.problem;

end