function [rho,sOut,p_,solProblem] = findRho(dt,A,Crho1_,Crho3_,initRegion,Kp,Kd,Er,ar)

checkDependency('yalmip');
monomialOrder = 2;
% Er = 0.0;
% ar = 15;
maxKp = max(max(Kp));
maxKd = max(max(Kd));

N = size(Crho1_,2);

rho = sdpvar(N+1,1);
p_ = sdpvar(6*(N+1),1);

e = sdpvar(6,1);
epbar = sdpvar(1,1);
edbar = sdpvar(1,1);

constraints = [];
cost = 0;
sVars = [];
sOut = [];

%%
k=1;

Ppv = sdpvar(1,1);
Pv = sdpvar(1,1);
S = sdpvar(6,6);

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

V = e'*P*e;
Vdot = e'*(P*A+A'*P)*e ...
       + Er*(maxKp*epbar + maxKd*edbar + ar)*(Ppv*epbar + Pv*edbar) ...
       + e'*(Pnext - P)*e / dt;
rhodot = (rho(k+1) - rho(k))/dt;

[Lrho1,Crho1] =     polynomial([e;epbar;edbar],monomialOrder);
[Lrho2,Crho2] =     polynomial(e,monomialOrder+1);       
[Lrho3,Crho3] =     polynomial(e,monomialOrder);       
[Lep,Cep] =         polynomial([e;epbar;edbar],monomialOrder);
[Led,Ced] =         polynomial([e;epbar;edbar],monomialOrder);
[Lepsign,Cepsign] = polynomial([e;epbar;edbar],monomialOrder);
[Ledsign,Cedsign] = polynomial([e;epbar;edbar],monomialOrder);

Lrho1 = replace(Lrho1,Crho1,value(Crho1_(:,k)));
Lrho3 = replace(Lrho3,Crho3,value(Crho3_(:,k)));

c1 = sos(rhodot - Vdot ...
         - Lrho1*(rho(k) - V) ...
         - Lep*(epbar^2 - e(1:3)'*e(1:3)) ...
         - Led*(edbar^2 - e(4:6)'*e(4:6)) ...
         - Lepsign*(epbar) ...
         - Ledsign*(edbar) ...
         );
c2 = sos(Lepsign);
c3 = sos(Ledsign);

c4 = sos(rho(k) - V ...
         -Lrho2*(1-e'*initRegion*e));
c5 = sos(Lrho2);

c6 = sos(1 - e'*S*e ...
         -Lrho3*(rho(k) - V));
     
c7 = S >= 0;

c8 = sos(Ppv - P(1,4));
c9 = sos(Ppv + P(1,4));
c10 = sos(Ppv - P(3,6));
c11 = sos(Ppv + P(3,6));

c12 = sos(Pv - P(4,4));
c13 = sos(Pv + P(4,4));
c14 = sos(Pv - P(6,6));
c15 = sos(Pv + P(6,6));
   
constraints = [constraints c1 c2 c3 c4 c5 c6 c7 ...
               c8 c9 c10 c11 c12 c13 c14 c15];

vars = [Ppv;Pv;Crho2;Cep;Ced;Cepsign;Cedsign]; 
sVars = [sVars;S(:)];
sOut = [sOut S(:)];

cost = cost + geomean(S);

for k=2:N
    Ppv = sdpvar(1,1);
    Pv = sdpvar(1,1);
    S = sdpvar(6,6);

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

    V = e'*P*e;
    Vdot = e'*(P*A+A'*P)*e ...
           + Er*(maxKp*epbar + maxKd*edbar + ar)*(Ppv*epbar + Pv*edbar) ...
           + e'*(Pnext - P)*e / dt;
    rhodot = (rho(k+1) - rho(k))/dt;
   
    [Lrho1,Crho1] =     polynomial([e;epbar;edbar],monomialOrder);
    [Lrho3,Crho3] =     polynomial(e,monomialOrder);       
    [Lep,Cep] =         polynomial([e;epbar;edbar],monomialOrder);
    [Led,Ced] =         polynomial([e;epbar;edbar],monomialOrder);
    [Lepsign,Cepsign] = polynomial([e;epbar;edbar],monomialOrder);
    [Ledsign,Cedsign] = polynomial([e;epbar;edbar],monomialOrder);

    Lrho1 = replace(Lrho1,Crho1,value(Crho1_(:,k)));
    Lrho3 = replace(Lrho3,Crho3,value(Crho3_(:,k)));

    c1 = sos(rhodot - Vdot ...
             - Lrho1*(rho(k) - V) ...
             - Lep*(epbar^2 - e(1:3)'*e(1:3)) ...
             - Led*(edbar^2 - e(4:6)'*e(4:6)) ...
             - Lepsign*(epbar) ...
             - Ledsign*(edbar) ...
             );
    c2 = sos(Lepsign);
    c3 = sos(Ledsign);

    c6 = sos(1 - e'*S*e ...
             -Lrho3*(rho(k) - V));

    c7 = S >= 0;

    c8 = sos(Ppv - P(1,4));
    c9 = sos(Ppv + P(1,4));
    c10 = sos(Ppv - P(3,6));
    c11 = sos(Ppv + P(3,6));

    c12 = sos(Pv - P(4,4));
    c13 = sos(Pv + P(4,4));
    c14 = sos(Pv - P(6,6));
    c15 = sos(Pv + P(6,6));
   
    constraints = [constraints c1 c2 c3 c6 c7 ...
                   c8 c9 c10 c11 c12 c13 c14 c15];

    vars = [vars;Ppv;Pv;Crho2;Cep;Ced;Cepsign;Cedsign]; 
    sVars = [sVars;S(:)];
    sOut = [sOut S(:)];

    cost = cost + geomean(S);
end

sol = solvesos(constraints,-cost,[],[vars;sVars;rho;p_]);
solProblem = sol.problem;

rho = value(rho);
sOut = value(sOut);
p_ = value(p_);

end
