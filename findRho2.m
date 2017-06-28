function [rho,sOut,p_,solProblem] = findRho2(dt,A,Crho1_,Crho3_,initRegion,unc)

monomialOrder = 2;

N = size(Crho1_,2);

rho = sdpvar(N,1);
p_ = sdpvar(3*N,1);

e = sdpvar(2,1);
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

s = sdpvar(3,1);
S = [s(1) s(2);
     s(2) s(3)];

p = p_(1:3);
P = [p(1) p(2);
     p(2) p(3)];

p = p_(3*k+1:3*(k+1));
Pnext = [p(1) p(2);
         p(2) p(3)];

V = e'*P*e;
Vdot = e'*(P*A+A'*P)*e ...
       + 2*unc*(Ppv*epbar + Pv*edbar) ...
       + e'*(Pnext - P)*e / dt;
rhodot = (rho(k+1) - rho(k))/dt;

[Lrho1,Crho1] =     polynomial([e;epbar;edbar],monomialOrder);
[Lrho2,Crho2] =     polynomial(e,monomialOrder);       
[Lrho3,Crho3] =     polynomial(e,monomialOrder);       
[Lep,Cep] =         polynomial([e;epbar;edbar],monomialOrder);
[Led,Ced] =         polynomial([e;epbar;edbar],monomialOrder);
[Lepsign,Cepsign] = polynomial([e;epbar;edbar],monomialOrder);
[Ledsign,Cedsign] = polynomial([e;epbar;edbar],monomialOrder);

Lrho1 = replace(Lrho1,Crho1,value(Crho1_(:,k)));
Lrho3 = replace(Lrho3,Crho3,value(Crho3_(:,k)));

c1 = sos(rhodot - Vdot ...
         - Lrho1*(rho(k) - V) ...
         - Lep*(epbar^2 - e(1)*e(1)) ...
         - Led*(edbar^2 - e(2)*e(2)) ...
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

c8 = sos(Ppv - P(1,2));
c9 = sos(Ppv + P(1,2));

c10 = sos(Pv - P(2,2));
c11 = sos(Pv + P(2,2));
   
constraints = [constraints c1 c2 c3 c4 c5 c6 c7 c8 c9 c10 c11];

vars = [Ppv;Pv;Crho2;Cep;Ced;Cepsign;Cedsign]; 
sVars = [sVars;s];
sOut = [sOut s];

cost = cost + geomean(S);

for k=2:N
%     Ppv = sdpvar(1,1);
%     Pv = sdpvar(1,1);
    s = sdpvar(3,1);
    S = [s(1) s(2);
         s(2) s(3)];

    p = p_(3*(k-1)+1:3*k);

    P = [p(1) p(2);
         p(2) p(3)];

    V = e'*P*e;
    
    if k ~= N
        p = p_(3*k+1:3*(k+1));
                
        Pnext = [p(1) p(2);
                 p(2) p(3)];

        Vdot = e'*(P*A+A'*P)*e ...
            + 2*unc*(Ppv*epbar + Pv*edbar) ...
            + e'*(Pnext - P)*e / dt;
        rhodot = (rho(k+1) - rho(k))/dt;
    else
        Vdot = e'*(P*A+A'*P)*e ...
            + 2*unc*(Ppv*epbar + Pv*edbar);
        rhodot = 0;        
    end
    
    [Lrho1,Crho1] =     polynomial([e;epbar;edbar],monomialOrder);
    [Lrho3,Crho3] =     polynomial( e,monomialOrder);       
    [Lep,Cep] =         polynomial([e;epbar;edbar],monomialOrder);
    [Led,Ced] =         polynomial([e;epbar;edbar],monomialOrder);
    [Lepsign,Cepsign] = polynomial([e;epbar;edbar],monomialOrder);
    [Ledsign,Cedsign] = polynomial([e;epbar;edbar],monomialOrder);

    Lrho1 = replace(Lrho1,Crho1,value(Crho1_(:,k)));
    Lrho3 = replace(Lrho3,Crho3,value(Crho3_(:,k)));

    c1 = sos(rhodot - Vdot ...
             - Lrho1*(rho(k) - V) ...
             - Lep*(epbar^2 - e(1)*e(1)) ...
             - Led*(edbar^2 - e(2)*e(2)) ...
             - Lepsign*(epbar) ...
             - Ledsign*(edbar) ...
             );
    c2 = sos(Lepsign);
    c3 = sos(Ledsign);

    c6 = sos(1 - e'*S*e ...
             -Lrho3*(rho(k) - V));

    c7 = S >= 0;

    c8 = sos(Ppv - P(1,2));
    c9 = sos(Ppv + P(1,2));

    c10 = sos(Pv - P(2,2));
    c11 = sos(Pv + P(2,2));
   
    constraints = [constraints c1 c2 c3 c6 c7 c8 c9 c10 c11];

    vars  = [vars;Ppv;Pv;Cep;Ced;Cepsign;Cedsign]; 
    sVars = [sVars;s];
    sOut  = [sOut s];

    cost = cost + geomean(S);
end

sol = solvesos(constraints,-cost,[],[vars;sVars;rho;p_]);
solProblem = sol.problem;

rho = value(rho);
sOut = value(sOut);
p_ = value(p_);

end
