function [rho_,s_,p_,sol] = findPrho(A,Cvdot_,Cs_,Kp,Kv,Er,ar,unc)

monomialOrder = 2;

maxKp = max(max(Kp));
maxKv = max(max(Kv));

e = sdpvar(6,1);
rho = sdpvar(1,1);
p = sdpvar(6,1);
s = sdpvar(6,1);
Ppv = sdpvar(1,1);
Pv = sdpvar(1,1);
epbar = sdpvar(1,1);
evbar = sdpvar(1,1);

P = [p(1)   0    0  p(3)   0    0;
       0  p(1)   0    0  p(3)   0;
       0    0  p(2)   0    0  p(4);
     p(3)   0    0  p(5)   0    0;
       0  p(3)   0    0  p(5)   0;
       0    0  p(4)   0    0  p(6)];

S = [s(1)   0    0  s(3)   0    0;
       0  s(1)   0    0  s(3)   0;
       0    0  s(2)   0    0  s(4);
     s(3)   0    0  s(5)   0    0;
       0  s(3)   0    0  s(5)   0;
       0    0  s(4)   0    0  s(6)];

V = e'*P*e;
Vdot = e'*(P*A+A'*P)*e ...
       + 2*(unc + Er*(unc + maxKp*epbar + maxKv*evbar + ar))*(Ppv*epbar + Pv*evbar);

[Lvdot,Cvdot] =     polynomial([e;epbar;evbar],monomialOrder);
[Lep,Cep] =         polynomial([e;epbar;evbar],monomialOrder);
[Lev,Cev] =         polynomial([e;epbar;evbar],monomialOrder);
[Lepsign,Cepsign] = polynomial([e;epbar;evbar],monomialOrder);
[Levsign,Cevsign] = polynomial([e;epbar;evbar],monomialOrder);
[Ls,Cs] =           polynomial([e;epbar;evbar],monomialOrder);

Lvdot = replace(Lvdot,Cvdot,value(Cvdot_));
Ls = replace(Ls,Cs,value(Cs_));

c1 = sos(rho - V ... 
         - Lvdot*(-Vdot) ...
         - Lep*(epbar^2 - e(1:3)'*e(1:3)) ...
         - Lev*(evbar^2 - e(4:6)'*e(4:6)) ...
         - Lepsign*(epbar) ...
         - Levsign*(evbar) ...
         );
     
c2 = sos(Lepsign);
c3 = sos(Levsign);

c4 = sos(1 - e'*S*e ...
         - Ls*(rho - V));
c6 = S >= 0;
c7 = P >= 0;

c8 = sos(Ppv - P(1,4));
c9 = sos(Ppv + P(1,4));
c10 = sos(Ppv - P(3,6));
c11 = sos(Ppv + P(3,6));

c12 = sos(Pv - P(4,4));
c13 = sos(Pv + P(4,4));
c14 = sos(Pv - P(6,6));
c15 = sos(Pv + P(6,6));

constraints = [c1,c2,c3,c4,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15];
vars = [p;s;Ppv;Pv;rho;Cep;Cev;Cepsign;Cevsign];

sol = solvesos(constraints,-geomean(S),[],vars);

rho_ = value(rho);
s_ = value(s);
p_ = value(p);

end