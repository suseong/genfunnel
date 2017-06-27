clear all
close all
clc

%% params
Er = 0.03;
ar = 9.8 + 5;
unc = 1.8;

Kp = diag([10 10 15]);
Kv = diag([4 4 6]);
A = [zeros(3,3) eye(3); -Kp -Kv];
P = lyap(A',-eye(6));
P = P/P(1,1);

maxKp = max(max(Kp));
maxKv = max(max(Kv));

Ppv_ = max(max(P(1:3,4:6)));
Pv_  = max(max(P(4:6,4:6)));

%% variables
rho = sdpvar(1,1);
epbar = sdpvar(1,1);
evbar = sdpvar(1,1);
e = sdpvar(6,1);

%%
V = e'*P*e;
Vdot = e'*(P*A+A'*P)*e ...
       + 2*(unc + Er*(unc + maxKp*epbar + maxKv*evbar + ar))*(Ppv_*epbar + Pv_*evbar);
   
monomialOrder = 2;

[Lvdot,Cvdot] =     polynomial([e;epbar;evbar],monomialOrder);
[Lep,Cep] =         polynomial([e;epbar;evbar],monomialOrder);
[Lev,Cev] =         polynomial([e;epbar;evbar],monomialOrder);
[Lepsign,Cepsign] = polynomial([e;epbar;evbar],monomialOrder);
[Levsign,Cevsign] = polynomial([e;epbar;evbar],monomialOrder);

c1 = sos(rho - V ...
         - Lvdot*(-Vdot) ...
         - Lep*(epbar^2 - e(1:3)'*e(1:3)) ...
         - Lev*(evbar^2 - e(4:6)'*e(4:6)) ...
         - Lepsign*(epbar) ...
         - Levsign*(evbar) ...
         );
     
c2 = sos(Lvdot);
c3 = sos(Lepsign);
c4 = sos(Levsign);

constraints = [c1 c3 c4];
vars = [rho;Cvdot;Cep;Cev;Cepsign;Cevsign];
sol = solvesos(constraints,rho,[],vars);

%%
s = sdpvar(6,1);
S = [s(1)   0    0  s(3)   0    0;
       0  s(1)   0    0  s(3)   0;
       0    0  s(2)   0    0  s(4);
     s(3)   0    0  s(5)   0    0;
       0  s(3)   0    0  s(5)   0;
       0    0  s(4)   0    0  s(6)];

[Ls,Cs] = polynomial([e;epbar;evbar],monomialOrder);
c1 = sos(1 - e'*S*e ...
         - Ls*(value(rho) - V) ...
        );
c2 = sos(Ls);
sol = solvesos([c1 c2],-geomean(S),[],[s;Cs]);

%%
[a,b,c,d] = findPrho(A,Cvdot,Cs,Kp,Kv,Er,ar,unc);













