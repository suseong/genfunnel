clear all
close all
clc

%% params
Kp = 10;
Kv = 4;
unc = 0.0;
dt = 0.01;
initRegion = diag(1./[0.2 0.3].^2);

A = [0 1; -Kp -Kv];
P = lyap(A',-eye(2));
P = P/P(1,1)*initRegion(1,1);

% P = initRegion;

Ppv_ = P(1,2);
Pv_  = P(2,2);

%% variables
rho = sdpvar(1,1);
rhodot = sdpvar(1,1);
epbar = sdpvar(1,1);
evbar = sdpvar(1,1);
Ppv = sdpvar(1,1);
Pv = sdpvar(1,1);
e = sdpvar(2,1);

%%
V = e'*P*e;
Vdot = e'*(P*A+A'*P)*e + 2*unc*(Ppv_*epbar + Pv_*evbar);
   
N = 200;
   
monomialOrder = 2;

[L_init,coeff_init] = polynomial(e,monomialOrder);

c1 = sos(L_init);
c2 = sos(rho - V ...
         - L_init*(1-e'*initRegion*e));
     
sol = solvesos([c1 c2],rho,[],coeff_init);
rhoInit = value(rho);

rhoTemp = rhoInit;
rhoCont = rhoInit;

%%
for i=1:N
    [Lrho,Crho] =       polynomial([e;epbar;evbar],monomialOrder);
    [Lep,Cep] =         polynomial([e;epbar;evbar],monomialOrder);
    [Led,Ced] =         polynomial([e;epbar;evbar],monomialOrder);
    [Lepsign,Cepsign] = polynomial([e;epbar;evbar],monomialOrder);
    [Ledsign,Cedsign] = polynomial([e;epbar;evbar],monomialOrder);
    
    vars = [Crho;Cep;Ced;Cepsign;Cedsign];
    vars = [vars;rhodot;Ppv;Pv];
    
    c1 = sos(rhodot - Vdot ...
             - Lrho*(rhoTemp - V) ...
             - Lep*(epbar^2 - e(1)*e(1)) ...
             - Led*(evbar^2 - e(2)'*e(2)) ...
             - Lepsign*(epbar) ...
             - Ledsign*(evbar) ...
             );
    c2 = sos(Lepsign);
    c3 = sos(Ledsign);

    c4 = sos(Ppv - P(1,2));
    c5 = sos(Ppv + P(1,2));
    
    c6 = sos(Pv - P(2,2));
    c7 = sos(Pv + P(2,2));
        
    constraints = [c1 c2 c3 c4 c5 c6 c7];
    sol = solvesos(constraints,rhodot,[],vars);

    rhoTemp = rhoTemp + value(rhodot)*dt;
    rhoCont(i+1) = rhoTemp;
    disp([num2str(i),' ',num2str(value(rhodot))]);
end

%%
ts = linspace(0,dt*(N-1),N);
rhodot = diff(rhoCont)/dt;

P_temp = [];
Q_temp = [];

for k = 1:N+1
    P_temp{k} = P;
    Q_temp{k} = P*A + A'*P;
end

%%
[coeffL1,coeffL3,S_] = findL2(P_temp,Q_temp,rhoCont,rhodot,unc);
chk = size(find(sum(coeffL1) == 0),2);

%%
[rho,sVars,p,solProblem] = findRho2(dt,A,coeffL1,coeffL3,initRegion,unc);

%%
ang = -pi:0.2:pi;
for jj = 100:N-2
    figure(101);clf;
    hold on
    p = sVars(:,jj);
    P = [p(1) p(2);
         p(2) p(3)];

    kk = 1;
    invp1 = inv(sqrtm(P));
    p2 = initRegion;
    invp2 = inv(sqrtm(p2));
    for k=1:length(ang)
        xx = invp1*[cos(ang(k));sin(ang(k))];
        plot(xx(1),xx(2),'.','markersize',15)
        yy = invp2*[cos(ang(k));sin(ang(k))];
%         plot(yy(1),yy(2),'*','markersize',15)        
        my_phase(yy,jj);
    end
   axis([-0.1 0.1 -0.2 0.2]);
   axis equal
   grid on
   jj
   pause(0.01);
end


%%
ang = -pi:0.2:pi;
for jj = 198:198
    figure(101);clf;
    hold on
    p = sVars(:,jj);
    P = [p(1) p(2);
         p(2) p(3)];

    kk = 1;
    invp1 = inv(sqrtm(P));
    p2 = initRegion;
    invp2 = inv(sqrtm(p2));
    for k=1:length(ang)
        xx = invp1*[cos(ang(k));sin(ang(k))];
        plot(xx(1),xx(2),'.','markersize',15)
        yy = invp2*[-cos(ang(k));-sin(ang(k))];
%         plot(yy(1),yy(2),'*','markersize',15)        
        my_phase(yy,jj);
    end
   axis([-0.5 0.5 -2 2]);
   axis equal
   grid on
   jj
   pause(0.01);
end

