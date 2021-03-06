clear all
% close all
clc

e = sdpvar(6,1);
dt = 0.05;

% initRegion = diag(1./([1 1 1 2 2 2]/2).^2);
initRegion = diag(1./([1 1 1 2 2 2]/5).^2);

%%
for kkk = 1:1

% Er = 0.052;
Er = 0.052;
ar = 15;

Kp = diag([10 10 15]);
Kd = diag([4 4 5]);
A = [zeros(3,3) eye(3); -Kp -Kd];    

% P = lyap(A',-eye(6));
p = [2.1141 3.5526 0.1433 0.1358 0.1632 0.1696];
P = [p(1)   0    0  p(3)   0    0;
    0  p(1)   0    0  p(3)   0;
    0    0  p(2)   0    0  p(4);
    p(3)   0    0  p(5)   0    0;
    0  p(3)   0    0  p(5)   0;
    0    0  p(4)   0    0  p(6)];

P = P/P(1,1)*initRegion(1,1);

N = 225;

unc = 2.5;

%%
rho = sdpvar(1,1);
rhodot = sdpvar(1,1);
epbar = sdpvar(1,1);
edbar = sdpvar(1,1);
Ppv = sdpvar(1,1);
Pv = sdpvar(1,1);

maxKp = max(max(Kp));                                                     
maxKd = max(max(Kd));

V = e'*P*e;
Vdot = e'*(P*A+A'*P)*e ...
       + 2*(unc + Er*(maxKp*epbar + maxKd*edbar + ar))*(Ppv*epbar + Pv*edbar);

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
    [Lrho,Crho] =       polynomial([e;epbar;edbar],monomialOrder);
    [Lep,Cep] =         polynomial([e;epbar;edbar],monomialOrder);
    [Led,Ced] =         polynomial([e;epbar;edbar],monomialOrder);
    [Lepsign,Cepsign] = polynomial([e;epbar;edbar],monomialOrder);
    [Ledsign,Cedsign] = polynomial([e;epbar;edbar],monomialOrder);
    
    vars = [Crho;Cep;Ced;Cepsign;Cedsign];
    vars = [vars;rhodot;Ppv;Pv];
    
    c1 = sos(rhodot - Vdot ...
             - Lrho*(rhoTemp - V) ...
             - Lep*(epbar^2 - e(1:3)'*e(1:3)) ...
             - Led*(edbar^2 - e(4:6)'*e(4:6)) ...
             - Lepsign*(epbar) ...
             - Ledsign*(edbar) ...
             );
    c2 = sos(Lepsign);
    c3 = sos(Ledsign);

    c4 = sos(Ppv - P(1,4));
    c5 = sos(Ppv + P(1,4));
    c6 = sos(Ppv - P(3,6));
    c7 = sos(Ppv + P(3,6));
    
    c8 = sos(Pv - P(4,4));
    c9 = sos(Pv + P(4,4));
    c10 = sos(Pv - P(6,6));
    c11 = sos(Pv + P(6,6));
        
    constraints = [c1 c2 c3 c4 c5 c6 c7 c8 c9 c10 c11];
    sol = solvesos(constraints,rhodot,[],vars);

    rhoTemp = rhoTemp + value(rhodot)*dt;
    rhoCont(i+1) = rhoTemp;
    i
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
[coeffL1,coeffL3,S_] = findL(dt,P_temp,Q_temp,rhoCont,rhodot,Kp,Kd,Er,ar,unc);

chk = size(find(sum(coeffL1) == 0),2);

[rho,sVars,p,solProblem] = findRho(dt,A,coeffL1,coeffL3,initRegion,Kp,Kd,Er,ar,unc);

%%
    SSS{kkk} = sVars;
    RHO{kkk} = rho;
    PPP{kkk} = p;
% end
end

%%
ang = linspace(-pi,pi,40);
for jj = 1:50
    figure(101);clf;
    axis equal
    hold on
    p = SSS{1}(:,jj);
% p = zeros(6,1);
    P = [p(1)   0    0  p(3)   0    0;
           0  p(1)   0    0  p(3)   0;
           0    0  p(2)   0    0  p(4);
         p(3)   0    0  p(5)   0    0;
           0  p(3)   0    0  p(5)   0;
           0    0  p(4)   0    0  p(6)];

    kk = 1;
    p1 = [P(kk,kk) P(kk,kk+3);P(kk+3,kk) P(kk+3,kk+3)];
    invp1 = inv(sqrtm(p1));
%     p2 = [initRegion(kk,kk) initRegion(kk,kk+3);initRegion(kk+3,kk) initRegion(kk+3,kk+3)];
%     invp2 = inv(sqrtm(p2));
    for k=1:length(ang)
        xx = invp1*[cos(ang(k));sin(ang(k))];
        plot(xx(1),xx(2),'.','markersize',15)
%         yy = invp1*[cos(ang(k));sin(ang(k))];
%         my_phase(xx,300);
    end
   axis([-0.5 0.5 -2 2]*4);
%    axis equal
   grid on
   jj
   pause(0.1);
end

%%
p = PPP{1};
p_ = reshape(p,6,length(p)/6);
rho_ = [rho';rho';rho';rho';rho';rho'];
p_ = p_./rho_;
pdiff = diff(p_')';
pdiffnorm = [];

for k=1:size(pdiff,2)
   pdiffnorm(k) = norm(pdiff(:,k));
end

%%
figure(11111)
subplot(2,1,1)
plot(rho)
subplot(2,1,2)
plot(pdiffnorm(1:end-1))


%%
ang = linspace(-pi,pi,40);
for jj = 1:225
    figure(101);clf;
    axis equal
    hold on
    p = SSS{1}(:,jj);

    P = [p(1)   0    0  p(3)   0    0;
           0  p(1)   0    0  p(3)   0;
           0    0  p(2)   0    0  p(4);
         p(3)   0    0  p(5)   0    0;
           0  p(3)   0    0  p(5)   0;
           0    0  p(4)   0    0  p(6)];

    kk = 1;
    p1 = [P(kk,kk) P(kk,kk+3);P(kk+3,kk) P(kk+3,kk+3)];
    invp1 = inv(sqrtm(p1));
%     p2 = [initRegion(kk,kk) initRegion(kk,kk+3);initRegion(kk+3,kk) initRegion(kk+3,kk+3)];
%     invp2 = inv(sqrtm(p2));
    for k=1:length(ang)
        xx = invp1*[cos(ang(k));sin(ang(k))];
        plot(xx(1),xx(2),'.','markersize',15)
%         yy = invp1*[cos(ang(k));sin(ang(k))];
%         my_phase(xx,300);
    end
   
    p = SSS_{1}(:,jj);

    P = [p(1)   0    0  p(3)   0    0;
           0  p(1)   0    0  p(3)   0;
           0    0  p(2)   0    0  p(4);
         p(3)   0    0  p(5)   0    0;
           0  p(3)   0    0  p(5)   0;
           0    0  p(4)   0    0  p(6)];

    kk = 1;
    p1 = [P(kk,kk) P(kk,kk+3);P(kk+3,kk) P(kk+3,kk+3)];
    invp1 = inv(sqrtm(p1));
%     p2 = [initRegion(kk,kk) initRegion(kk,kk+3);initRegion(kk+3,kk) initRegion(kk+3,kk+3)];
%     invp2 = inv(sqrtm(p2));
    for k=1:length(ang)
        xx = invp1*[cos(ang(k));sin(ang(k))];
        plot(xx(1),xx(2),'.','markersize',15)
%         yy = invp1*[cos(ang(k));sin(ang(k))];
%         my_phase(xx,300);
    end
   axis([-0.5 0.5 -2 2]*4);
%    axis equal
   grid on
   jj
   
   pause(0.01);
end





