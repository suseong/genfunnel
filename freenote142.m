clear all
close all

checkDependency('yalmip');

% rho = .1;
Er = 0.1;
B = 15;

ep = sdpvar(3,1);
ev = sdpvar(3,1);
% E = sdpvar(2,1);

e = [ep;ev];

[L1,coeff1] = polynomial(e,3);
[L_init,coeff_init] = polynomial(e,3);

rhodot = sdpvar(1,1);
rho = sdpvar(1,1);
p = sdpvar(5,1);

initRegion = diag(1./[1 1 1 2 2 2].^2);
A = [zeros(3,3) eye(3); -diag([10 10 15]) -diag([4 4 6])];    

% P1 = 0.9030 0.1794 0.0970;
% P1 = 0.9437 0.1272 0.0563;
% R2 = [cosd(ang(k)) -sind(ang(k));sind(ang(k)) cosd(ang(k))];
% mainAxis = [4 1];
% P1 = R2'*diag(mainAxis)^2/norm(mainAxis)^2*R2;

P = lyap(A',-0.1*eye(6));
P = P/P(1,1);

%%
% P(1:3,1:3) = eye(3); P(1:3,4:6) = eye(3); 
% P(4:6,1:3) = eye(3); P(4:6,4:6) = eye(3);

% P(1:3,1:3) = 0.903*eye(3); P(1:3,4:6) = 0.18*eye(3); 
% P(4:6,1:3) = 0.18*eye(3); P(4:6,4:6) = 0.097*eye(3);

V = 0.5*e'*P*e;
Vdot = 0.5*e'*(P*A+A'*P)*e;

c2 = sos(L_init);
c3 = sos(rho - V - L_init*(1-e'*initRegion*e));
sol = solvesos([c2 c3],rho,[],coeff_init);

rho_ = value(rho);
c1 = sos(rhodot - Vdot - L1*(rho_ - V));
sol = solvesos(c1,rhodot,[],coeff1);

% c2 = sos(L_init);
% c3 = sos(rho_ - V - L_init*(1-e'*initRegion*e));
% sol = solvesos([c1 c2 c3],rhodot,[],[coeff1;coeff_init]);

[L2,coeff2] = polynomial(e,3);
L2 = replace(L2,coeff2,value(coeff1));

[L3,coeff3] = polynomial(e,3);
L3 = replace(L3,coeff3,value(coeff_init));

P_ = [   1    0    0 p(2)    0    0;
         0    1    0    0 p(2)    0;
         0    0 p(1)   0     0 p(3);
      p(2)    0    0 p(4)    0    0;
         0 p(2)    0    0 p(4)    0;
         0    0 p(3)    0    0 p(5)];

V_ = 0.5*e'*P_*e;
Vdot_ = 0.5*e'*(P_*A+A'*P_)*e;

c4 = sos(rhodot - Vdot_ - L2*(rho - V_));
c5 = sos(rho - V - L3*(1-e'*initRegion*e));

sol = solvesos([c4 c5],rhodot,[],[p;rho]);

%%
P = [   1    0    0        value(p(2))    0    0;
        0    1    0           0 value(p(2))    0;
        0    0 value(p(1))   0     0 value(p(3));
     value(p(2))    0    0 value(p(4))    0    0;
        0 value(p(2))    0    0 value(p(4))    0;
        0    0 value(p(3))    0    0 value(p(5))];
    
V = 0.5*e'*P*e;
Vdot = 0.5*e'*(P*A+A'*P)*e;

c2 = sos(L_init);
c3 = sos(rho - V - L_init*(1-e'*initRegion*e));
sol = solvesos([c2 c3],rho,[],coeff_init);

rho_ = value(rho);
c1 = sos(rhodot - Vdot - L1*(rho_ - V));
c2 = sos(L_init);
c3 = sos(rho_ - V - L_init*(1-e'*initRegion*e));
sol = solvesos([c1 c2 c3],rhodot,[],[coeff1;coeff_init]);

[L2,coeff2] = polynomial(e,3);
L2 = replace(L2,coeff2,value(coeff1));

[L3,coeff3] = polynomial(e,3);
L3 = replace(L3,coeff3,value(coeff_init));

P_ = [   1    0    0 p(2)    0    0;
         0    1    0    0 p(2)    0;
         0    0 p(1)   0     0 p(3);
      p(2)    0    0 p(4)    0    0;
         0 p(2)    0    0 p(4)    0;
         0    0 p(3)    0    0 p(5)];

V_ = 0.5*e'*P_*e;
Vdot_ = 0.5*e'*(P_*A+A'*P_)*e;

c4 = sos(rhodot - Vdot_ - L2*(rho - V_));
c5 = sos(rho - V - L3*(1-e'*initRegion*e));

sol = solvesos([c4 c5],rhodot,[],[p;rho]);


% c1 = sos(rhodot - Vdot - L1*(rho - V));
% sol = solvesos(c1,rhodot,[],coeff1);
% 
% [L2,coeff2] = polynomial(e,3);
% L2 = replace(L2,coeff2,value(coeff1));
% 
% P_ = [   1    0    0 p(2)    0    0;
%         0    1    0    0 p(2)    0;
%         0    0 p(1)   0     0 p(3);
%      p(2)    0    0 p(4)    0    0;
%         0 p(2)    0    0 p(4)    0;
%         0    0 p(3)    0    0 p(5)];
% 
% V_ = 0.5*e'*P_*e;
% Vdot_ = 0.5*e'*(P_*A+A'*P_)*e;
% 
% c2 = sos(rhodot - Vdot_ - L2*(rho - V_));
% sol = solvesos(c2,rhodot,[],p);
value(rhodot)

%%
X = -2:0.5:2;
dX = -2:0.5:2;

figure(142);clf;
hold on
for k=1:length(X)
    for j=1:length(dX)
        IC = [X(k) dX(j)];
        my_phase(IC);
        axis equal;
        pause(0.01);
    end
end

ang = 30;
R2 = [cosd(ang) -sind(ang);sind(ang) cosd(ang)];
mainAxis = [1 100];
P1 = R2'*diag(mainAxis)/norm(mainAxis)*R2;

%%
figure(100);clf;
P = reshape(double(sVars),6,6);
for kk=1:3
subplot(1,3,kk)
p1 = [P(kk,kk) P(kk,kk+3);P(kk+3,kk) P(kk+3,kk+3)]/(value(rho));
p2 = [initRegion(kk,kk) initRegion(kk,kk+3);initRegion(kk+3,kk) initRegion(kk+3,kk+3)];
invp1 = inv(sqrtm(p1));
invp2 = inv(sqrtm(p2));
ang = -pi:0.02:pi;
hold on
for k=1:length(ang)
   xx = invp1*[cos(ang(k));sin(ang(k))];
   plot(xx(1),xx(2),'.') 
   yy = invp2*[cos(ang(k));sin(ang(k))];
   plot(yy(1),yy(2),'.')    
end
axis equal

end




