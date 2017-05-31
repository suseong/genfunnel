checkDependency('yalmip');

rho = 0.1;
rhodot = 10;
Er = 0.1;
B = 15;

ep = sdpvar(3,1);
ev = sdpvar(3,1);
E = sdpvar(2,1);

e = [ep;ev];

[L1,coeff1_] = polynomial(e,2);
[L2,coeff2_] = polynomial(e,2);

% p = sdpvar(5,1);
% P = [   1    0    0 p(2)    0    0;
%         0    1    0    0 p(2)    0;
%         0    0 p(1)   0     0 p(3);
%      p(2)    0    0 p(4)    0    0;
%         0 p(2)    0    0 p(4)    0;
%         0    0 p(3)    0    0 p(5)];
P = eye(6);
A = [zeros(3,3) eye(3); -diag([10 10 20]) -diag([4 4 6])];    
    
V = 0.5*e'*P*e;
Vdot = 0.5*e'*(P*A+A'*P)*e;

c1 = sos(rhodot - Vdot - L1*(rho - V));

sol = solvesos(c1,[],[],coeff1_);





