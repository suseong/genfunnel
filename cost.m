function [g,dg] = cost(dt,x,u)

R = diag([1 1 1 1]);
% Q = diag([0 0 0 0 0 0 5 5 5]);
g = u'*R*u;
dg = [zeros(1,1+size(x,1)),2*u'*R];
% dg = [2*x'*Q,2*u'*R];

end
