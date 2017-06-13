function my_phase(IC,k)
[~,X] = ode45(@EOM,[0 0.01*k],IC);
u = X(:,1);
w = X(:,2);
plot(u(end),w(end),'x')
% plot(u,w)
xlabel('u')
ylabel('w')
grid
end

function dX = EOM(t, y)
dX = zeros(2,1);
u  = y(1);
w  = y(2);
% kp = 10; kd = 4;
kp = 15; kd = 6;
dX = [w; -kp*u-kd*w];
end