function my_phase(IC)
[~,X] = ode45(@EOM,[0 10],IC);
u = X(:,1);
w = X(:,2);
plot(u,w,'r')
xlabel('u')
ylabel('w')
grid
end

function dX = EOM(t, y)
dX = zeros(2,1);
u  = y(1);
w  = y(2);
A  = 1;
B  = 1;
kp = 10; kd = 4;
dX = [w; -kp*u-kd*w];
end