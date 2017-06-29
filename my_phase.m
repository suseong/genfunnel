function my_phase(IC,k)

t = 0.01*k;

[~,X] = ode45(@EOM1,[0 t],IC);
u = X(:,1);
w = X(:,2);
plot(u(end),w(end),'*')

% [~,X] = ode45(@EOM2,[0 t],IC);
% u = X(:,1);
% w = X(:,2);
% plot(u(end),w(end),'*')

xlabel('u')
ylabel('w')
grid

axis([-0.5 0.5 -2 2]);
axis equal
grid on
% pause(0.01);

end

function dX = EOM1(t, y)
dX = zeros(2,1);
u  = y(1);
w  = y(2);
kp = 10; kd = 4;
% kp = 15; kd = 6;
% if t > 2
    dX = [w; -kp*(u)-kd*(w)];
% else
%     dX = [w; -kp*(u)-kd*(w)-1];
% end
% sig = sign(dX(2));
% dX(2) = dX(2) + sig;
end

function dX = EOM2(t, y)
dX = zeros(2,1);
u  = y(1);
w  = y(2);
kp = 10; kd = 4;
% kp = 15; kd = 6;
% if t > 2
    dX = [w; -kp*(u)-kd*(w)-0.5];
% else
%     dX = [w; -kp*(u)-kd*(w)+1];
% end
% sig = sign(dX(2));
% dX(2) = dX(2) - sig;
end