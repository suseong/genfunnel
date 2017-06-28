function your_phase()

IC = [];
xx = -2:0.5:2;
yy = -2:0.5:2;

for k=1:length(xx)
    for j = 1:length(yy)
        IC = [IC [xx(k);yy(j)]];
    end
end

IC = IC';
hold on
for ii = 1:length(IC(:,1))
    [~,X] = ode45(@EOM,[0 50],IC(ii,:));
    u = X(:,1);
    w = X(:,2);
    plot(u,w,':k')
end
xlabel('u')
ylabel('w')
grid
end
function dX = EOM(t, y)
dX = zeros(2,1);
u  = y(1);
w  = y(2);
kp = 10;
kd = 4;
dX = [w+0.05;-kp*(u+0.3)-kd*(w+0.05)];
end