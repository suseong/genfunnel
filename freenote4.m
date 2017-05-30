clear
close all
clc

p = quadRPG();
numTraj = 11;

figure(1)
hold on

for idx = 1:numTraj
    xs = [0 0 0    2 0 0    0 0 0;
          0 0 0    2 0 0    0 0 0;
          0 0 0    2 0 0    0 0 0;
          0 0 0    2 0 0    0 0 0;
          0 0 0    2 0 0    0 0 0;
          0 0 0    2 0 0    0 0 0;
          0 0 0    2 0 0    0 0 0;
          0 0 0    2 0 0    0 0 0;
          0 0 0    2 0 0    0 0 0;
          0 0 0    2 0 0    0 0 0;          
          0 0 0    2 0 0    0 0 0];
      
    xf = [2 0.0 0    2 0 0    0 0 0;
          2 0.3 0    2 0 0    0 0 0;
          2 0.6 0    2 0 0    0 0 0;
          2 0.9 0    2 0 0    0 0 0;          
          2 1.2 0    2 0 0    0 0 0;
          2 1.5 0    2 0 0    0 0 0;          
          2 -0.3 0    2 0 0    0 0 0;
          2 -0.6 0    2 0 0    0 0 0;
          2 -0.9 0    2 0 0    0 0 0;          
          2 -1.2 0    2 0 0    0 0 0;          
          2 -1.5 0    2 0 0    0 0 0];
    duration = [1 1 1 1 1 1 1 1 1 1 1];

    [xtj,utj] = trajGen(p,xs(idx,:),xf(idx,:),duration(idx));
    xtraj{idx} = xtj;
    utraj{idx} = utj;
    
    tt_ = xtj.getBreaks;
    tt = 0:0.02:tt_(end);
    
    states = xtj.eval(tt);
    states = [states xtj.eval(tt_(end))];
    plot(states(2,:),states(1,:))
end

axis([-2 2 0 2])
xlabel('y')
ylabel('x')
axis equal

%%
idx = 2;

ts = xtraj{idx}.getBreaks();
tt = ts;

states = [];
inputs = [];
for i=1:size(tt,2)
    states = [states xtraj{idx}.eval(tt(i))];
    inputs = [inputs utraj{idx}.eval(tt(i))];
end

figure(1)
for i=1:9
   subplot(9,1,i)
   plot(tt,states(i,:))
   axis([ts(1) ts(end) -1 1])
   axis 'auto y'
end

figure(2)
for i=1:4
   subplot(4,1,i)
   plot(tt,inputs(i,:))
   axis([ts(1) ts(end) -1 1])
   axis 'auto y'
end

figure
plot(states(1,:),states(2,:))
axis equal

%%
Q = diag([10*ones(6,1); 50*ones(3,1)]);
R = 1*eye(4);

for idx = 1:numTraj
    xtj = xtraj{idx}; utj = utraj{idx};
    [c_,V__]=tvlqr1(p,xtj,utj,Q,R,Q);  % c, V_
    poly_ = taylorApprox(feedback(p,c_),xtj,utj,2); % poly
    c{idx} = c_; V_{idx} = V__; poly{idx} = poly_;
end

%%
t = msspoly('t',1);
ts = [];
for idx = 1:numTraj
    V0 = V_{idx};
    G = QuadraticLyapunovFunction(V0.getFrame,Q);

    poly_idx = poly{idx}.inStateFrame(V0.getFrame);
    x = V0.getFrame.getPoly;

    ts{idx} = V0.S.getBreaks;
    N = length(ts{idx});

    for i=1:floor(N)
       V{idx,i} = V0.getPoly(ts{idx}(i));
       f{idx,i} = poly_idx.getPolyDynamics(ts{idx}(i));
       f{idx,i} = subs(f{idx,i},poly_idx.getInputFrame.getPoly,[0;0;0;0]);
       dVdt{idx,i} = V0.getPolyTimeDeriv(ts{idx}(i));
       Vdot{idx,i} = diff(V{idx,i},x)*f{idx,i} + dVdt{idx,i};
    end
end

%%
for idx = 1:numTraj
    dts{idx} = diff(ts{idx});
    rho{idx} = 1*exp(-0.5*(ts{idx}-ts{idx}(1))/(ts{idx}(end)-ts{idx}(1)));
    rhodot{idx} = diff(rho{idx})./dts{idx};

    for i=1:length(ts{idx})-1
      m{idx}(i)=sampleCheck_(x,V{idx,i},Vdot{idx,i},rho{idx}(i),rhodot{idx}(i));
    end
%     if max(m{idx}) > 0
%         disp('failed to find feasible rho')
%     end
end

%%
for idx = 5:numTraj
    V_idx = []; Vdot_idx = [];
    for zz = 1:length(ts{idx}) 
        V_idx{zz} = V{idx,zz};
        Vdot_idx{zz} = Vdot{idx,zz};
    end
    N = length(ts{idx});
    rho_prev = zeros(N,1);

    for i=1:1
        L1 = findL_(x,V_idx,Vdot_idx,double(rho{idx}),double(rhodot{idx}));
        [rho_,rhodot_] = findRho_(x,V_idx,Vdot_idx,L1,dts{idx},rho{idx}(1));
        rho{idx} = double(rho_);rhodot{idx} = double(rhodot_);
        eps = norm(double(rho{idx}) - rho_prev);
        rho_prev = double(rho{idx});
        eps
    end
end

%%
idx = 1;
odefun{idx} = @(t,x)p.dynamics(t,x,utraj{idx}.eval(t)+c{idx}.D.eval(t)*(x-xtraj{idx}.eval(t)));
% xtraj_pos = xtraj{idx}.eval(ts{idx}(end));
% xtraj_pos(4:9) = zeros(6,1);
% idx = 2;
% odefun{idx} = @(t,x)p.dynamics(t,x,utraj{idx}.eval(t)+c{idx}.D.eval(t)*(x-xtraj{idx}.eval(t)-xtraj_pos));
% xtraj_pos = xtraj{idx}.eval(ts{idx}(end));
% xtraj_pos(4:9) = zeros(6,1);
% idx = 3;
% odefun{idx} = @(t,x)p.dynamics(t,x,utraj{idx}.eval(t)+c{idx}.D.eval(t)*(x-xtraj{idx}.eval(t)-xtraj_pos));

for idx = 2:numTraj
    xtraj_pos = zeros(9,1);
    for j=1:idx-1
        xtraj_pos = xtraj_pos + xtraj{j}.eval(ts{j}(end));
    end
    xtraj_pos(4:9) = zeros(6,1);
%     idx = 2;
    odefun{idx} = @(t,x)p.dynamics(t,x,utraj{idx}.eval(t)+c{idx}.D.eval(t)*(x-xtraj{idx}.eval(t)-xtraj_pos));
end
opts = odeset('RelTol',1e-6,'MaxStep',0.01); 

% refTraj = [];
% ttt = 0:0.03:ts(end);
% for i=1:length(ttt)
%     refTraj(:,i) = xtraj.eval(ttt(i));
% end
% sol = ode45(odefun,[ts(1) ts(end)],[0;0;0;2;0;0;0;0;0],opts);

%%
% figure(1)
% for i=1:9
%     subplot(9,1,i)
%     hold on
%     plot(ttt,refTraj(i,:),'g')
%     plot(sol.x,sol.y(i,:),':')
% end

%%
ang = -pi:0.3:pi+0.1;
xx = [];xx_ = [];yy = [];yy_ = [];zz = [];zz_ = [];

for idx = 1:numTraj
    N = length(ts{idx});
    for i=1:N-1
        x = V0.getFrame.getPoly;
        V_quad=.5*double(subs(diff(diff(V{idx,i},x)',x),x,0*x))/double(rho{idx}(i));
        invV = inv(V_quad);
        for j=1:length(ang)
            for k=1:length(ang)
                a = [sin(ang(j))*sin(ang(k));sin(ang(j))*cos(ang(k));cos(ang(j));zeros(6,1)];
                lambda = sqrt(a'*invV*a);
                temp = 1/lambda*invV*a;
                xx{idx,i}(j,k) = temp(1); yy{idx,i}(j,k) = temp(2); zz{idx,i}(j,k) = temp(3);
                xx_{idx,i}(j,k) = -temp(1); yy_{idx,i}(j,k) = -temp(2); zz_{idx,i}(j,k) = -temp(3);
            end
        end
    end
end

%%
figure(100);clf;
hold on
idx = 1;
N = length(ts{idx});
for i=2:N-1
    pos = xtraj{idx}.eval(ts{idx}(i));
    surf([xx{idx,i} xx_{idx,i}]+pos(1),[yy{idx,i} yy_{idx,i}]+pos(2),[zz{idx,i} zz_{idx,i}]+pos(3),'FaceColor',[0.8 0.8 0.8],'FaceAlpha',0.05,'LineStyle','none');
end
i=1;
pos = xtraj{idx}.eval(ts{idx}(i));
surf([xx{idx,i} xx_{idx,i}]+pos(1),[yy{idx,i} yy_{idx,i}]+pos(2),[zz{idx,i} zz_{idx,i}]+pos(3),'FaceColor',[0.8 0.8 0.1],'FaceAlpha',0.1,'LineStyle','none');

for idx = 2:numTraj

    N = length(ts{idx});
    for i=2:N-1        
        pos = xtraj{idx}.eval(ts{idx}(i));
        for j = 1:idx-1
            pos = pos + xtraj{j}.eval(ts{j}(end));
        end
        surf([xx{idx,i} xx_{idx,i}]+pos(1),[yy{idx,i} yy_{idx,i}]+pos(2),[zz{idx,i} zz_{idx,i}]+pos(3),'FaceColor',[0.8 0.8 0.8],'FaceAlpha',0.05,'LineStyle','none');
    end
    i=1;
    pos = xtraj{idx}.eval(ts{idx}(i));
    for j = 1:idx-1
        pos = pos + xtraj{j}.eval(ts{j}(end));
    end
    surf([xx{idx,i} xx_{idx,i}]+pos(1),[yy{idx,i} yy_{idx,i}]+pos(2),[zz{idx,i} zz_{idx,i}]+pos(3),'FaceColor',[0.8 0.8 0.1],'FaceAlpha',0.1,'LineStyle','none');
end
xlabel('x')
ylabel('y')
zlabel('z')
axis equal

%%
for k=1:100
    i = floor(rand(1)*size(xx{1,1},1));
    j = floor(rand(1)*size(xx{1,1},1));
    if i == 0
        i = 1;
    elseif j == 0
        j = 1;
    end
    x = V0.getFrame.getPoly;
    V_quad=.5*double(subs(diff(diff(V{1,1},x)',x),x,0*x))/double(rho{1}(1));    
    a = 10*randn(9,1);
    r = a'*V_quad*a;
    a = a / sqrt(r) + xtraj{1}.eval(0);
    
    idx = 1;
    sol{idx} = ode45(odefun{idx},[ts{idx}(1) ts{idx}(end)],a,opts);
    plot3(sol{idx}.y(1,:),sol{idx}.y(2,:),sol{idx}.y(3,:));
    for idx = 2:numTraj
        sol{idx} = ode45(odefun{idx},[ts{idx}(1) ts{idx}(end)],sol{idx-1}.y(:,end),opts);
        plot3(sol{idx}.y(1,:),sol{idx}.y(2,:),sol{idx}.y(3,:));
    end
    pause(0.1)
end
axis equal

%%
for idx = 1:numTraj-1
    N = length(ts{idx});
%     V_in = 0.5*double(diff(diff(V{idx,N}.getPoly(ts{idx}(1)),x)',x))/double(rho{idx}(1));
    V_1 = 0.5*double(diff(diff(V{idx,N},x)',x))/double(rho{idx}(end));
    V_2 = 0.5*double(diff(diff(V{idx+1,1},x)',x))/double(rho{idx}(1));

    x_1 = xtraj{idx}.eval(ts{idx}(end));
    x_2 = xtraj{idx+1}.eval(ts{idx}(1));

    B_ = inv(chol(V_2));
    b_ = -V_1*x_1;
    c_ = x_1'*V_1*x_1-1;

    l = 1e10;
    chad = [l-c_-b_'*inv(V_1)*b_ zeros(1,9) (x_2+inv(V_1)*b_)';
            zeros(9,1) l*eye(9) B_;
            x_2+inv(V_1)*b_ B_ inv(V_1)];

    composability = (eig(chad));
    minval = min(composability);
    if minval < 0
        disp('not composable');
    end
end

%%
% save test1.mat xtraj utraj Q R c V Vdot poly rho rhodot ts numTraj x;















