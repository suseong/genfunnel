clear
close all
clc

% 0  0  0  0  0  0
% 4  0  0  4 -2  0
% 12 0  0  4  2  0
% 20 0  0  4 -2  0

%  0  0  0  /  0  0  0
%  2 -2  0  /  0 -2  0
%  0 -4  0  / -2  0  0
% -2 -2  0  /  0  2  0
%  0  0  0  /  2  0  0

p = quadRPG();
x0 = Point(getStateFrame(p));

wypt = [ 0  0  0  0  0  0;
         2 -2  0  0 -2  0;
         0 -4  0 -2  0  0;
        -2 -2  0  0  2  0;
         0  0  0  2  0  0;
         2 -2  0  0 -2  0;
         0 -4  0 -2  0  0;
        -2 -2  0  0  2  0;
         0  0  0  0  0  0];

desStateTimeBuff = [];
desStateBuff = [];

%%
for k=1:8
x0.x = wypt(k,1);
x0.y = wypt(k,2);
x0.dx =  wypt(k,4);
x0.dy =  wypt(k,5);

xf = x0;
xf.x = wypt(k+1,1);
xf.y = wypt(k+1,2);
xf.dx = wypt(k+1,4);
xf.dy = wypt(k+1,5);

% x0.x = 0;
% x0.y = 0;
% x0.dx = 0;
% x0.dy = -2;
% xf = x0;
% xf.x = 3;
% xf.y = 0;
% xf.dx = 0;
% xf.dy = 2;

u0 = double([-1.0*p.m*p.g;0;0;0]);

N = 50;
minimum_duration = .1;
maximum_duration = 2;

prog = DircolTrajectoryOptimization(p,N,[minimum_duration maximum_duration]);  
prog = prog.addStateConstraint(ConstantConstraint(double(x0)),1);
prog = prog.addInputConstraint(ConstantConstraint(u0),1);
prog = prog.addStateConstraint(ConstantConstraint(double(xf)),N);
prog = prog.addInputConstraint(ConstantConstraint(u0),N);
prog = prog.addRunningCost(@cost);
prog = prog.addFinalCost(@finalCost);

tf0 = 1;                    % initial guess at duration 
traj_init.x = PPTrajectory(foh([0,tf0],[double(x0),double(xf)]));
traj_init.u = ConstantTrajectory(u0);

info=0;
while (info~=1)
  tic
  [xtraj,utraj,z,F,info] = prog.solveTraj(tf0,traj_init);
  toc
end

if k > 1
   desStateTimeBuff = [desStateTimeBuff xtraj.getBreaks+desStateTimeBuff(end)];
else
   desStateTimeBuff = xtraj.getBreaks;
end

desStateBuff = [desStateBuff xtraj.eval(xtraj.getBreaks)];

disp('-----------------------')

Q = diag([20*ones(6,1); 30*ones(3,1)]);
R = 1*eye(4);

[c,V_]=tvlqr1(p,xtraj,utraj,Q,R,Q);

i = k;
generateCSV(xtraj.getBreaks,xtraj.pp.coefs,['refTraj',num2str(i),'.csv']);
generateCSV(utraj.getBreaks,utraj.pp.coefs,['inputTraj',num2str(i),'.csv']);
generateCSV(c.D.getBreaks,c.D.pp.coefs,['controlTraj',num2str(i),'.csv']);
end


figure(1000)
for k=1:9
   subplot(9,1,k)
   plot(desStateTimeBuff,desStateBuff(k,:))
end

%%
t = 0:0.01:16;
tj = [];dtj = [];
for i=1:length(t)
   x = 2*sin(t(i)*2*pi/8); dx = 2*2*pi/8*cos(t(i)*2*pi/8);
   y = -2*(cos(t(i)*2*pi/8)-1); dy = 2*2*pi/8*(sin(t(i)*2*pi/8));
   
   tj = [tj [x;y]];
   dtj = [dtj [dx;dy]];
end

%%
figure(98)
subplot(4,1,1)
plot(t,tj(1,:))
subplot(4,1,2)
plot(t,tj(2,:))
subplot(4,1,3)
plot(t,dtj(1,:))
subplot(4,1,4)
plot(t,dtj(2,:))

%%
clear all
close all

p = quadRPG();

wypt = [   0    0    0   0   0   0   0   0   0;
        0.25    0    0 0.5   0   0   0   0   0;
        0.75  0.3    0 0.5   0   0   0   0   0;
        1.25  0.9    0 0.5   0   0   0   0   0;
        1.75  0.9    0 0.5   0   0   0   0   0;
        2.25  0.3    0 0.5   0   0   0   0   0;
        2.75    0    0 0.5   0   0   0   0   0;
        3.00    0    0   0   0   0   0   0   0];

desStateTimeBuff = [];
desStateBuff = [];
inputBuff = [];
trajNum = 7;

for k=1:trajNum
tf = 1;
delT = 0.05;

xs = wypt(k,:);xf = wypt(k+1,:);
[xtraj,utraj] = refTraj2(xs,xf,tf,delT,p);

if k > 1
   desStateTimeBuff = [desStateTimeBuff xtraj.getBreaks+desStateTimeBuff(end)];
else
   desStateTimeBuff = xtraj.getBreaks;
end

desStateBuff = [desStateBuff xtraj.eval(xtraj.getBreaks)];
inputBuff = [inputBuff utraj.eval(utraj.getBreaks)];

disp('-----------------------')

Q = diag([20*ones(6,1); 40*ones(3,1)]);
R = 1*eye(4);

[c,V_]=tvlqr1(p,xtraj,utraj,Q,R,Q);

i = k;
generateCSV(xtraj.getBreaks,xtraj.pp.coefs,['refTraj',num2str(i),'.csv']);
generateCSV(utraj.getBreaks,utraj.pp.coefs,['inputTraj',num2str(i),'.csv']);
generateCSV(c.D.getBreaks,c.D.pp.coefs,['controlTraj',num2str(i),'.csv']);

end

figure(1001)
for k=1:9
   subplot(9,1,k)
   plot(desStateTimeBuff,desStateBuff(k,:))
   grid on
end

figure(1002)
for k=1:4
   subplot(4,1,k)
   plot(desStateTimeBuff,inputBuff(k,:))
   grid on
end

%%
clear
close
clc

for j = 1:1
for k=1:5000
    y_list(:,k) = diag([2 2 5 5 10 10])*(rand(6,1)-0.5);   
end

stPts{j} = y_list([1 3 5],:);
edPts{j} = y_list([2 4 6],:);

stPts{j}(:,1) = zeros(3,1);

% figure(94);clf
% hold on
% plot3(stPts(1,:),stPts(2,:),stPts(3,:),'sq')
% plot3(edPts(1,:),edPts(2,:),edPts(3,:),'*')

% for k=1:length(stPts)
%    plot([stPts(1,k),edPts(1,k)],[stPts(2,k),edPts(2,k)],':k','linewidth',.1)
% end

id = 1;
idPack{j} = 1;
for i=1:10
    prevEdPts = edPts{j}(:,id);
    stFromPrevEdPts = stPts{j} - prevEdPts;
    for kk = 1:length(stPts{j})
        normPack(kk) = norm(stFromPrevEdPts(2:3,kk));
    end
    [~,id] = min(normPack);
    idPack{j}(i+1) = id;
end

% for kk=1:length(idPack)
%    plot3([stPts(1,idPack(kk)),edPts(1,idPack(kk))],[stPts(2,idPack(kk)),edPts(2,idPack(kk))],[stPts(3,idPack(kk)),edPts(3,idPack(kk))],'-k','linewidth',.1)
% end
% 
% for kk=1:length(idPack)-1
%    plot3([stPts(1,idPack(kk+1)),edPts(1,idPack(kk))],[stPts(2,idPack(kk+1)),edPts(2,idPack(kk))],[stPts(3,idPack(kk+1)),edPts(3,idPack(kk))],':k','linewidth',.1)
% end

end

desStateTimeBuff = [];
desStateBuff = [];
inputBuff = [];

p = quadRPG();

for k=1:length(idPack{1})
    xs = [stPts{1}(1,idPack{1}(k)) stPts{2}(1,idPack{2}(k)) 0 stPts{1}(2,idPack{1}(k)) stPts{2}(2,idPack{2}(k)) 0 stPts{1}(3,idPack{1}(k)) stPts{2}(3,idPack{2}(k)) 0];
    if(k == length(idPack{1}))
        xf = [edPts{1}(1,idPack{1}(k)) edPts{2}(1,idPack{2}(k)) 0 0 0 0 0 0 0];
    else
        xf = [edPts{1}(1,idPack{1}(k)) edPts{2}(1,idPack{2}(k)) 0 edPts{1}(2,idPack{1}(k)) edPts{2}(2,idPack{2}(k)) 0 edPts{1}(3,idPack{1}(k)) edPts{2}(3,idPack{2}(k)) 0];
    end
    
    tf = 2;
    delT = 0.05;

%     xs = wypt(k,:);
%     xf = wypt(k+1,:);
    [xtraj,utraj] = refTraj1(xs,xf,tf,delT,p);
    size(xtraj.pp.coefs)

    if k > 1
       desStateTimeBuff = [desStateTimeBuff xtraj.getBreaks+desStateTimeBuff(end)];
    else
       desStateTimeBuff = xtraj.getBreaks;
    end

    desStateBuff = [desStateBuff xtraj.eval(xtraj.getBreaks)];
    inputBuff = [inputBuff utraj.eval(utraj.getBreaks)];

    disp('-----------------------')

    Q = diag([20*ones(6,1); 40*ones(3,1)]);
    R = 1*eye(4);

    [c,V_]=tvlqr1(p,xtraj,utraj,Q,R,Q);

%     i = k;
%     generateCSV(xtraj.getBreaks,xtraj.pp.coefs,['refTraj',num2str(i),'.csv']);
%     generateCSV(utraj.getBreaks,utraj.pp.coefs,['inputTraj',num2str(i),'.csv']);
%     generateCSV(c.D.getBreaks,c.D.pp.coefs,['controlTraj',num2str(i),'.csv']);
end

figure(1001)
for k=1:9
   subplot(9,1,k)
   plot(desStateTimeBuff,desStateBuff(k,:))
end

figure(1002)
for k=1:4
   subplot(4,1,k)
   plot(desStateTimeBuff,inputBuff(k,:))
end

%%
% ts = xtraj.getBreaks();
% % tt = ts;
% tt = 0:delT:ts(end);
% ts = xtraj.getBreaks;

% tf = xtraj.getBreaks;
% tt = 0:delT:tf;
% ts = tt;
% tt = xtraj.getBreaks;
% % 
% states = [];
% inputs = [];
% for i=1:size(tt,2)
%     states = [states xtraj.eval(tt(i))];
%     inputs = [inputs utraj.eval(tt(i))];
% end
% 
% figure(1)
% for i=1:9
%    subplot(9,1,i)
%    plot(tt,states(i,:))
%    axis([tt(1) tt(end) -1 1])
%    axis 'auto y'
% end
% 
% figure(2)
% for i=1:4
%    subplot(4,1,i)
%    plot(tt,inputs(i,:))
%    axis([tt(1) tt(end) -1 1])
%    axis 'auto y'
% end
% 
% figure(3)
% plot3(states(2,:),states(1,:),states(3,:))
% axis equal

% %%
% stateName = p.getStateFrame.getCoordinateNames;
% ts = xtraj.getBreaks();
% 
% docNode = com.mathworks.xml.XMLUtils.createDocument('drake');
% 
% xtraj_node = docNode.createElement('xtraj');
% docNode.getDocumentElement.appendChild(xtraj_node);
% 
% entry_node = docNode.createElement(['t','0']);
% docNode.getDocumentElement.appendChild(entry_node);
% 
% num_node = docNode.createElement('num');
% num_text = docNode.createTextNode(num2str(length(ts)));
% num_node.appendChild(num_text);
% entry_node.appendChild(num_node);
% xtraj_node.appendChild(entry_node);
% 
% for z = 1:length(ts)-1
%     entry_node = docNode.createElement(['t',num2str(z)]);
%     docNode.getDocumentElement.appendChild(entry_node);
%     time_node = docNode.createElement('t');
%     time_text = docNode.createTextNode(num2str(ts(z)));
%     time_node.appendChild(time_text);
%     entry_node.appendChild(time_node);
%     
%     for i=1:9
%         sn = stateName{i};
%         for j=1:6
%             trajName = [sn,num2str(j)];
%             thisElement = docNode.createElement(trajName);
%             thisElement.appendChild(docNode.createTextNode(sprintf('%d',xtraj.pp.coefs(9*(z-1)+i,j))));
%             entry_node.appendChild(thisElement);            
%         end
%     end
%     xtraj_node.appendChild(entry_node);
% end
% xmlFileName = ['aaa','.xml'];
% xmlwrite(xmlFileName,docNode);
% type(xmlFileName);
% 
%%
tt = xtraj.getBreaks;
for i=1:length(tt)
    a = xtraj.dderiv(tt(i));
    a = a(1:3);
    u = utraj.eval(tt(i));
    ang = xtraj.eval(tt(i));
    ang = ang(7:9);
    R = rot(ang(1),ang(2),ang(3));
    e3 = [0;0;1];
    m = p.m;
    (a - p.g*e3 - u(1)/p.m*R*e3)'
    
    phi = ang(1);theta = ang(2);psi = ang(3);
    invQ = [ 1   (sin(theta)*sin(phi))/cos(theta)  (cos(phi)*sin(theta))/cos(theta);
             0                           cos(phi)                         -sin(phi);
             0                sin(phi)/cos(theta)              cos(phi)/cos(theta)];
   
    xdot = xtraj.deriv(tt(i));
    (xdot(7:9) - invQ*u(2:4))'
end

%%
odefun = @(t,x)p.dynamics(t,x,utraj.eval(t));
opts = odeset('RelTol',1e-6,'MaxStep',0.005);
sol = ode45(odefun,[ts(1) ts(end)],[0;0;0;2;0;0;0;0;0],opts);

figure(15);clf;
for i=1:9
   subplot(9,1,i)
   hold on
   plot(tt,states(i,:),'g','linewidth',3)
   plot(sol.x,sol.y(i,:),'k')
end

%%
% Q = diag([20*ones(6,1); 40*ones(3,1)]);
% R = 1*eye(4);
% 
% [c,V_]=tvlqr1(p,xtraj,utraj,Q,R,Q);

%%
poly = taylorApprox(feedback(p,c),xtraj,utraj,2);

t = msspoly('t',1);

num_x = p.getNumStates();
num_xc = p.getNumContStates();
num_u = p.getNumInputs();

V0 = V_;
G = QuadraticLyapunovFunction(V0.getFrame,Q);

poly = poly.inStateFrame(V0.getFrame);
x = V0.getFrame.getPoly;

ts = V0.S.getBreaks;
N = length(ts);

for i=1:floor(N)
   V{i} = V0.getPoly(ts(i));
   f{i} = poly.getPolyDynamics(ts(i));
   f{i} = subs(f{i},poly.getInputFrame.getPoly,[0;0;0;0]);
   dVdt{i} = V0.getPolyTimeDeriv(ts(i));
   Vdot{i} = diff(V{i},x)*f{i} + dVdt{i};
end

%%
i = k;
generateCSV(xtraj.getBreaks,xtraj.pp.coefs,['refTraj',num2str(i),'.csv']);
generateCSV(utraj.getBreaks,utraj.pp.coefs,['inputTraj',num2str(i),'.csv']);
generateCSV(c.D.getBreaks,c.D.pp.coefs,['controlTraj',num2str(i),'.csv']);

%%
V_quad=.5*double(subs(diff(diff(V{1},x)',x),x,0*x));

pE = 0.2; vE = 0.2; aE = 0.1;
E = [pE pE pE vE vE vE aE aE aE]';
initRegion = inv(diag(E)*diag(E));

minEig = -1e5;
rho_first = 0.5;

while minEig < 0
   eigs = eig(initRegion - V_quad / rho_first);
   minEig = min(eigs);
   if minEig < 0
       rho_first = rho_first + 0.1;
   end
end

dts = diff(ts);
rho = rho_first*exp(-0.0*(ts-ts(1))/(ts(end)-ts(1))); 

rhodot = diff(rho)./dts;

for i=1:length(ts)-1
  m(i)=sampleCheck_(x,V{i},Vdot{i},rho(i),rhodot(i));
end

if max(m) > 0
    disp('failed to find feasible rho')
end

%%
V_0 = V;
Vdot_0 = Vdot;
rho_prev = zeros(N,1);
eps = 1;

solProb = 0;
while eps  > 0.05 && solProb == 0
    [L1,L2] = findL_(x,V,Vdot,double(rho),double(rhodot),initRegion);
    [rho_temp,rhodot_temp,solProb] = findRho_(x,V,Vdot,L1,L2,dts,initRegion);
    if solProb == 0
        rho = double(rho_temp); rhodot = double(rhodot_temp);
    end
    eps
    eps = norm(double(rho) - rho_prev);
    rho_prev = double(rho);
    eps
    i = i+1
    pause(1);
end

%%
for i=1:N-1
  m(i)=sampleCheck_(x,V{i},Vdot{i},double(rho(i)),double(rhodot(i)));
end

if max(m) > 0
    disp('failed to find feasible rho')
end


%%
% for i=1:length(D)
%     D_{i} = double(D{i});
% end

%%
odefun = @(t,x)p.dynamics(t,x,utraj.eval(t)+c.D.eval(t)*(x-xtraj.eval(t))+c.y0.eval(t));
opts = odeset('RelTol',1e-6,'MaxStep',0.01);

refTraj = [];
ttt = 0:0.03:ts(end);
for i=1:length(ttt)
    refTraj(:,i) = xtraj.eval(ttt(i));
end
sol = ode45(odefun,[ts(1) ts(end)],[0;0;0;3;0;0;0;0;0],opts);

%%
figure(1)
for i=1:9
    subplot(9,1,i)
    hold on
    plot(ttt,refTraj(i,:),'g')
    plot(sol.x,sol.y(i,:),':')
end

%%
ang = -pi:0.3:pi+0.2;
xx = [];xx_ = [];yy = [];yy_ = [];zz = [];zz_ = [];
for i=1:N-1
    x = V0.getFrame.getPoly;
    V_quad=.5*double(subs(diff(diff(V{i},x)',x),x,0*x))/double(rho(i));
    invV = inv(V_quad);
    for j=1:length(ang)
        for k=1:length(ang)
            a = [sin(ang(j))*sin(ang(k));sin(ang(j))*cos(ang(k));cos(ang(j));zeros(6,1)];
            lambda = sqrt(a'*invV*a);
            temp = 1/lambda*invV*a;
            xx{i}(j,k) = temp(1); yy{i}(j,k) = temp(2); zz{i}(j,k) = temp(3);
            xx_{i}(j,k) = -temp(1); yy_{i}(j,k) = -temp(2); zz_{i}(j,k) = -temp(3);
        end
    end
end

%%
figure(100);clf;
hold on
% for i=2:N-1
%     pos = xtraj.eval(ts(i));
%     surf([xx{i} xx_{i}]+pos(1),[yy{i} yy_{i}]+pos(2),[zz{i} zz_{i}]+pos(3),'FaceColor',[0.8 0.8 0.8],'FaceAlpha',0.05,'LineStyle','none');
% end
i=N-1;
pos = xtraj.eval(ts(i));
surf([xx{i} xx_{i}]+pos(1),[yy{i} yy_{i}]+pos(2),[zz{i} zz_{i}]+pos(3),'FaceColor',[0.8 0.8 0.1],'FaceAlpha',0.1,'LineStyle','none');

xlabel('x')
ylabel('y')
zlabel('z')
axis equal

%%
pts = [];
for i=1:size(xx{1},1)
    for j=1:size(xx{1},2)
        chad = [xx{N-1}(i,j);yy{N-1}(i,j);zz{N-1}(i,j)];
        plot3(chad(1),chad(2),chad(3),'*')
        pts = [pts chad]; 
    end
end

[AAA,ccc] = MinVolEllipse(pts,0.001);

for i=1:length(pts)
   aaa = pts(:,i);
   bbb = aaa'*AAA*aaa;
   if bbb > 1
       bbb
   end
end



%%
for k=1:100
    i = floor(rand(1)*size(xx{1},1));
    j = floor(rand(1)*size(xx{1},1));
    if i == 0
        i = 1;
    elseif j == 0
        j = 1;
    end
    x = V0.getFrame.getPoly;
    V_quad=.5*double(subs(diff(diff(V{1},x)',x),x,0*x))/double(rho(1));    
    %     a = [xx{1}(i,j)*0.9 yy{1}(i,j)*0.9 zz{1}(i,j)*0.9 2 zeros(1,5)];
    a = 10*randn(9,1);
    r = a'*V_quad*a;
    a = a / sqrt(r) + xtraj.eval(0);
    sol = ode45(odefun,[ts(1) ts(end)],a,opts);
    plot3(sol.y(1,:),sol.y(2,:),sol.y(3,:));
    pause(0.1)
end
axis equal

%%
V_in   = 0.5*double(diff(diff(V{1},x)',x))/double(rho(1));
V_exit = 0.5*double(diff(diff(V{end},x)',x))/double(rho(end));

%%
% for i=1:length(V)
%     V_quad=.5*double(subs(diff(diff(V{i},x)',x),x,0*x))/double(rho(i));
%     [u,v,d] = svd(V_quad); % V_quad = u*v*u'
%     V_3d = u(1:3,:)*v*u(1:3,:)';
%     ch = chol(V_3d);
%     pt = xtraj.eval(ts(i));
%     for j=1:size(xs,1)
%         for k=1:size(xs,1)
%             temp = inv(ch)*[xs(j,k);ys(j,k);zs(j,k)];
%             
%             xe(j,k) = temp(1)+pt(1); ye(j,k) = temp(2)+pt(2); ze(j,k) = temp(3)+pt(3);
%         end
%     end
%     surf(xe,ye,ze,'FaceColor',[0.8 0.8 0.8],'FaceAlpha',0.1,'LineStyle','none');
% end
% 
% axis equal

%%
% figure(3)
% temp_v = rT.deriv(tt);
% temp_a = rT.dderiv(tt);
% for i=1:3
%    subplot(3,2,2*i-1)
%    plot(tt,temp_v(i,:));    
%    subplot(3,2,2*i)
%    plot(tt,temp_a(i,:));
% end


