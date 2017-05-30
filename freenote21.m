clear
close all
clc

p = simplePlant();
% 
% wypt = [   0    0    0   0   0   0   0   0   0;
%         0.25    0    0 0.5   0   0   0   0   0;
%         0.75  0.3    0 0.5   0   0   0   0   0;
%         1.25  0.9    0 0.5   0   0   0   0   0;
%         1.75  0.9    0 0.5   0   0   0   0   0;
%         2.25  0.3    0 0.5   0   0   0   0   0;
%         2.75    0    0 0.5   0   0   0   0   0;
%         3.00    0    0   0   0   0   0   0   0];
% 
% desStateTimeBuff = [];
% desStateBuff = [];
% inputBuff = [];
% trajNum = 7;

%%
for k=6:trajNum
% for k=1:1

tf = 1;
delT = 0.05;

xs = wypt(k,:);xf = wypt(k+1,:);
[xtraj,utraj] = refTraj2(xs,xf,tf,delT,p);

% if k > 6
%    desStateTimeBuff = [desStateTimeBuff xtraj.getBreaks+desStateTimeBuff(end)];
% else
%    desStateTimeBuff = xtraj.getBreaks;
% end
% 
% desStateBuff = [desStateBuff xtraj.eval(xtraj.getBreaks)];
% inputBuff = [inputBuff utraj.eval(utraj.getBreaks)];

disp('-----------------------')

Q = diag([20*ones(6,1); 40*ones(3,1)]);
R = 1*eye(4);

[c,V_]=tvlqr1(p,xtraj,utraj,Q,R,Q);
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
V = []; Vdot = [];
for i=1:floor(N)
   V{i} = V0.getPoly(ts(i));
   f{i} = poly.getPolyDynamics(ts(i));
   f{i} = subs(f{i},poly.getInputFrame.getPoly,[0;0;0;0]);
   dVdt{i} = V0.getPolyTimeDeriv(ts(i));
   Vdot{i} = diff(V{i},x)*f{i} + dVdt{i};
end
V_quad=.5*double(subs(diff(diff(V{1},x)',x),x,0*x));

pE = 0.2; vE = 0.2; aE = 0.2;
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
cnt = 1;
solProb = 0;
while eps  > 0.05 && solProb == 0
    [L1,L2] = findL_(x,V,Vdot,double(rho),double(rhodot),initRegion);
    [rho_temp,rhodot_temp,solProb] = findRho_(x,V,Vdot,L1,L2,dts,initRegion);
    if solProb == 0
        rho = double(rho_temp); rhodot = double(rhodot_temp);
    end
    eps = norm(double(rho) - rho_prev);
    rho_prev = double(rho);
    eps
    cnt = cnt + 1
end
%%
ang = -pi:0.5:pi+0.1;
xx = [];xx_ = [];yy = [];yy_ = [];zz = [];zz_ = [];
ptsOnEllipsoid = [];
coverEllipsoid = [];
centerOfEllipsoid = [];
dataForRVIZ = [];
for i=1:N-1
    x = V0.getFrame.getPoly;
    V_quad=.5*double(subs(diff(diff(V{i},x)',x),x,0*x))/double(rho(i));
    invV = inv(V_quad);
    ptsOnEllipsoid_ = [];
    for j=1:length(ang)
        for kk=1:length(ang)
            a = [sin(ang(j))*sin(ang(kk));sin(ang(j))*cos(ang(kk));cos(ang(j));zeros(6,1)];
            lambda = sqrt(a'*invV*a);
            temp = 1/lambda*invV*a;
            xx{i}(j,kk) = temp(1); yy{i}(j,kk) = temp(2); zz{i}(j,kk) = temp(3);
            xx_{i}(j,kk) = -temp(1); yy_{i}(j,kk) = -temp(2); zz_{i}(j,kk) = -temp(3);
            ptsOnEllipsoid_ = [ptsOnEllipsoid_ temp(1:3) -temp(1:3)];
        end
    end
    states = xtraj.eval(ts(i));
%     ptsOnEllipsoid{i} = ptsOnEllipsoid_ + states(1:3);
    ptsOnEllipsoid_ = ptsOnEllipsoid_ + states(1:3);
    [coverEllipsoid{i},centerOfEllipsoid{i}] = MinVolEllipse(ptsOnEllipsoid_,0.001);
    
    [u,v,~] = svd(coverEllipsoid{i});
    quat = rotmat2quat(u');
    scale = diag(inv(sqrt(v)));
    rviz = [centerOfEllipsoid{i}' scale' quat'];
    dataForRVIZ = [dataForRVIZ; rviz];
end

%%
csvwrite(['funnel',num2str(k),'.csv'],dataForRVIZ);
generateCSV(xtraj.getBreaks,xtraj.pp.coefs,['refTraj',num2str(k),'.csv']);
generateCSV(utraj.getBreaks,utraj.pp.coefs,['inputTraj',num2str(k),'.csv']);
generateCSV(c.D.getBreaks,c.D.pp.coefs,['controlTraj',num2str(k),'.csv']);

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

% %%
% figure(100);clf;
% hold on
% i=N-1;
% pos = xtraj.eval(ts(i));
% surf([xx{i} xx_{i}]+pos(1),[yy{i} yy_{i}]+pos(2),[zz{i} zz_{i}]+pos(3),'FaceColor',[0.8 0.8 0.1],'FaceAlpha',0.1,'LineStyle','none');
% 
% xlabel('x')
% ylabel('y')
% zlabel('z')
% axis equal
% 
% %%
% pts = [];
% for i=1:size(xx{1},1)
%     for j=1:size(xx{1},2)
%         chad = [xx{N-1}(i,j);yy{N-1}(i,j);zz{N-1}(i,j)];
%         plot3(chad(1),chad(2),chad(3),'*')
%         pts = [pts chad]; 
%     end
% end
% 
% [AAA,ccc] = MinVolEllipse(pts,0.001);
% 
% for i=1:length(pts)
%    aaa = pts(:,i);
%    bbb = aaa'*AAA*aaa;
%    if bbb > 1
%        bbb
%    end
% end
% 
% 
% 
% %%
% for k=1:100
%     i = floor(rand(1)*size(xx{1},1));
%     j = floor(rand(1)*size(xx{1},1));
%     if i == 0
%         i = 1;
%     elseif j == 0
%         j = 1;
%     end
%     x = V0.getFrame.getPoly;
%     V_quad=.5*double(subs(diff(diff(V{1},x)',x),x,0*x))/double(rho(1));    
%     %     a = [xx{1}(i,j)*0.9 yy{1}(i,j)*0.9 zz{1}(i,j)*0.9 2 zeros(1,5)];
%     a = 10*randn(9,1);
%     r = a'*V_quad*a;
%     a = a / sqrt(r) + xtraj.eval(0);
%     sol = ode45(odefun,[ts(1) ts(end)],a,opts);
%     plot3(sol.y(1,:),sol.y(2,:),sol.y(3,:));
%     pause(0.1)
% end
% axis equal
% 
% %%
% V_in   = 0.5*double(diff(diff(V{1},x)',x))/double(rho(1));
% V_exit = 0.5*double(diff(diff(V{end},x)',x))/double(rho(end));