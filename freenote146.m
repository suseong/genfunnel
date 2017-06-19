clear all
close all
clc

load funnel_3deg.mat;

%%
figure(10)
for i=1:length(Rho)
    subplot(8,1,i)
    plot(Rho{i})
end

validFinal = 180;

%%
P = [];
for i=1:length(Lyap)
    for j=1:validFinal
        p = Lyap{i}(6*(j-1)+1:6*j);
        P_ = [p(1)   0    0  p(3)   0    0;
               0  p(1)   0    0  p(3)   0;
               0    0  p(2)   0    0  p(4);
             p(3)   0    0  p(5)   0    0;
               0  p(3)   0    0  p(5)   0;
               0    0  p(4)   0    0  p(6)];
        P{i,j} = P_;
    end
end

%%
ang = -pi:0.2:pi+0.1;
xx = []; yy = []; zz = [];
xx_ = []; yy_ = []; zz_ = [];

chad = 5;
for acc = 1:length(Lyap)
    for num = 1:floor(validFinal / chad)
        p = P{acc,num*chad}/ Rho{i}(j);
        invp = inv(p);      
        for j=1:length(ang)
            for k=1:length(ang)
                a = [sin(ang(j))*sin(ang(k));sin(ang(j))*cos(ang(k));cos(ang(j));zeros(3,1)];
                lambda = sqrt(a'*invp*a);
                temp = 1/lambda*invp*a;
                xx{acc,num}(j,k) = temp(1); yy{acc,num}(j,k) = temp(2); zz{acc,num}(j,k) = temp(3);
%                 xx_{acc,num}(j,k) = -temp(1); yy_{acc,num}(j,k) = -temp(2); zz_{acc,num}(j,k) = -temp(3);
            end
        end
    end
end
p = [];

%%
figure(20);clf;
for acc = 1:length(Lyap)
    subplot(1,8,acc);
    hold on
    for num = 1:floor(validFinal / chad)
%         surf([xx{acc,num} xx_{acc,num}],[yy{acc,num} yy_{acc,num}],[zz{acc,num} zz_{acc,num}]+num*0.05,'FaceColor',[0.8 0.8 0.8],'FaceAlpha',0.01*num,'LineStyle','none');
        xlabel('x')
        ylabel('y')
        zlabel('time')
        surf(xx{acc,num},yy{acc,num},zz{acc,num}+num*0.05,'FaceColor',[0.8 0.8 0.8],'FaceAlpha',0.01*num,'LineStyle','none');
        axis equal
        view(30,20)
    end
end

%%
T  =  [0  0.7  1.0  1.3  2.0]*2.0;

p(1,:) = [0 2  4 6 0];
p(2,:) = [0 1 -1 0 0];
p(3,:) = [1.8 1.6 1.4 1.2 1.8];

% T  =  [0  0.7  1.0  1.3  2.0]*2;
% 
% p(1,:) = [0 3 4 3 0];
% p(2,:) = [1 -1 0 1 -1];
% p(3,:) = [1.8 1.6 1.4 1.4 1.4];

opts = qpOASES_options('mpc');
tic;
[cost,T,coeffs_x,coeffs_y,coeffs_z] = calc_min_cost4(T,p,opts);
toc;

n = 8;
i=1; a_c_x = coeffs_x(n*(i-1)+1:n*i);
i=2; b_c_x = coeffs_x(n*(i-1)+1:n*i);
i=3; c_c_x = coeffs_x(n*(i-1)+1:n*i);
i=4; d_c_x = coeffs_x(n*(i-1)+1:n*i);

i=1; a_c_y = coeffs_y(n*(i-1)+1:n*i);
i=2; b_c_y = coeffs_y(n*(i-1)+1:n*i);
i=3; c_c_y = coeffs_y(n*(i-1)+1:n*i);
i=4; d_c_y = coeffs_y(n*(i-1)+1:n*i);

i=1; a_c_z = coeffs_z(n*(i-1)+1:n*i);
i=2; b_c_z = coeffs_z(n*(i-1)+1:n*i);
i=3; c_c_z = coeffs_z(n*(i-1)+1:n*i);
i=4; d_c_z = coeffs_z(n*(i-1)+1:n*i);

Ta = T(1):0.01:T(2);
Tb = T(2):0.01:T(3);
Tc = T(3):0.01:T(4);
Td = T(4):0.01:T(5);

p1_x = calc_trj_(a_c_x,Ta);
p2_x = calc_trj_(b_c_x,Tb);
p3_x = calc_trj_(c_c_x,Tc);
p4_x = calc_trj_(d_c_x,Td);

p1_y = calc_trj_(a_c_y,Ta);
p2_y = calc_trj_(b_c_y,Tb);
p3_y = calc_trj_(c_c_y,Tc);
p4_y = calc_trj_(d_c_y,Td);

p1_z = calc_trj_(a_c_z,Ta);
p2_z = calc_trj_(b_c_z,Tb);
p3_z = calc_trj_(c_c_z,Tc);
p4_z = calc_trj_(d_c_z,Td);

figure(30);clf;
% subplot(2,1,1)
hold on
plot3([p1_x p2_x p3_x p4_x],[p1_y p2_y p3_y p4_y],[p1_z p2_z p3_z p4_z],'-k');
for i=1:5
    plot3(p(1,i),p(2,i),p(3,i),'*','markersize',10)
end
axis equal;
view(30,20);
xlabel('x');
ylabel('y');
zlabel('z');

acc1_x = calc_a_(a_c_x,Ta);
acc2_x = calc_a_(b_c_x,Tb);
acc3_x = calc_a_(c_c_x,Tc);
acc4_x = calc_a_(d_c_x,Td);

acc1_y = calc_a_(a_c_y,Ta);
acc2_y = calc_a_(b_c_y,Tb);
acc3_y = calc_a_(c_c_y,Tc);
acc4_y = calc_a_(d_c_y,Td);

acc1_z = calc_a_(a_c_z,Ta);
acc2_z = calc_a_(b_c_z,Tb);
acc3_z = calc_a_(c_c_z,Tc);
acc4_z = calc_a_(d_c_z,Td);

pos = [p1_x(1:end-1) p2_x(1:end-1) p3_x(1:end-1) p4_x;
       p1_y(1:end-1) p2_y(1:end-1) p3_y(1:end-1) p4_y;
       p1_z(1:end-1) p2_z(1:end-1) p3_z(1:end-1) p4_z];


acc = [acc1_x(1:end-1) acc2_x(1:end-1) acc3_x(1:end-1) acc4_x;
       acc1_y(1:end-1) acc2_y(1:end-1) acc3_y(1:end-1) acc4_y;
       acc1_z(1:end-1) acc2_z(1:end-1) acc3_z(1:end-1) acc4_z];

acc_norm = 0;
for k=1:length(acc)
   acc_norm(k) = norm(acc(:,k)); 
end
% subplot(2,1,2)
% plot([Ta(1:end-1) Tb(1:end-1) Tc(1:end-1) Td],acc_norm);
% xlabel('sec');
% ylabel('acc');

%%
funnel = [];
funnel{1} = P{1,validFinal};
prevFunnelIdx = 1;
prevTimeIdx = validFinal;
infos = [1;validFinal];
tt = [Ta(1:end-1) Tb(1:end-1) Tc(1:end-1) Td];

for k = 2 : length(tt)
    ar_chad = acc_norm(k);
    funnelIdx = floor(ar_chad / 2) + 1;
    if funnelIdx == prevFunnelIdx
       timeIdx = prevTimeIdx + 1;
       if timeIdx > validFinal
          timeIdx = validFinal;
       end
       funnel{k} = P{funnelIdx,timeIdx}; 
       prevTimeIdx = timeIdx;
    else
        cnt = 0;
        while(1)
            Pchad = P{funnelIdx,validFinal - cnt};
            Pdiff = Pchad - funnel{k-1};
            eigchad = eig(Pdiff);
            if(sum(find(eigchad < 0)) == 21)
                break;
            end
            cnt = cnt + 1;
        end
        funnel{k} = P{funnelIdx,validFinal - cnt};   
        prevTimeIdx = validFinal - cnt;
    end
    
    prevFunnelIdx = funnelIdx;
    infos(:,k) = [funnelIdx;prevTimeIdx];
end

%%
figure(23);
clf;
for k = 1:length(tt)
    acc = infos(1,k); num = floor(infos(2,k)/5);
%     surf([xx{acc,num} xx_{acc,num}]+pos(1,k),[yy{acc,num} yy_{acc,num}]+pos(2,k),[zz{acc,num} zz_{acc,num}]+pos(3,k),'FaceColor',[0.8 0.8 0.8],'FaceAlpha',0.1,'LineStyle','none');
    surf([xx{acc,num}]+pos(1,k),[yy{acc,num}]+pos(2,k),[zz{acc,num}]+pos(3,k),'FaceColor',[0.8 0.8 0.8],'FaceAlpha',0.1,'LineStyle','none');
    hold on
end
grid off
axis equal
xlabel('x')
ylabel('y')
zlabel('z')

%%
bag = rosbag('pos2.bag');
pose = bag.timeseries.Data;
simPos = pose(:,1:3);

%%
plot3(simPos(:,1),simPos(:,2),simPos(:,3),'linewidth',3);





