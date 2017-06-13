clear all
close all
clc

T  =  [0  0.7  1.0  1.3  2.0]*2.5;

% p(1,:) = [0 3 4 3 0];
% p(2,:) = [1 -1 0 1 -1];
% p(3,:) = [1.8 1.6 1.4 1.4 1.4];

p(1,:) = [0 2  4 6 0];
p(2,:) = [0 1 -1 0 0];
p(3,:) = [1.8 1.6 1.4 1.2 1.8];

opts = qpOASES_options('mpc');
tic;
[cost,T,coeffs_x,coeffs_y,coeffs_z] = calc_min_cost4(T,p,opts);
toc;

%%
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

%%
figure(1)
hold on
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

plot3([p1_x p2_x p3_x p4_x],[p1_y p2_y p3_y p4_y],[p1_z p2_z p3_z p4_z],'-k');

for i=1:5
    plot3(p(1,i),p(2,i),p(3,i),'*','markersize',10)
end

grid on
axis equal

time = T(2:end);
c_x = [a_c_x b_c_x c_c_x d_c_x];
c_y = [a_c_y b_c_y c_c_y d_c_y];
c_z = [a_c_z b_c_z c_c_z d_c_z];

genCSV(time,c_x,c_y,c_z,'coeffs.csv')

%%
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

%%
figure(2)

subplot(2,3,1)
plot([Ta Tb Tc Td],[p1_x p2_x p3_x p4_x])
subplot(2,3,2)
plot([Ta Tb Tc Td],[p1_y p2_y p3_y p4_y])
subplot(2,3,3)
plot([Ta Tb Tc Td],[p1_z p2_z p3_z p4_z])

subplot(2,3,4)
plot([Ta Tb Tc Td],[acc1_x acc2_x acc3_x acc4_x])
subplot(2,3,5)
plot([Ta Tb Tc Td],[acc1_y acc2_y acc3_y acc4_y])
subplot(2,3,6)
plot([Ta Tb Tc Td],[acc1_z acc2_z acc3_z acc4_z])

%%
acc = [acc1_x acc2_x acc3_x acc4_x;acc1_y acc2_y acc3_y acc4_y;acc1_z acc2_z acc3_z acc4_z];

for k=1:length(acc)
   acc_norm(k) = norm(acc(:,k)); 
end

figure(3)
plot([Ta Tb Tc Td],acc_norm);

%%
load('shells_5deg.mat')

%%
funnel = [];
validFinal = 140;
funnel{1} = reshape(double(shells{1}(:,validFinal)),6,6);
prevFunnelIdx = 1;
prevTimeIdx = validFinal;
infos = [1;validFinal];

for k = 2 : length([Ta Tb Tc Td])
    ar_chad = 0.98 + acc_norm(k);
    funnelIdx = floor(ar_chad / 2) + 1;
    if funnelIdx == prevFunnelIdx
       timeIdx = prevTimeIdx + 1;
       if timeIdx > validFinal
          timeIdx = validFinal;
       end
       funnel{k} = reshape(shells{funnelIdx}(:,timeIdx),6,6); 
       prevTimeIdx = timeIdx;
    else
        cnt = 0;
        while(1)
            Pchad = reshape(shells{funnelIdx}(:,validFinal - cnt),6,6);
            Pdiff = Pchad - funnel{k-1};
            eigchad = eig(Pdiff);
            if(sum(find(eigchad < 0)) == 21)
                break;
            end
            cnt = cnt + 1;
        end
        funnel{k} = reshape(shells{funnelIdx}(:,validFinal - cnt),6,6);   
        prevTimeIdx = validFinal - cnt;
    end
    
    prevFunnelIdx = funnelIdx;
    infos(:,k) = [funnelIdx;prevTimeIdx];
end

%%
ang = -pi:0.1:pi+0.1;

xx = []; yy = []; zz = [];
xx_ = []; yy_ = []; zz_ = [];

for acc = 1:length(shells)
    for num = 1:29
        P = reshape(double(shells{acc}(:,num*5)),6,6);
        invP = inv(P);      
        for j=1:length(ang)
            for k=1:length(ang)
                a = [sin(ang(j))*sin(ang(k));sin(ang(j))*cos(ang(k));cos(ang(j));zeros(3,1)];
                lambda = sqrt(a'*invP*a);
                temp = 1/lambda*invP*a;
                xx{acc,num}(j,k) = temp(1); yy{acc,num}(j,k) = temp(2); zz{acc,num}(j,k) = temp(3);
                xx_{acc,num}(j,k) = -temp(1); yy_{acc,num}(j,k) = -temp(2); zz_{acc,num}(j,k) = -temp(3);
            end
        end
    end
end

%%
px = [p1_x p2_x p3_x p4_x];
py = [p1_y p2_y p3_y p4_y];
pz = [p1_z p2_z p3_z p4_z];

figure(25);
clf;
for k = 1:length([Ta Tb Tc Td])
    acc = infos(1,k); num = floor(infos(2,k)/5)+1;
%     surf([xx{acc,num} xx_{acc,num}]+px(k),[yy{acc,num} yy_{acc,num}]+py(k),[zz{acc,num} zz_{acc,num}]+pz(k),'FaceColor',[0.8 0.8 0.8],'FaceAlpha',0.1,'LineStyle','none');
    surf([xx{acc,num}]+px(k),[yy{acc,num}]+py(k),[zz{acc,num}]+pz(k),'FaceColor',[0.8 0.8 0.8],'FaceAlpha',0.1,'LineStyle','none');
    hold on
    axis equal
end
view(0,90)



















