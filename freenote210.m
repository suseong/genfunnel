
% k=11
% PPP_{k} = PPP{1};
% RHO_{k} = RHO{1};
% SSS_{k} = SSS{1};

%%
% clf;
% for j = 1:11
%     figure(1);

%%

figure(8)
subplot(1,3,1)
hold on
ang = -pi:0.2:pi;
p = SSS_{6}(:,100);
P = [p(1)   0    0  p(3)   0    0;
    0  p(1)   0    0  p(3)   0;
    0    0  p(2)   0    0  p(4);
    p(3)   0    0  p(5)   0    0;
    0  p(3)   0    0  p(5)   0;
    0    0  p(4)   0    0  p(6)];

kk = 1;
p1 = [P(kk,kk) P(kk,kk+3);P(kk+3,kk) P(kk+3,kk+3)];
invp1 = inv(sqrtm(p1));
for k=1:length(ang)
    xx = invp1*[cos(ang(k));sin(ang(k))];
    plot(xx(1),xx(2),'.','markersize',15)
end

%%
% for j=1:11
% p = SSS_{5}(:,100);
% P = [p(1)   0    0  p(3)   0    0;
%     0  p(1)   0    0  p(3)   0;
%     0    0  p(2)   0    0  p(4);
%     p(3)   0    0  p(5)   0    0;
%     0  p(3)   0    0  p(5)   0;
%     0    0  p(4)   0    0  p(6)];
% exy(j) = sqrt(1/p(1));
% ez(j) =  sqrt(1/p(2));
% end
% 
% figure(1);clf
% subplot(3,1,1)
% plot(exy,'-*')
% subplot(3,1,2)
% plot(ez,'-*')
% subplot(3,1,3)
% hold on
% plot(exy,'-*r')
% plot(ez,'-*b')
%     
%     for kk = 2:3
%         p1 = [P(kk,kk) P(kk,kk+3);P(kk+3,kk) P(kk+3,kk+3)];
%         invp = inv(sqrtm(p1));
%         
%         subplot(1,2,kk-1)
%         axis equal
%         hold on
%         for k=1:length(ang)
%             xx = invp*[cos(ang(k));sin(ang(k))];
%             plot(xx(1),xx(2),'.','markersize',5)
%         end
%         axis([-0.5 0.5 -2 2]*4);
%         grid on
%     end
% %     keyboard
% end

%%
ang = -pi:0.3:pi;

xx = []; yy = []; zz = [];
xx_ = []; yy_ = []; zz_ = [];

for acc = 1:11
    for num = 1:50
        p = SSS_{acc}(:,(num-1)*2+1);
        P = [p(1)   0    0  p(3)   0    0;
            0  p(1)   0    0  p(3)   0;
            0    0  p(2)   0    0  p(4);
            p(3)   0    0  p(5)   0    0;
            0  p(3)   0    0  p(5)   0;
            0    0  p(4)   0    0  p(6)];
        
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
figure(200);clf
for acc = 1:11
    subplot(1,11,acc);
    for num = 1:50
        surf([xx{acc,num} xx_{acc,num}],[yy{acc,num} yy_{acc,num}],[zz{acc,num} zz_{acc,num}]+num*0.2,'FaceColor',[0.8 0.8 0.8],'FaceAlpha',0.8,'LineStyle','none');
        hold on
        axis([-1 1 -1 1 -2 12]*1.2)
        axis equal
        view(20,30)
%         pause(0.1)
    end
end

%%
% syms ex ey ez evx evy evz real
% syms epbar edbar real
% e = [ex ey ez evx evy evz]';

p = SSS_{5}(:,100);
P = [p(1)   0    0  p(3)   0    0;
    0  p(1)   0    0  p(3)   0;
    0    0  p(2)   0    0  p(4);
    p(3)   0    0  p(5)   0    0;
    0  p(3)   0    0  p(5)   0;
    0    0  p(4)   0    0  p(6)];
P = P*RHO_{3}(100);

p_ = SSS_{3}(:,101);
P_ = [p_(1)   0    0  p_(3)    0     0;
        0  p_(1)   0     0  p_(3)    0;
        0     0  p_(2)   0     0  p_(4);
      p_(3)   0    0  p_(5)    0     0;
        0  p_(3)   0     0  p_(5)    0;
        0     0  p_(4)   0     0  p_(6)];
P_ = P_*RHO_{3}(101);

ep = [sqrt(1/p(1)) sqrt(1/p(1)) sqrt(1/p(2))]';
ed = [sqrt(1/p(5)) sqrt(1/p(5)) sqrt(1/p(6))]';

epbar = norm(0);
edbar = norm(ed);
e = [0;0;0;ed];

rho = RHO_{3}(100);
Kp = diag([10 10 15]);
Kd = diag([4 4 5]);
maxKp = 15;
maxKd = 5;
A = [zeros(3,3) eye(3); -Kp -Kd];
unc = 1.0;
Er = 0.05;
ar = 15;
Ppv = max(max(P(1:3,4:6)));
Pv = max(max(P(4:6,4:6)));

Vdot = e'*(P*A+A'*P)*e ...
       + 2*(unc + Er*(maxKp*epbar + maxKd*edbar + ar))*(Ppv*epbar + Pv*edbar) ...
       + e'*(P_-P)/0.05*e;

%%
% figure(3);clf
% acc = 4;
% for num = 1:50
%     surf([xx{acc,num} xx_{acc,num}],[yy{acc,num} yy_{acc,num}],[zz{acc,num} zz_{acc,num}]+num*0.2,'FaceColor',[0.8 0.8 0.8],'FaceAlpha',0.8,'LineStyle','none');
%     hold on
%     axis([-1 1 -1 1 -2 12])
%     axis equal
%     view(20,30)
% end
