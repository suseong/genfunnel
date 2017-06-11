close all
clc

%%
ang = -pi:0.1:pi+0.1;

figure(1);clf;
for jj = 1:29
    hold on
    P = reshape(double(shells{1}(:,jj*5)),6,6);
    kk = 3;
    p1 = [P(kk,kk) P(kk,kk+3);P(kk+3,kk) P(kk+3,kk+3)];
    invp1 = inv(sqrtm(p1));
    for k=1:length(ang)
       x(:,k) = invp1*[cos(ang(k));sin(ang(k))];
    end
    plot(x(2,:),x(1,:))
%    axis([-0.5 0.5 -1.5 1.5]);
   axis equal
%    pause(0.1);
end

%%
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
figure(2);clf
for acc = 1:length(shells)
    subplot(1,8,acc);
    for num = 1:29
        surf([xx{acc,num} xx_{acc,num}],[yy{acc,num} yy_{acc,num}],[zz{acc,num} zz_{acc,num}]+num*0.05,'FaceColor',[0.8 0.8 0.8],'FaceAlpha',0.03*num,'LineStyle','none');
        hold on
        axis equal
%         pause(0.1)
    end
end

%%
for k = 1:130
   P = reshape(double(shells{1}(:,k)),6,6);
   P_ = reshape(double(shells{1}(:,k+10)),6,6);
   
   chad = P - P_;
   diff(:,k) = eig(chad);
end

%%
xs = [0;0;0;2;0;0];
xf = [0;1.8;0;-2;0;0];

% xs = [0;1.8;0;-2;0;0];
% xf = [-1;1.8;0;0;0;0];

tf = 1;
delT = tf / 80;
ts = linspace(0,tf,1+80);

[xtraj,utraj] = refTraj2(xs,xf,tf,delT);

states = [];
for i=1:length(ts)
   chad = xtraj.deriv(ts(i));
   states(:,i) = xtraj.eval(ts(i));
   ar(i) = norm(chad(4:6));
end

%%
figure(32);clf;
for i=1:9
    subplot(10,1,i)
    plot(ts,states(i,:))
end
subplot(10,1,10)
plot(ts,ar)

%%
funnel = [];
validFinal = 140;
funnel{1} = reshape(double(shells{1}(:,validFinal)),6,6);
prevFunnelIdx = 1;
prevTimeIdx = validFinal;
infos = [1;validFinal];

% funnel{1} = reshape(double(shells{1}(:,72)),6,6);
% prevFunnelIdx = 1;
% prevTimeIdx = 72;
% infos = [1;72];

for k = 2 : 81
    ar_chad = 0.98 + ar(k);
    funnelIdx = floor(ar_chad / 2) + 1;
    if funnelIdx == prevFunnelIdx
       timeIdx = prevTimeIdx + 1;
       if timeIdx > validFinal
          timeIdx = validFinal;
       end
       funnel{k} = reshape(double(shells{funnelIdx}(:,timeIdx)),6,6); 
       prevTimeIdx = timeIdx;
    else
        cnt = 0;
        while(1)
            Pchad = reshape(double(shells{funnelIdx}(:,validFinal - cnt)),6,6);
            Pdiff = Pchad - funnel{k-1};
            eigchad = eig(Pdiff);
            if(sum(find(eigchad < 0)) == 21)
                break;
            end
            cnt = cnt + 1;
        end
        funnel{k} = reshape(double(shells{funnelIdx}(:,validFinal - cnt)),6,6);   
        prevTimeIdx = validFinal - cnt;
    end
    
    prevFunnelIdx = funnelIdx;
    infos(:,k) = [funnelIdx;prevTimeIdx];
end

%%
figure(23);
clf;
for k = 1:81
    acc = infos(1,k); num = floor(infos(2,k)/5)+1;
    surf([xx{acc,num} xx_{acc,num}]+states(1,k),[yy{acc,num} yy_{acc,num}]+states(2,k),[zz{acc,num} zz_{acc,num}]+states(3,k),'FaceColor',[0.8 0.8 0.8],'FaceAlpha',0.1,'LineStyle','none');
    hold on
    axis equal
end

%%






