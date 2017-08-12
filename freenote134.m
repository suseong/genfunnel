clear
clc
% % close
% 
data = [];
breakList = [];

% for kk = 1:100
kk = 1;

zinit =  [10*(rand(1,2)-0.5) 0];
zfinal = [10*(rand(1,2)-0.5) 0];

disp([num2str(zfinal(1) - zinit(1)),' ',num2str(zinit(2)),' ',num2str(zfinal(2))])
% x-axis : [0 0 0] -> [5 0 0]
% y-axis : [-2 0 0] -> [0 1 0];
% z-axis : [0 0 0] -> [1 0 0];

% zinit = [0 0 0];
% zfinal = [5 0 0];

% zinit = [-2 0 0];
% zfinal = [0 1 0];

% zinit = [0 0 0];
% zfinal = [1 0 0];
%%
amax = 20;
g = 9.8;

a_ = 10:-0.02:0.05;
% a_ = 1.5:-0.02:0.5;
cont = [];
num_ = [];
actCont = [];
signal = [];

% jerk max : (z_min + 9.8) x omega_max / sqrt(3)
% jerk max : a_min x omega_max / sqrt(3)

for k=1:length(a_)
    chad = a_(k);
    [ztsq,zpos,zacc,ziter,zact,num,isOkay] = calc_mintime_traj(zinit,zfinal,[8 chad]);
    if isOkay == 0
        breakList{length(breakList)+1} = {zinit,zfinal};
        k = k-1;
        break;
    end
    cont(k) = ztsq(end);
    num_(k) = num;
    actCont(k) = zact;
%     if or(ztsq(end) > 30, isempty(ztsq))
    if isempty(ztsq)
        k = k-1;
        break;
    end
end

% if isOkay == 0
%     continue;
% else

if find(diff(diff(cont(1:k))) < -0.1)
    disp('bad case 1 found')
    signal1 = 1;
else
    signal1 = 0;
end

p0 = zinit(1);
pf = zfinal(1);
v0 = zinit(2);
vf = zfinal(2);

delT = (pf - p0) / (v0 + 0.5*(vf - v0));
acc = (vf - v0) / delT;

cs1 = and(pf-p0>0,and(v0>0,vf>0));
cs2 = and(pf-p0<0,and(v0<0,vf<0));

if or(cs1,cs2)
    disp('bad case 2 found')
    signal2 = 1;
else
    signal2 = 0;
end

% if and(or(signal1,signal2),abs(acc) < a_(1))
    figure(110);clf;
    subplot(3,1,1)
    hold on;
    plot(a_(1:k),cont(1:k),'.');
%     plot(abs([acc acc]),[-1 5],':');
    xlabel('acc')
    ylabel('final time')
    subplot(3,1,2)
    hold on
%     plot(a_(1:k),mod(floor(2*actCont(1:k)/10)/2,1)*2,'o');
    plot(a_(1:k),num_(1:k),'.');
    grid on
    xlabel('acc')
    ylabel('act num')
    subplot(3,1,3)
%     plot(a_(1:k),exp(1./cont(1:k)));    
    plot(a_(1:k),actCont,'*');    
    pause(0.1);
% end

% end
% end

%%
% zinit
% zfinal

% figure(110);clf;
% subplot(3,1,1)
% plot(a_(1:k),cont(1:k),'.');
% xlabel('acc')
% ylabel('final time')
% subplot(3,1,2)
% plot(a_(1:k-1),diff(cont(1:k)),'.');
% xlabel('acc')
% ylabel('final time diff')
% subplot(3,1,3)
% plot(a_(1:k-2),diff(diff(cont(1:k))),'.');
% xlabel('acc')
% ylabel('final time diff')



% figure(104);clf;
% subplot(4,1,1)
% plot(a_(1:k),exp(1./cont(1:k)),'.k')
% xlabel('acc')
% ylabel('final time')
% hold on
% % plot([a_(1:200:2200) a_(length(cont))],exp(1./[cont(1:400:2400) cont(end)]),'*');
% plot(a_(floor(linspace(1,k,11))),exp(1./cont(floor(linspace(1,k,11)))),'*');

% subplot(4,1,2)
% plot(a_(1:k-1),actCont(1:k-1),'*')

% pft = polyfit(a_(floor(linspace(1,k,10))),exp(1./cont(floor(linspace(1,k,10)))),6);
% mu
% pft = polyfit(a_(1:400:2401),exp(1./cont(1:400:2401)),5);
% pft = polyfit([a_(1:200:2400) a_(length(cont))],exp(1./[cont(1:200:2400) cont(end)]),5);

% cs = pft;
% xx = a_(k-1):0.1:20;
% xx = a_;
% chad = 1./log((cs(1).*xx.^6+cs(2).*xx.^5+cs(3).*xx.^4+cs(4).*xx.^3+cs(5).*xx.^2+cs(6).*xx.^1+cs(7).*xx.^0));
% mm = max(abs(cont - chad(1:k)));

% if (mm > 2)
%     signal = 0;
% else
%     signal = 1;
% end

% data{kk} = {zinit,zfinal,signal};
% kk
% end
% subplot(4,1,2)
% hold on
% plot(a_(1:k),cont,'.k')
% plot(xx(1:k),chad(1:k),'g');
% % plot(xx,1./log((cs(1).*xx.^5+cs(2).*xx.^4+cs(3).*xx.^3+cs(4).*xx.^2+cs(5).*xx.^1+cs(6).*xx.^0)),'g')
% % plot(xx,1./log((cs(1).*xx.^6+cs(2).*xx.^5+cs(3).*xx.^4+cs(4).*xx.^3+cs(5).*xx.^2+cs(6).*xx.^1+cs(7).*xx.^0)),'g')
% % % plot(xx,1./(cs(10).*xx.^9+cs(9).*xx.^8+cs(8).*xx.^7+cs(7).*xx.^6+cs(6).*xx.^5+cs(5).*xx.^4+cs(4).*xx.^3+cs(3).*xx.^2+cs(2).*xx.^1+cs(1).*xx.^0),'g')
% xlabel('acc')
% ylabel('final time')
% 
% subplot(4,1,3)
% plot(xx(1:k),chad(1:k) - cont,'k');
% xlabel('acc')
% ylabel('diff')
% 
% subplot(4,1,4)
% plot(a_(1:k),num_,'.k')

% end

%%
% figure(100);clf;
% hold on
% for k = 1:length(data)
%     if size(data{k,1}) > 0
%         if data{k,3} == 0
%             %         plot3(data{k}{1}(1) - data{k}{2}(1),data{k}{1}(2),data{k}{2}(2),'o')
%         else
%             plot3(data{k,1}(1) - data{k,2}(1),data{k,1}(2),data{k,2}(2),'.','markersize',30)
%         end
%     end
% end
% % 
% xlabel('pos diff')
% ylabel('init vel')
% zlabel('final vel')
% 
% figure(101);clf;
% hold on
% for k = 1:length(data)
%     if size(data{k,1}) > 0
%         if data{k,3} == 0
%               plot3(data{k,1}(1) - data{k,2}(1),data{k,1}(2),data{k,2}(2),'o')
%         else
% %             plot3(data{k,1}(1) - data{k,2}(1),data{k,1}(2),data{k,2}(2),'.','markersize',30)
%         end
%     end
% end
% % 
% xlabel('pos diff')
% ylabel('init vel')
% zlabel('final vel')


