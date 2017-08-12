clear;
clf;

% % x-axis : [0 0 0] -> [5 0 0]
% % y-axis : [-2 0 0] -> [0 1 0];
% % z-axis : [0 0 0] -> [1 0 0];
% 
% zinit{1} = [0 0 0];
% zfinal{1} = [5 0 0];
% 
% zinit{2} = [-2 0 0];
% zfinal{2} = [0 1 0];
% 
% zinit{3} = [0 0 0];
% zfinal{3} = [1 0 0];
% 
% % acc_x_1 = 11.4 m/s^2 (1.752 sec), 
% % acc_y_1 = 1.8 m/s^2 (1.763 sec),
% % acc_z_1 = 1.4 m/s^2 (1.763 sec),

% x-axis : [5 0 0] -> [0 0 0]
% y-axis : [0 1 0] -> [2 0 0];
% z-axis : [1 0 0] -> [0 0 0];

zinit{1} = [5 0 0];
zfinal{1} = [0 0 0];

zinit{2} = [0 1 0];
zfinal{2} = [2 0 0];

zinit{3} = [1 0 0];
zfinal{3} = [0 0 0];

% acc_x_2 = 11.4 m/s^2 (1.752 sec), 
% acc_y_2 = 1.8 m/s^2 (1.763 sec),
% acc_z_2 = 1.4 m/s^2 (1.763 sec),

for j=1:3

a_ = 15:-0.2:0.4;
cont = [];
num_ = [];
actCont = [];
signal = [];

% jerk max : (z_min + 9.8) x omega_max / sqrt(3)
% jerk max : a_min x omega_max / sqrt(3)

for k=1:length(a_)
    chad = a_(k);
    if j == 1
%         zfinal_ = zfinal{j};
%         zfinal_(3) = max(-abs(chad),-9);
%         [ztsq,zpos,zacc,ziter,zact,num,isOkay] = calc_mintime_traj(zinit{j},zfinal_,[20 chad]);
        zinit_ = zinit{j};
        zinit_(3) = max(-abs(chad),-9);
        [ztsq,zpos,zacc,ziter,zact,num,isOkay] = calc_mintime_traj(zinit_,zfinal{j},[20 chad]);
    else
        [ztsq,zpos,zacc,ziter,zact,num,isOkay] = calc_mintime_traj(zinit{j},zfinal{j},[20 chad]);
    end
    cont(k) = ztsq(end);
    num_(k) = num;
    actCont(k) = zact;
    if(ztsq(end) > 40)
        break;
    end
end

pft = polyfit(a_(floor(linspace(1,k,10))),exp(1./cont(floor(linspace(1,k,10)))),6);
cs = pft; coeffs{j} = pft;
xx = a_;
chad = 1./log((cs(1).*xx.^6+cs(2).*xx.^5+cs(3).*xx.^4+cs(4).*xx.^3+cs(5).*xx.^2+cs(6).*xx.^1+cs(7).*xx.^0));
% chad = cs(1).*xx.^6+cs(2).*xx.^5+cs(3).*xx.^4+cs(4).*xx.^3+cs(5).*xx.^2+cs(6).*xx.^1+cs(7).*xx.^0;

figure(110);
subplot(3,3,3*j-2)
hold on;
plot(a_(1:k),cont(1:k),'.');
xlabel('acc (m/s^2)')
ylabel('final time')

subplot(3,3,3*j-1)
hold on
plot(a_,chad,'-k');
plot(a_(1:k),cont(1:k),'o');
xlabel('acc (m/s^2)')
ylabel('final time')
legend('regressed','computed')

subplot(3,3,3*j)
hold on
plot(a_(1:k),floor(actCont(1:k)/5)/2,'o');
plot(a_(1:k),num_,'.');
axis([0 10 0 5])
xlabel('acc (m/s^2)')
ylabel('profile')


end

