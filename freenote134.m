clear
% close

zinit =  [5*randn(1,2) 0];
zfinal = [5*randn(1,2) 0];

amax = 20;
g = 9.8;

a_ = 20:-0.05:0.1;

for k=1:length(a_)
    chad = a_(k);
    [ztsq,zpos,zacc,ziter,zact] = calc_mintime_traj(zinit,zfinal,[30 chad]);
    cont(k) = ztsq(end);
    actCont(k) = zact;
    if(ztsq(end) > 40)
        break;
    end
end

%%
zinit
zfinal

figure(104);clf;
subplot(3,1,1)
plot(a_(1:k-1),exp(1./cont(1:k-1)),'.k')
xlabel('acc')
ylabel('final time')
hold on
% plot([a_(1:200:2200) a_(length(cont))],exp(1./[cont(1:400:2400) cont(end)]),'*');
plot(a_(floor(linspace(1,k-1,10))),exp(1./cont(floor(linspace(1,k-1,10)))),'*');

subplot(3,1,2)
plot(a_(1:k-1),actCont(1:k-1),'*')


pft = polyfit(a_(floor(linspace(1,k-1,10))),exp(1./cont(floor(linspace(1,k-1,10)))),6);
% pft = polyfit(a_(1:400:2401),exp(1./cont(1:400:2401)),5);
% pft = polyfit([a_(1:200:2400) a_(length(cont))],exp(1./[cont(1:200:2400) cont(end)]),5);

cs = pft;
xx = a_(k-1):0.1:20;

subplot(3,1,3)
hold on
plot(a_(1:k-1),cont(1:k-1),'.k')
% plot(xx,1./log((cs(1).*xx.^5+cs(2).*xx.^4+cs(3).*xx.^3+cs(4).*xx.^2+cs(5).*xx.^1+cs(6).*xx.^0)),'g')
% plot(xx,1./log((cs(1).*xx.^6+cs(2).*xx.^5+cs(3).*xx.^4+cs(4).*xx.^3+cs(5).*xx.^2+cs(6).*xx.^1+cs(7).*xx.^0)),'g')
% % plot(xx,1./(cs(10).*xx.^9+cs(9).*xx.^8+cs(8).*xx.^7+cs(7).*xx.^6+cs(6).*xx.^5+cs(5).*xx.^4+cs(4).*xx.^3+cs(3).*xx.^2+cs(2).*xx.^1+cs(1).*xx.^0),'g')

xlabel('acc')
ylabel('final time')



