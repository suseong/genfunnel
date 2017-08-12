clear
clc

zinit =  [10*(rand(1,2)-0.5) 10];
zfinal = [10*(rand(1,2)-0.5) -10];

amax = 20;
g = 9.8;

a_max = sqrt(amax^2 - g^2);
% a_ = flip(linspace(0.2,a_max,20));
a_ = flip([0.2 0.3 0.4 0.5 0.7 1.0 2.0 4.0 6.0 9.0 12.0 15.0 18.0]);

for k=1:length(a_)
    chad = a_(k);
    [ztsq,zpos,zacc,ziter,zact,num,isOkay] = calc_mintime_traj(zinit,zfinal,[20 chad]);
    if isempty(ztsq)
        k = k-1;
        break;
    end
    cont(k) = ztsq(end);
    num_(k) = num;
    actCont(k) = zact;
end
a_ = a_(1:k);

figure(110);clf;
subplot(3,1,1)
hold on
plot(a_(1:k),cont,'o');
subplot(3,1,2)
hold on
plot(a_(1:k),num_,'o');
axis([0 a_(1) -2 2])

%%
idx = 0;
for k = 1:length(cont)-1
    if num_(k)*num_(k+1) < 0
        idx = k;

        maxAcc = a_(idx);
        maxTf = cont(idx);
        maxNum = num_(idx);
        
        minAcc = a_(idx+1);
        minTf = cont(idx+1);
        minNum = num_(idx+1);
        
        break;
    end
end

% bisection to find the discontinuous point

if idx > 0
    acc_ = (minAcc + maxAcc)/2;
    
    for k=1:5
        [ztsq,zpos,zacc,ziter,zact,num,isOkay] = calc_mintime_traj(zinit,zfinal,[20 acc_]);
        
        if abs(ztsq(end) - maxTf) < abs(ztsq(end) - minTf)
%         if num == num_(idx)
           maxAcc = acc_;
           maxTf = ztsq(end);
           maxNum = num_(idx);
           maxActNum = zact;
        else
           minAcc = acc_;
           minTf = ztsq(end);
           minNum = num_(idx+1);
           minActNum = zact;
        end

        cont = [cont minTf maxTf];
        num_ = [num_ minNum maxNum];
        a_ = [a_ minAcc maxAcc];
        
        subplot(3,1,1)
        plot(acc_,ztsq(end),'*')    
        subplot(3,1,2)
        plot(acc_,num,'*')
        
        acc_ = (minAcc + maxAcc)/2;
    end

    cont = [cont minTf maxTf];
    num_ = [num_ minNum maxNum];
    a_ = [a_ minAcc maxAcc];
    actCont = [actCont minActNum maxActNum];
    
    subplot(3,1,1)
    plot(a_(end),cont(end),'ksq','markersize',13)
    plot(a_(end-1),cont(end-1),'rsq','markersize',13)
    subplot(3,1,2)
    plot(a_(end),num_(end),'ksq','markersize',13)
    plot(a_(end-1),num_(end-1),'rsq','markersize',13)
end

% subplot(3,1,3)
% plot(a_,actCont,'o')

%%
bigAccIdx = find(num_ == num_(1));
smallAccIdx = find(num_ == -1*num_(1));

if ~isempty(bigAccIdx)
    accMat = [];
    pft = polyfit(a_(bigAccIdx),exp(1./cont(bigAccIdx)),min(6,length(bigAccIdx)-1));
    accTemp = linspace(min(a_(bigAccIdx)),max(a_(bigAccIdx)),100);
    for k = 1:length(pft)
        accMat(k,:) = accTemp.^(length(pft)-k);
    end
    tfTemp = pft*accMat;
    subplot(3,1,1)
    plot(accTemp,1./log(tfTemp),'-g');
    
end

if ~isempty(smallAccIdx)
    accMat = [];
    pft = polyfit(a_(smallAccIdx),exp(1./cont(smallAccIdx)),min(6,length(smallAccIdx)-1));
    accTemp = linspace(min(a_(smallAccIdx)),max(a_(smallAccIdx)),100);
    for k = 1:length(pft)
        accMat(k,:) = accTemp.^(length(pft)-k);
    end
    tfTemp = pft*accMat;
    subplot(3,1,1)
    plot(accTemp,1./log(tfTemp),'-r');    
end










