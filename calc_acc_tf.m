function [ACC,TF,possibleRange] = calc_acc_tf(init,final,figNum)

a_ = flip([0.02   0.1   0.2   0.3   0.4   0.5   0.6   0.8   1.0   1.5   2.0   3.0   4.0   5.0   7.0   10.0]);

for k=1:length(a_)
    chad = a_(k);
    [tsq,~,~,~,~,num,~] = calc_mintime_traj_(init,final,[20 chad]);
    if isempty(tsq)
        k = k-1;
        break;
    end
    cont(k) = tsq(end);
    num_(k) = num;
end
a_ = a_(1:k);

idx = 0;
for k = 1:length(cont)-1
    if num_(k)*num_(k+1) < 0
        idx = k;
        break;
    end
end

if idx > 0
    additionalAcc = linspace(a_(idx),a_(idx+1),10);
    for k = 2:length(additionalAcc)-1
        [tsq,~,~,~,~,num,~] = calc_mintime_traj(init,final,[20 additionalAcc(k)]);
        if ~isempty(tsq)
            a_ = [a_ additionalAcc(k)];
            cont = [cont tsq(end)];
            num_ = [num_ num];
        else
            break;
        end
    end
end

[~,sortIdx] = sort(a_);
ACC = a_(sortIdx);
TF = cont(sortIdx);

figure(figNum);clf;
hold on
plot(ACC,TF,'o');
plot(ACC,TF,'-k');

%% checking singular range 
dTf_dAcc = diff(TF)./diff(ACC);
ddTf_dAcc = diff(dTf_dAcc);

checkSingular = [];
for k=1:length(ddTf_dAcc)-1
    checkSingular(k) = ddTf_dAcc(k)*ddTf_dAcc(k+1);
end

if ~isempty(checkSingular)
    [~,singularIdx] = min(checkSingular);
    isSingular = 0;
    if(checkSingular(singularIdx) < -10)
        isSingular = 1;
        k = singularIdx + 1;
        plot([ACC(k) ACC(k)],[0 20],':','linewidth',2);
        plot([ACC(k+1) ACC(k+1)],[0 20],':','linewidth',2);
        possibleRange{1} = [min(TF) TF(k+1)];
        possibleRange{2} = [TF(k) max(TF)];
    else
        possibleRange{1} = [min(TF) max(TF)];
    end
else
    possibleRange{1} = [min(TF) max(TF)];
end

end











