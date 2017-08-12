function Acc = get_acc(accList,tfList,figNum,Tf,drawFigure)

for k=1:length(tfList)-1
    if and(Tf <= tfList(k), Tf >= tfList(k+1))
        Acc = accList(k) + (tfList(k) - Tf)/(tfList(k) - tfList(k+1)) * (accList(k+1) - accList(k));
        break;
    end
end

if drawFigure
    figure(figNum)
    hold on
    plot(Acc,Tf,'x','markersize',15,'linewidth',5);
end

end