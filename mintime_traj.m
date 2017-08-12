function [tsq,pos,acc,act] = mintime_traj(init,final)

maxAcc = 20;
    
xinit  =  init{1}; yinit  =  init{2}; zinit  =  init{3};
xfinal = final{1}; yfinal = final{2}; zfinal = final{3};

%%
[xacc,xtf,xrange] = calc_acc_tf(xinit,xfinal,100);
[yacc,ytf,yrange] = calc_acc_tf(yinit,yfinal,101);
[zacc,ztf,zrange] = calc_acc_tf(zinit,zfinal,102);

tfRange = [];
minTf = 100;

cnt = 1;
for j = 1:length(xrange)
    for k = 1:length(yrange)
        for l = 1:length(zrange)
            range(1) = max([xrange{j}(1) yrange{k}(1) zrange{l}(1)]);
            range(2) = min([xrange{j}(2) yrange{k}(2) zrange{l}(2)]);
            if(range(1) < range(2))
                xacc_ = get_acc(xacc,xtf,100,range(2),0);
                yacc_ = get_acc(yacc,ytf,101,range(2),0);
                zacc_ = get_acc(zacc,ztf,102,range(2),0);
                if sqrt(xacc_^2 + yacc_^2 + (zacc_+9.8)^2) <= maxAcc                    
                    tfRange{cnt} = [range(1) range(2)];
                    minTf = min([minTf,range(1)]);
                    cnt = cnt+1;
                end
            end
        end
    end
end
    
validTfRegion = [100 200];

if ~isempty(tfRange)
    for j = 1:length(tfRange)
        if validTfRegion(1) > tfRange{j}(1)
            validTfRegion = tfRange{j};
        end
    end
    
    xMaxAcc = get_acc(xacc,xtf,100,validTfRegion(1),1);
    yMaxAcc = get_acc(yacc,ytf,101,validTfRegion(1),1);
    zMaxAcc = get_acc(zacc,ztf,102,validTfRegion(1),1);
    
    xMinAcc = get_acc(xacc,xtf,100,validTfRegion(2),1);
    yMinAcc = get_acc(yacc,ytf,101,validTfRegion(2),1);
    zMinAcc = get_acc(zacc,ztf,102,validTfRegion(2),1);
    
    %%
    xtsq = calc_mintime_traj(xinit,xfinal,[20 xMaxAcc]);
    ytsq = calc_mintime_traj(yinit,yfinal,[20 yMaxAcc]);
    ztsq = calc_mintime_traj(zinit,zfinal,[20 zMaxAcc]);
    
    disp(['final times are : ',num2str(xtsq(end)),', ',num2str(ytsq(end)),', ',num2str(ztsq(end))]);
    % disp(['max final time is : ',num2str(max([xtsq(end) ytsq(end) ztsq(end)]))]);
    
    if xtsq(end) ~= max([xtsq(end) ytsq(end) ztsq(end)])
        xAccRefined = find_valid_acc(xinit,xfinal,[20 0],[xMinAcc xMaxAcc],max([xtsq(end) ytsq(end) ztsq(end)]));
    else
        xAccRefined = xMaxAcc;
    end
    
    if ytsq(end) ~= max([xtsq(end) ytsq(end) ztsq(end)])
        yAccRefined = find_valid_acc(yinit,yfinal,[20 0],[yMinAcc yMaxAcc],max([xtsq(end) ytsq(end) ztsq(end)]));
    else
        yAccRefined = yMaxAcc;
    end
    
    if ztsq(end) ~= max([xtsq(end) ytsq(end) ztsq(end)])
        zAccRefined = find_valid_acc(zinit,zfinal,[20 0],[zMinAcc zMaxAcc],max([xtsq(end) ytsq(end) ztsq(end)]));
    else
        zAccRefined = zMaxAcc;
    end
    
    [xtsq,xpos,xacc,~,~,xnum] = calc_mintime_traj(xinit,xfinal,[20 xAccRefined]);
    [ytsq,ypos,yacc,~,~,ynum] = calc_mintime_traj(yinit,yfinal,[20 yAccRefined]);
    [ztsq,zpos,zacc,~,~,znum] = calc_mintime_traj(zinit,zfinal,[20 zAccRefined]);
    
    disp(['refined final times are : ',num2str(xtsq(end)),', ',num2str(ytsq(end)),', ',num2str(ztsq(end))]);
    
    [xp1,~] = show_pv(xtsq,xacc,xinit,[20 xAccRefined],-xnum,100);
    [yp1,~] = show_pv(ytsq,yacc,yinit,[20 yAccRefined],-ynum,100);
    [zp1,~] = show_pv(ztsq,zacc,zinit,[20 zAccRefined],-znum,100);
    
    tsq{1} =  xtsq; tsq{2} =  ytsq; tsq{3} =  ztsq;
    pos{1} =  xpos; pos{2} =  ypos; pos{2} =  zpos;
    acc{1} =  xacc; acc{2} =  yacc; acc{3} =  zacc;
    act{1} = -xnum; act{2} = -ynum; act{3} = -znum;
    
    figure(150);clf;
    subplot(3,2,[1 3 5])
    hold on
    plot3(xp1,yp1,zp1,'-k','linewidth',2)
    xlabel('x')
    ylabel('y')
    zlabel('z')
    axis equal
    box on
    subplot(3,2,2)
    plot(xp1)
    subplot(3,2,4)
    plot(yp1)
    subplot(3,2,6)
    plot(zp1)
    
    
else
    tsq = []; pos = []; acc = []; act = [];
    disp('no valid solution found');
end



