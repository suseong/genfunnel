clear
close
clc

xinit = [0,2,0];
xfinal = [2,0,0];
yinit = [0,0,0];
yfinal = [2,2,0];

actLimit = [10,5];

[xtsq,xpos,xacc,xiter,xact] = calc_mintime_traj(xinit,xfinal,actLimit);
[ytsq,ypos,yacc,yiter,yact] = calc_mintime_traj(yinit,yfinal,actLimit);

%%
[xp,xv] = show_pv(xtsq,xacc,xinit,actLimit,xact,100);
[yp,yv] = show_pv(ytsq,yacc,yinit,actLimit,yact,100);

figure(1)
subplot(3,1,1)
plot(xp,yp)
axis equal

subplot(3,1,2)
tt = linspace(0,xtsq(end),100);
plot(tt,xv)

subplot(3,1,3)
tt = linspace(0,ytsq(end),100);
plot(tt,yv)

%%
% pts = trackInfo;
% 
% figure(1)
% hold on
% for k=1:length(pts)-1
%    plot3([pts{k}(1) pts{k+1}(1)],[pts{k}(2) pts{k+1}(2)],[pts{k}(3) pts{k+1}(3)],'-')
%    plot3(pts{k}(1),pts{k}(2),pts{k}(3),'*')
% end
% axis equal
