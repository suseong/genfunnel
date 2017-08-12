clear
clf

xinit1 = [0 0 0];  yinit1 = [-2 0 0]; zinit1 = [0 0 0];
xfinal1 = [5 0 -9]; yfinal1 = [0 1 0]; zfinal1 = [1 0 0];

xinit2 = [5 0 -9]; yinit2 = [0 1 0]; zinit2 = [1 0 0];
xfinal2 = [0 0 0]; yfinal2 = [2 0 0]; zfinal2 = [0 0 0];

[xtsq1,xpos1,xacc1,xiter1,xact1,~,~] = calc_mintime_traj(xinit1,xfinal1,[20 11.4]);
[ytsq1,ypos1,yacc1,yiter1,yact1,~,~] = calc_mintime_traj(yinit1,yfinal1,[20 1.8]);
[ztsq1,zpos1,zacc1,ziter1,zact1,~,~] = calc_mintime_traj(zinit1,zfinal1,[20 1.4]);

[xtsq2,xpos2,xacc2,xiter2,xact2,~,~] = calc_mintime_traj(xinit2,xfinal2,[20 11.4]);
[ytsq2,ypos2,yacc2,yiter2,yact2,~,~] = calc_mintime_traj(yinit2,yfinal2,[20 1.8]);
[ztsq2,zpos2,zacc2,ziter2,zact2,~,~] = calc_mintime_traj(zinit2,zfinal2,[20 1.4]);

[xp1,xv1] = show_pv(xtsq1,xacc1,xinit1,[20 11.4],1,100);
[yp1,yv1] = show_pv(ytsq1,yacc1,yinit1,[20 1.8],1,100);
[zp1,zv1] = show_pv(ztsq1,zacc1,zinit1,[20 1.4],1,100);

[xp2,xv2] = show_pv(xtsq2,xacc2,xinit2,[20 11.4],-1,100);
[yp2,yv2] = show_pv(ytsq2,yacc2,yinit2,[20 1.8],1,100);
[zp2,zv2] = show_pv(ztsq2,zacc2,zinit2,[20 1.4],-1,100);

%%
figure(13);clf;
hold on
plot3(xp1,yp1,zp1,'-k','linewidth',2)
plot3(xp2,yp2,zp2,'-k','linewidth',2)
xlabel('x')
ylabel('y')
zlabel('z')
axis equal
box on