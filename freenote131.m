clear all
close all
% clc

xinit =  [ 5,-2, 0];
xfinal = [ 2, 1, 0];
% xinit =  [ 0, 3, 0];
% xfinal = [ 2, 2, 0];

yinit =  [ 0, 3, 0];
yfinal = [ 2, 2, 0];

zinit =  [-2, 3, 0];
zfinal = [ 0, -2, 0];
% zinit =  [ 0, 3, 0];
% zfinal = [ 2, 2, 0];

amax = 20;
g = 9.8;

az = 0.2;
ax = 0.5;

%%
for k=1:10000

azmax = az*(amax - g);
axmax = ax*sqrt(amax^2 - (azmax + g)^2);
aymax = sqrt(amax^2 - axmax^2 - (azmax + g)^2);

[xtsq,xpos,xacc,xiter,xact] = calc_mintime_traj(xinit,xfinal,[40 axmax]);
[ytsq,ypos,yacc,yiter,yact] = calc_mintime_traj(yinit,yfinal,[40 aymax]);
[ztsq,zpos,zacc,ziter,zact] = calc_mintime_traj(zinit,zfinal,[40 azmax]);
    
xtf = xtsq(end);
ytf = ytsq(end);
ztf = ztsq(end);

kk = 1;

kaz = 0.001;
kax = 0.001;

if and( xtf < ztf, ytf < ztf)
    az = az + kaz;
elseif and( xtf > ztf, ytf > ztf)
    az = az - kaz;
else
    if (xtf < ytf)
        ax = ax - kax;
    else
        ax = ax + kax;
    end
end

ax = min(max(ax,0.0001),0.9999);
az = min(max(az,0.0001),0.9999);

disp([num2str(xtf),'   ',num2str(ytf),'   ',num2str(ztf),'   ',num2str(axmax),'   ',num2str(aymax),'   ',num2str(azmax),'   ',num2str(ax),'   ',num2str(az)]);

if and(abs(xtf - ytf) < 1e-3, abs(ytf - ztf) < 1e-3)
    break;
end

end

%%
[xp,xv] = show_pv(xtsq,xacc,xinit,[40 axmax],xact,100);
[yp,yv] = show_pv(ytsq,yacc,yinit,[40 aymax],yact,100);
[zp,zv] = show_pv(ztsq,zacc,zinit,[40 azmax],yact,100);

figure(13)
hold on
plot3(xp(1),yp(1),zp(1),'*')
plot3(xp,yp,zp)

xlabel('x')
ylabel('y')
zlabel('z')

axis equal;


