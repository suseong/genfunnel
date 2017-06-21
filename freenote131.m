clear all
close all
clc

xinit = [0,0,0];
xfinal = [3,1,0];
yinit = [0,0,0];
yfinal = [1,4,0];
zinit = [0,0,0];
zfinal = [5,1,0];

amax = 20;
g = 9.8;

az = 0.5;
ax = 0.5;

%%
for k=1:100

azmax = az*(amax - g);
axmax = ax*sqrt(amax^2 - (azmax + g)^2);
aymax = sqrt(amax^2 - axmax^2 - (azmax + g)^2);

[xtsq,xpos,xacc,xiter,xact] = calc_mintime_traj(xinit,xfinal,[8 axmax]);
[ytsq,ypos,yacc,yiter,yact] = calc_mintime_traj(yinit,yfinal,[8 aymax]);
[ztsq,zpos,zacc,ziter,zact] = calc_mintime_traj(zinit,zfinal,[8 azmax]);
    
xtf = xtsq(end);
ytf = ytsq(end);
ztf = ztsq(end);

kk = 1;

if and( xtf < ztf, ytf < ztf)

    if or(abs(ztf - ytf) < 0.03, abs(ztf - xtf) < 0.03)
        kaz = 0.1;
    else
        kaz = 0.1;
    end    
    az = az + kaz*abs(ztf - ytf)^kk;
    
elseif and( xtf > ztf, ytf > ztf)

    if or(abs(ztf - ytf) < 0.03, abs(ztf - xtf) < 0.03)
        kaz = 0.1;
    else
        kaz = 0.1;
    end
    az = az - kaz*abs(ztf - ytf)^kk;
    
else

    if abs(xtf - ytf) < 0.03
        kax = 0.1;
    else
        kax = 0.1;
    end
    
    if (xtf < ytf)
        ax = ax - kax*abs(xtf-ytf)^kk;
    else
        ax = ax + kax*abs(xtf-ytf)^kk;
    end
end 

disp([num2str(xtf),' ',num2str(ytf),' ',num2str(ztf)]);

if and(abs(xtf - ytf) < 1e-3, abs(ytf - ztf) < 1e-3)
    break;
end

end

%%
[xp,xv] = show_pv(xtsq,xacc,xinit,[8 axmax],xact,100);
[yp,yv] = show_pv(ytsq,yacc,yinit,[8 aymax],yact,100);
[zp,zv] = show_pv(ztsq,zacc,zinit,[8 azmax],yact,100);

figure(13)
plot3(xp,yp,zp)
axis equal;






