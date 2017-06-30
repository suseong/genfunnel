clear all
close all
% clc

xinit =  [ 5,-2, 0];
xfinal = [ 0, 0, 0];
% xinit =  [ 0, 3, 0];
% xfinal = [ 2, 2, 0];

yinit =  [ -10, -3, 0];
yfinal = [ 0, 0, 0];

zinit =  [-2, 3, 0];
zfinal = [ 0, 0, 0];
% zinit =  [ 0, 3, 0];
% zfinal = [ 2, 2, 0];

amax = 20;
g = 9.8;

az = 0.2;
ax = 0.5;

%%
for k=1:100
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
    
    kaz = 0.0001;
    kax = 0.05;

    if or(and( xtf > ztf, ytf < ztf), and( xtf < ztf, ytf > ztf))
        smaller = [];
        bigger = [];
        
        ax_ = 0.01:0.1:0.99; ax_ = [ax_ 0.999];
        
        for kk = 1:length(ax_)                
            azmax = az*(amax - g);
            axmax = ax_(kk)*sqrt(amax^2 - (azmax + g)^2);
            aymax = sqrt(amax^2 - axmax^2 - (azmax + g)^2);
            
            [xtsq,xpos,xacc,xiter,xact] = calc_mintime_traj(xinit,xfinal,[40 axmax]);
            [ytsq,ypos,yacc,yiter,yact] = calc_mintime_traj(yinit,yfinal,[40 aymax]);
            
            xtf = xtsq(end);
            ytf = ytsq(end);

            chad(kk) = xtf - ytf;
        end
        
        for kk = 1:length(ax_)-1
            if(chad(kk)*chad(kk+1) < 0)
                if chad(kk) < chad(kk+1)
                    smaller = ax_(kk); bigger = ax_(kk+1);
                else
                    smaller = ax_(kk+1); bigger = ax_(kk);
                end
                break;
            end
        end
        
        for kk = 1:30
            ax_ = (smaller + bigger)/2;
            
            azmax = az*(amax - g);
            axmax = ax_*sqrt(amax^2 - (azmax + g)^2);
            aymax = sqrt(amax^2 - axmax^2 - (azmax + g)^2);
            
            [xtsq,xpos,xacc,xiter,xact] = calc_mintime_traj(xinit,xfinal,[40 axmax]);
            [ytsq,ypos,yacc,yiter,yact] = calc_mintime_traj(yinit,yfinal,[40 aymax]);
            
            xtf = xtsq(end);
            ytf = ytsq(end);
            
            if xtf - ytf <= 0 
                smaller = ax_;
            else
                bigger = ax_;
            end
            
            if or(abs(xtf - ytf) < 1e-3, abs(bigger - smaller) < 1e-3)
                ax = ax_;
                break;
            end
        end        
%         cf = max(abs(ztf-xtf),abs(ztf-ytf));
%         az = az + cf*kaz;
    elseif and( xtf > ztf, ytf > ztf)
        
        smaller = [];
        bigger = [];
        
        az_ = 0.01:0.02:0.2; 
%         az_ = [az_ 0.999];
        
        for kk = 1:length(az_)                
            azmax = az_(kk)*(amax - g);
            axmax = ax*sqrt(amax^2 - (azmax + g)^2);
            aymax = sqrt(amax^2 - axmax^2 - (azmax + g)^2);
            
            if ytf >= xtf
                [tsq,~,~,~,~] = calc_mintime_traj(xinit,xfinal,[40 axmax]);
                [ztsq,zpos,zacc,ziter,zact] = calc_mintime_traj(zinit,zfinal,[40 azmax]);
            else
                [tsq,~,~,~,~] = calc_mintime_traj(yinit,yfinal,[40 aymax]);
                [ztsq,zpos,zacc,ziter,zact] = calc_mintime_traj(zinit,zfinal,[40 azmax]);                
            end
            
            tf = tsq(end);
            ztf = ztsq(end);

            chad(kk) = tf - ztf;
        end
        
        for kk = 1:length(az_)-1
            if(chad(kk)*chad(kk+1) < 0)
                if chad(kk) < chad(kk+1)
                    smaller = az_(kk); bigger = az_(kk+1);
                else
                    smaller = az_(kk+1); bigger = az_(kk);
                end
                break;
            end
        end
        
        for kk = 1:30
            az_ = (smaller + bigger)/2;
            
            azmax = az_*(amax - g);
            axmax = ax*sqrt(amax^2 - (azmax + g)^2);
            aymax = sqrt(amax^2 - axmax^2 - (azmax + g)^2);

            if ytf >= xtf
                [tsq,~,~,~,~] = calc_mintime_traj(xinit,xfinal,[40 axmax]);
                [ztsq,zpos,zacc,ziter,zact] = calc_mintime_traj(zinit,zfinal,[40 azmax]);
            else
                [tsq,~,~,~,~] = calc_mintime_traj(yinit,yfinal,[40 aymax]);
                [ztsq,zpos,zacc,ziter,zact] = calc_mintime_traj(zinit,zfinal,[40 azmax]);                
            end
            
            tf = tsq(end);
            ztf = ztsq(end);
            
            if tf - ztf <= 0 
                smaller = az_;
            else
                bigger = az_;
            end
            
            if or(abs(tf - ztf) < 1e-3, abs(smaller - bigger) < 1e-3)
                az = az_;
                break;
            end
        end        
                
%         cf = max(abs(ztf-xtf),abs(ztf-ytf));
%         az = az - cf*kaz;
    else
        smaller = [];
        bigger = [];
        
        az_ = 0.1:0.1:0.99; az_ = [az_ 0.999];
        
        for kk = 1:length(az_)                
            azmax = az_(kk)*(amax - g);
            axmax = ax*sqrt(amax^2 - (azmax + g)^2);
            aymax = sqrt(amax^2 - axmax^2 - (azmax + g)^2);
            
            if ytf >= xtf
                [tsq,~,~,~,~] = calc_mintime_traj(yinit,yfinal,[40 aymax]);
                [ztsq,zpos,zacc,ziter,zact] = calc_mintime_traj(zinit,zfinal,[40 azmax]);
            else
                [tsq,~,~,~,~] = calc_mintime_traj(xinit,xfinal,[40 axmax]);
                [ztsq,zpos,zacc,ziter,zact] = calc_mintime_traj(zinit,zfinal,[40 azmax]);                
            end
            
            tf = tsq(end);
            ztf = ztsq(end);

            chad(kk) = tf - ztf;
        end
        
        for kk = 1:length(az_)-1
            if(chad(kk)*chad(kk+1) < 0)
                if chad(kk) < chad(kk+1)
                    smaller = az_(kk); bigger = az_(kk+1);
                else
                    smaller = az_(kk+1); bigger = az_(kk);
                end
                break;
            end
        end
        
        for kk = 1:30
            az_ = (smaller + bigger)/2;
            
            azmax = az_*(amax - g);
            axmax = ax*sqrt(amax^2 - (azmax + g)^2);
            aymax = sqrt(amax^2 - axmax^2 - (azmax + g)^2);

            if ytf >= xtf
                [tsq,~,~,~,~] = calc_mintime_traj(yinit,yfinal,[40 aymax]);
                [ztsq,zpos,zacc,ziter,zact] = calc_mintime_traj(zinit,zfinal,[40 azmax]);
            else
                [tsq,~,~,~,~] = calc_mintime_traj(xinit,xfinal,[40 axmax]);
                [ztsq,zpos,zacc,ziter,zact] = calc_mintime_traj(zinit,zfinal,[40 azmax]);                
            end
            
            tf = tsq(end);
            ztf = ztsq(end);
            
            if tf - ztf <= 0 
                smaller = az_;
            else
                bigger = az_;
            end
            
            if or(abs(tf - ztf) < 1e-3, abs(smaller - bigger) < 1e-3)
                az = az_;
                break;
            end
        end                
        
        

    end
    
    ax = min(max(ax,0.01),0.99);
    az = min(max(az,0.01),0.99);
    
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


