function [state_euler_traj,input_traj] = refTraj2(xs,xf,tf,delT)

ps = xs(1:3);
vs = xs(4:6);
as = zeros(3,1);
js = zeros(3,1);
% as = zeros(3,1);

pf = xf(1:3);
vf = xf(4:6);
af = zeros(3,1);
jf = zeros(3,1);
% af = zeros(3,1);

delp = [];
delv = [];
dela = [];
delj = [];

calcCoeff = 1/tf^7*[-20*tf^0  10*tf^1   -2*tf^2  1/6*tf^3;
                     70*tf^1 -34*tf^2 13/2*tf^3 -1/2*tf^4;
                    -84*tf^2  39*tf^3   -7*tf^4  1/2*tf^5;
                     35*tf^3 -15*tf^4  5/2*tf^5 -1/6*tf^6];

traj_coeffs = [];

for i=1:3
    delp(i,1) = pf(i) - ps(i) - vs(i)*tf - 0.5*as(i)*tf^2 - js(i)*tf^3/6;
    delv(i,1) = vf(i) - vs(i) - as(i)*tf - 0.5*js(i)*tf^2;
    dela(i,1) = af(i) - as(i) - js(i)*tf;
    delj(i,1) = jf(i) - js(i);
    temp = calcCoeff*[delp(i);delv(i);dela(i);delj(i)];
    alpha = temp(1);
    beta  = temp(2);
    gamma = temp(3);
    lambda = temp(4);
    traj_coefs(i,1) = alpha;
    traj_coefs(i,2) = beta;
    traj_coefs(i,3) = gamma;
    traj_coefs(i,4) = lambda;
    traj_coefs(i,5) = js(i)/6;
    traj_coefs(i,6) = as(i)/2;
    traj_coefs(i,7) = vs(i);
    traj_coefs(i,8) = ps(i);    
    traj_coefs(i+3,1) = 0;
    traj_coefs(i+3,2) = 7*alpha;
    traj_coefs(i+3,3) = 6*beta;
    traj_coefs(i+3,4) = 5*gamma;
    traj_coefs(i+3,5) = 4*lambda;
    traj_coefs(i+3,6) = js(i)/2;
    traj_coefs(i+3,7) = as(i);
    traj_coefs(i+3,8) = vs(i);    
end

stateTraj = PPTrajectory(mkpp([0;tf],traj_coefs,6));

tt = 0:delT:tf;

for i=1:length(tt)
    temp = stateTraj.dderiv(tt(i));
    T_ =  temp(1:3) - [0;0;-9.8];
    T_mag = 1.0*norm(T_);
    zb = T_ / T_mag;

    v = stateTraj.deriv(tt(i));
    yaw_des = 0;

    xc = [cos(yaw_des);sin(yaw_des);0];
    yb = skew2Mat(zb)*xc/norm(skew2Mat(zb)*xc);
    xb = skew2Mat(yb)*zb;
    Rb = [xb yb zb];

    rpy = rotmat2rpy(Rb);
    phi = rpy(1);theta = rpy(2);psi = rpy(3);
    Euler(:,i) = [phi;theta;psi];
    p = 0;q = 0;r = 0;
    input(:,i) = [T_mag;p;q;r];
end

for i=1:length(tt)-1
   R_ = rot(Euler(1,i+1),Euler(2,i+1),Euler(3,i+1));
   R =  rot(Euler(1,i),Euler(2,i),Euler(3,i));
   
   dR = (R_-R)/(tt(i+1)-tt(i));
   omg = R'*dR;
   p = omg(3,2);
   q = omg(1,3);
   r = omg(2,1);
   input(2,i) = p;input(3,i) = q;input(4,i) = r;
end
input(2,i+1) = p;input(3,i+1) = q;input(4,i+1) = r;

Euler_pp = foh(tt,Euler);
input_pp = foh(tt,input);

state_euler_coefs = [];
for i=1:length(tt)-1
    t1 = tt(i);
    for j=1:6
%         traj_coefs_ = traj_coefs(j,:);
        traj_coefs_ = zeros(1,8);
        traj_coefs_(1,8) = traj_coefs(j,:)*[    t1^7      t1^6      t1^5     t1^4     t1^3    t1^2   t1^1  t1^0]';
        traj_coefs_(1,7) = traj_coefs(j,:)*[  7*t1^6    6*t1^5    5*t1^4   4*t1^3   3*t1^2  2*t1^1   t1^0     0]';
        traj_coefs_(1,6) = traj_coefs(j,:)*[ 21*t1^5   15*t1^4   10*t1^3   6*t1^2   3*t1^1    t1^0      0     0]';
        traj_coefs_(1,5) = traj_coefs(j,:)*[ 35*t1^4   20*t1^3   10*t1^2   4*t1^1     t1^0       0      0     0]';
        traj_coefs_(1,4) = traj_coefs(j,:)*[ 35*t1^3   15*t1^2    5*t1^1     t1^0        0       0      0     0]';
        traj_coefs_(1,3) = traj_coefs(j,:)*[ 21*t1^2    6*t1^1      t1^0        0        0       0      0     0]';
        traj_coefs_(1,2) = traj_coefs(j,:)*[  7*t1^1      t1^0         0        0        0       0      0     0]';
        traj_coefs_(1,1) = traj_coefs(j,:)*[  1*t1^0         0         0        0        0       0      0     0]';

        state_euler_coefs = [state_euler_coefs;traj_coefs_];    
    end
    state_euler_coefs = [state_euler_coefs;[zeros(3,6) Euler_pp.coefs(3*i-2:3*i,:)]];
end

state_euler_traj = PPTrajectory(mkpp(tt,state_euler_coefs,9));
% state_euler_traj = state_euler_traj.setOutputFrame(plant.getStateFrame);
input_traj = PPTrajectory(input_pp);
% input_traj = input_traj.setOutputFrame(plant.getInputFrame);

end
