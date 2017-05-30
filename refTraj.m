function [state_euler_traj,input_traj] = refTraj(xs,xf,tf,delT,plant)

ps = xs(1:3);
vs = xs(4:6);
% as = xs(7:9);
as = zeros(3,1);

pf = xf(1:3);
vf = xf(4:6);
% af = xf(7:9);
af = zeros(3,1);

delp = [];
delv = [];
dela = [];

calcCoeff = 1/tf^5*[      720   -360*tf   60*tf^2;
                      -360*tf  168*tf^2  -24*tf^3;
                      60*tf^2  -24*tf^3    3*tf^4];

traj_coeffs = [];

for i=1:3
    delp(i,1) = pf(i) - ps(i) - vs(i)*tf - 0.5*as(i)*tf^2;
    delv(i,1) = vf(i) - vs(i) - as(i)*tf;
    dela(i,1) = af(i) - as(i);
    temp = calcCoeff*[delp(i);delv(i);dela(i)];
    alpha = temp(1);
    beta  = temp(2);
    gamma = temp(3);
    traj_coefs(i,1) = alpha/120;
    traj_coefs(i,2) = beta/24;
    traj_coefs(i,3) = gamma/6;
    traj_coefs(i,4) = as(i)/2;
    traj_coefs(i,5) = vs(i);
    traj_coefs(i,6) = ps(i);
    traj_coefs(i+3,1) = 0;
    traj_coefs(i+3,2) = alpha/24;
    traj_coefs(i+3,3) = beta/6;
    traj_coefs(i+3,4) = gamma/2;
    traj_coefs(i+3,5) = as(i);
    traj_coefs(i+3,6) = vs(i);
end

stateTraj = PPTrajectory(mkpp([0;tf],traj_coefs,6));

tt = 0:delT:tf;

for i=1:length(tt)
    temp = stateTraj.dderiv(tt(i));
    T_ =  temp(1:3) - [0;0;plant.g];
    T_mag = -1.0*norm(T_);
    zb = T_ / T_mag;

    yaw_des = 0;

    xc = [cos(yaw_des);sin(yaw_des);0];
    yb = skew2Mat(zb)*xc/norm(skew2Mat(zb)*xc);
    xb = skew2Mat(yb)*zb;
    Rb = [xb yb zb];

%     [psi,theta,phi] = dcm2angle(Rb');
    agls = rotmat2rpy(Rb);
    phi = agls(1);theta = agls(2);psi = agls(3);

    a = stateTraj.dderiv(tt(i));
    hw = 1/T_mag*(a(4:6) - zb'*a(4:6)*zb);

    p = -hw'*yb; q = hw'*xb; r = 0;

    Euler(:,i) = [phi;theta;psi];
    input(:,i) = [T_mag*plant.m;p;q;r];
end

att_coef = [];
for i=1:length(tt)-1
    phi = Euler(1,i);theta = Euler(2,i);psi = Euler(3,i);
    invQ = [ 1   (sin(theta)*sin(phi))/cos(theta)  (cos(phi)*sin(theta))/cos(theta);
             0                           cos(phi)                         -sin(phi);
             0                sin(phi)/cos(theta)              cos(phi)/cos(theta)];    
    omg = input(2:4,i);
    dEuler = invQ*omg;
    phi_ = Euler(1,i+1);theta_ = Euler(2,i+1);psi_ = Euler(3,i+1);
    invQ_ = [ 1   (sin(theta_)*sin(phi_))/cos(theta_)  (cos(phi_)*sin(theta_))/cos(theta_);
             0                           cos(phi_)                         -sin(phi_);
             0                sin(phi_)/cos(theta_)              cos(phi_)/cos(theta_)];    
    omg_ = input(2:4,i+1);
    dEuler_ = invQ_*omg_;
    
    for j=1:3
        d = Euler(j,i); c = dEuler(j);
        dt = tt(i+1) - tt(i);
        temp = 1/dt^3*[-2 dt;3*dt -dt^2];
        ab = temp*[Euler(j,i+1)-c*dt-d;omg_(j)-c];
        att_coef = [att_coef;[0 0 ab(1) ab(2) c d]];
    end
end

input_pp = foh(tt,input);

state_euler_coefs = [];
for i=1:length(tt)-1
    t1 = tt(i);
    for j=1:6

        traj_coefs_ = zeros(1,6);
        traj_coefs_(1,6) = traj_coefs(j,:)*[   t1^5   t1^4   t1^3   t1^2 t1^1 t1^0]';
        traj_coefs_(1,5) = traj_coefs(j,:)*[ 5*t1^4 4*t1^3 3*t1^2 2*t1^1    1    0]';
        traj_coefs_(1,4) = traj_coefs(j,:)*[10*t1^3 6*t1^2 3*t1^1      1    0    0]';
        traj_coefs_(1,3) = traj_coefs(j,:)*[10*t1^2 4*t1^1      1      0    0    0]';
        traj_coefs_(1,2) = traj_coefs(j,:)*[ 5*t1^1 1*t1^0      0      0    0    0]';
        traj_coefs_(1,1) = traj_coefs(j,:)*[ 1*t1^0      0      0      0    0    0]';
        state_euler_coefs = [state_euler_coefs;traj_coefs_];    
    end

    state_euler_coefs = [state_euler_coefs;att_coef(3*i-2:3*i,:)];
end

state_euler_traj = PPTrajectory(mkpp(tt,state_euler_coefs,9));
state_euler_traj = state_euler_traj.setOutputFrame(plant.getStateFrame);
input_traj = PPTrajectory(input_pp);
input_traj = input_traj.setOutputFrame(plant.getInputFrame);

end
