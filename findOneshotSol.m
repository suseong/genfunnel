function oneshotAcc = findOneshotAcc(init,final,input)

pos = final(1) - init(1);
pSign = sign(pos);
xf = pSign*pos;

v0 = pSign*init(2); a0 = pSign*init(3);
vf = pSign*final(2); af = pSign*final(3);

%% case acc
% 1-1
dt_ = (vf - (a0^2 + af^2 - 2*am^2 + 2*u*v0)/(2*u))/am;
err = -(12*u^2*v0^2 - 12*u^2*vf^2 - 8*a0^3*am + 8*af^3*am + 3*a0^4 - 3*af^4 + 6*a0^2*am^2 - 6*af^2*am^2 + 12*a0^2*u*v0 + 12*am^2*u*v0 + 12*af^2*u*vf + 12*am^2*u*vf - 24*am*u^2*x0 + 24*am*u^2*xf - 24*a0*am*u*v0 - 24*af*am*u*vf)/(24*am*u^2)

% 1-2
dt_ = (vf - (a0^2 - af^2 + 2*u*v0)/(2*u))/am;
err = -(12*u^2*v0^2 - 12*u^2*vf^2 - 8*a0^3*am + 8*af^3*am + 3*a0^4 - 3*af^4 + 6*a0^2*am^2 - 6*af^2*am^2 + 12*a0^2*u*v0 + 12*am^2*u*v0 - 12*af^2*u*vf - 12*am^2*u*vf - 24*am*u^2*x0 + 24*am*u^2*xf - 24*a0*am*u*v0 + 24*af*am*u*vf)/(24*am*u^2)

% 2-1
dt_ = (vf - (- a0^2 + af^2 + 2*u*v0)/(2*u))/am;
err = -(12*u^2*v0^2 - 12*u^2*vf^2 - 8*a0^3*am + 8*af^3*am + 3*a0^4 - 3*af^4 + 6*a0^2*am^2 - 6*af^2*am^2 - 12*a0^2*u*v0 - 12*am^2*u*v0 + 12*af^2*u*vf + 12*am^2*u*vf - 24*am*u^2*x0 + 24*am*u^2*xf + 24*a0*am*u*v0 - 24*af*am*u*vf)/(24*am*u^2)

% 2-2
dt_ = (vf - (- a0^2 - af^2 + 2*am^2 + 2*u*v0)/(2*u))/am;
err = (12*u^2*vf^2 - 12*u^2*v0^2 + 8*a0^3*am - 8*af^3*am - 3*a0^4 + 3*af^4 - 6*a0^2*am^2 + 6*af^2*am^2 + 12*a0^2*u*v0 + 12*am^2*u*v0 + 12*af^2*u*vf + 12*am^2*u*vf + 24*am*u^2*x0 - 24*am*u^2*xf - 24*a0*am*u*v0 - 24*af*am*u*vf)/(24*am*u^2)

%% case dec
% 1-1
dt_ = -(vf - (a0^2 + af^2 - 2*am^2 + 2*u*v0)/(2*u))/am;
err = (12*u^2*v0^2 - 12*u^2*vf^2 + 8*a0^3*am - 8*af^3*am + 3*a0^4 - 3*af^4 + 6*a0^2*am^2 - 6*af^2*am^2 + 12*a0^2*u*v0 + 12*am^2*u*v0 + 12*af^2*u*vf + 12*am^2*u*vf + 24*am*u^2*x0 - 24*am*u^2*xf + 24*a0*am*u*v0 + 24*af*am*u*vf)/(24*am*u^2)

% 1-2
dt_ = -(vf - (a0^2 - af^2 + 2*u*v0)/(2*u))/am;
err = (12*u^2*v0^2 - 12*u^2*vf^2 + 8*a0^3*am - 8*af^3*am + 3*a0^4 - 3*af^4 + 6*a0^2*am^2 - 6*af^2*am^2 + 12*a0^2*u*v0 + 12*am^2*u*v0 - 12*af^2*u*vf - 12*am^2*u*vf + 24*am*u^2*x0 - 24*am*u^2*xf + 24*a0*am*u*v0 - 24*af*am*u*vf)/(24*am*u^2)

% 2-1
dt_ = -(vf - (- a0^2 + af^2 + 2*u*v0)/(2*u))/am;
err = (12*u^2*v0^2 - 12*u^2*vf^2 + 8*a0^3*am - 8*af^3*am + 3*a0^4 - 3*af^4 + 6*a0^2*am^2 - 6*af^2*am^2 - 12*a0^2*u*v0 - 12*am^2*u*v0 + 12*af^2*u*vf + 12*am^2*u*vf + 24*am*u^2*x0 - 24*am*u^2*xf - 24*a0*am*u*v0 + 24*af*am*u*vf)/(24*am*u^2)

% 2-2
dt_ = -(vf - (- a0^2 - af^2 + 2*am^2 + 2*u*v0)/(2*u))/am;
err = -(12*u^2*vf^2 - 12*u^2*v0^2 - 8*a0^3*am + 8*af^3*am - 3*a0^4 + 3*af^4 - 6*a0^2*am^2 + 6*af^2*am^2 + 12*a0^2*u*v0 + 12*am^2*u*v0 + 12*af^2*u*vf + 12*am^2*u*vf - 24*am*u^2*x0 + 24*am*u^2*xf + 24*a0*am*u*v0 + 24*af*am*u*vf)/(24*am*u^2)


end