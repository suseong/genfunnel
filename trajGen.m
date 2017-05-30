function [xtraj,utraj] = trajGen(p,xs,xe,duration)

% p = quadRPG();

x0 = Point(getStateFrame(p));
x0.x = xs(1);
x0.y = xs(2);
x0.z = xs(3);
x0.dx = xs(4);
x0.dy = xs(5);
x0.dz = xs(6);
x0.phi = xs(7);
x0.theta = xs(8);
x0.psi = xs(9);

xf = x0;
xf.x = xe(1);
xf.y = xe(2);
xf.z = xe(3);
xf.dx = xe(4);
xf.dy = xe(5);
xf.dz = xe(6);
xf.phi = xe(7);
xf.theta = xe(8);
xf.psi = xe(9);

u0 = double([-1.0*p.m*p.g;0;0;0]);

N = 20;
minimum_duration = .1;
maximum_duration = duration;

prog = DircolTrajectoryOptimization(p,N,[minimum_duration maximum_duration]);  
prog = prog.addStateConstraint(ConstantConstraint(double(x0)),1);
prog = prog.addInputConstraint(ConstantConstraint(u0),1);
prog = prog.addStateConstraint(ConstantConstraint(double(xf)),N);
prog = prog.addInputConstraint(ConstantConstraint(u0),N);
prog = prog.addRunningCost(@cost);
prog = prog.addFinalCost(@finalCost);

tf0 = duration;                    % initial guess at duration 
traj_init.x = PPTrajectory(foh([0,tf0],[double(x0),double(xf)]));
traj_init.u = ConstantTrajectory(u0);

info=0;
while (info~=1)
  tic
  [xtraj,utraj,~,~,info] = prog.solveTraj(tf0,traj_init);
  toc
end

end