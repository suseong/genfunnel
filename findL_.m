function [coeff1, coeff2] = findL_(x,V,Vdot,rho,rhodot,initRegion)

checkDependency('yalmip');
N = length(V)-1;
coeff1 = [];
coeff2 = [];

parfor i = 1:N
  if i == 1
      xs = sdpvar(9,1);
      [L1,coeff1_] = polynomial(xs,2);
      [L2,coeff2_] = polynomial(xs,2);

      Vdots = msspoly2sdpvar(x,xs,(Vdot{i}));
      Vs = msspoly2sdpvar(x,xs,(V{i}));

      c1 = sos(rhodot(i) - Vdots - L1*(rho(i) - Vs));
      c2 = sos(L2);
      c3 = sos(rho(i) - Vs - L2*(1-xs'*initRegion*xs));
      constraints = [c1 c2 c3];
      vars = [coeff1_; coeff2_];
      
      sol = solvesos(constraints,[],[],vars);
      if sol.problem == 0
         coeff1(:,i) = double(coeff1_);
         coeff2(:,i) = double(coeff2_);
      end
  else
      xs = sdpvar(9,1);
      [L1,coeff1_] = polynomial(xs,2);

      Vdots = msspoly2sdpvar(x,xs,(Vdot{i}));
      Vs = msspoly2sdpvar(x,xs,(V{i}));

      c1 = sos(rhodot(i) - Vdots - L1*(rho(i) - Vs));
      constraints = c1;
      vars = coeff1_;

      sol = solvesos(constraints,[],[],vars);
      if sol.problem == 0
         coeff1(:,i) = double(coeff1_);
      end
  end
end

end