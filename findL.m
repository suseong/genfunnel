function [coeff1_,coeff3_,S_] = findL(P,Q,rho,rhodot,initRegion,Kp,Kd)

checkDependency('yalmip');
monomialOrder = 2;
N = length(rho)-1;
coeff1_ = [];
coeff3_ = [];
S_ = [];
Er = 0.08;

maxPpv = max(max(P(1:3,4:6)));
% maxPp = max(max(P(1:3,1:3)));
maxPv = max(max(P(4:6,4:6)));
maxKp = max(max(Kp));
maxKd = max(max(Kd));

parfor i=1:N
   if i==1
       e = sdpvar(6,1);
       S = sdpvar(6,6);
       epbar = sdpvar(1,1);
       edbar = sdpvar(1,1);
       epep = e(1:3)'*e(1:3);
       eded = e(4:6)'*e(4:6); 
       
       [L1,coeff1] = polynomial([e;epbar;edbar],monomialOrder);
       [L2,coeff2] = polynomial(e,monomialOrder);
       [L3,coeff3] = polynomial(e,monomialOrder);

       [L4,coeff4] = polynomial([e;epbar;edbar],monomialOrder);
       [L5,coeff5] = polynomial([e;epbar;edbar],monomialOrder);
       [L6,coeff6] = polynomial([e;epbar;edbar],monomialOrder);
       [L7,coeff7] = polynomial([e;epbar;edbar],monomialOrder);
       
       V = e'*P*e;
       Vdot = e'*Q*e ...
              + Er*(maxPpv*epbar+maxPv*edbar)*(maxKp*epbar+maxKd*edbar + 18);
   
       c1 = sos(rhodot(i) - Vdot ...
                - L1*(rho(i) - V) ...
                - L4*(epbar^2 - epep) ...
                - L5*(edbar^2 - eded) ...
                - L6*(edbar) ...
                - L7*(epbar));    
       c2 = sos(L2);
       c3 = sos(rho(i) - V - L2*(1-e'*initRegion*e));
       c4 = sos(L3);
       c5 = sos(1 - e'*S*e - L3*(rho(i) - V)); 
       c6 = sos(L6);
       c7 = sos(L7);
       c8 = S >= 0;
       
       vars = [coeff1;coeff2;coeff3;coeff4;coeff5;coeff6;coeff7;S(:)];

       sol = solvesos([c1 c2 c3 c4 c5 c6 c7 c8],-geomean(S),[],vars);

       if sol.problem == 0
           coeff1_(:,i) = value(coeff1);
           coeff3_(:,i) = value(coeff3);
           S_(:,i) = value(S(:));
       else
           coeff1_(:,i) = zeros(size(coeff1));
           coeff3_(:,i) = zeros(size(coeff3));           
           S_(:,i) = zeros(size(S(:)));           
       end
   else       
       e = sdpvar(6,1);
       S = sdpvar(6,6);
       epbar = sdpvar(1,1);
       edbar = sdpvar(1,1);
       epep = e(1:3)'*e(1:3);
       eded = e(4:6)'*e(4:6);
       
%        [L1,coeff1] = polynomial(e,monomialOrder);
%        [L3,coeff3] = polynomial(e,monomialOrder);

       [L1,coeff1] = polynomial([e;epbar;edbar],monomialOrder);
       [L3,coeff3] = polynomial(e,monomialOrder);

       [L4,coeff4] = polynomial([e;epbar;edbar],monomialOrder);
       [L5,coeff5] = polynomial([e;epbar;edbar],monomialOrder);
       [L6,coeff6] = polynomial([e;epbar;edbar],monomialOrder);
       [L7,coeff7] = polynomial([e;epbar;edbar],monomialOrder);
       
       V = e'*P*e;
       Vdot = e'*Q*e ...
              + Er*(maxPpv*epbar+maxPv*edbar)*(maxKp*epbar+maxKd*edbar + 18);
   
       c1 = sos(rhodot(i) - Vdot ...
                - L1*(rho(i) - V) ...
                - L4*(epbar^2 - epep) ...
                - L5*(edbar^2 - eded) ...
                - L6*(edbar) ...
                - L7*(epbar));    
       c4 = sos(L3);
       c5 = sos(1 - e'*S*e - L3*(rho(i) - V)); 
       c6 = sos(L6);
       c7 = sos(L7);
       c8 = S >= 0;
       
       vars = [coeff1;coeff3;coeff4;coeff5;coeff6;coeff7;S(:)];
       sol = solvesos([c1 c4 c5 c6 c7 c8],-geomean(S),[],vars);
         
       if sol.problem == 0
           coeff1_(:,i) = value(coeff1);
           coeff3_(:,i) = value(coeff3);
           S_(:,i) = value(S(:));
       else
           coeff1_(:,i) = zeros(size(coeff1));
           coeff3_(:,i) = zeros(size(coeff3)); 
           S_(:,i) = zeros(size(S(:)));           
       end
   end 
end

end