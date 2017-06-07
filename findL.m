function [coeff1_,coeff3_,S_] = findL(dt,P,Q,rho,rhodot,Kp,Kd)

checkDependency('yalmip');
monomialOrder = 2;
N = length(rho)-1;
coeff1_ = [];
coeff3_ = [];
S_ = [];

parfor i=1:N

       e = sdpvar(6,1);
       S = sdpvar(6,6);
       
       [L1,coeff1] = polynomial(e,monomialOrder);
       [L3,coeff3] = polynomial(e,monomialOrder);
              
       V = e'*P{i}*e;
       Vdot = e'*Q{i}*e ...
              + e'*(P{i+1} - P{i})*e / dt;
       c1 = sos(rhodot(i) - Vdot ...
                - L1*(rho(i) - V));
       c4 = sos(L3);
       c5 = sos(1 - e'*S*e ...
                - L3*(rho(i) - V));
       c8 = S >= 0;
       vars = [coeff1;coeff3;S(:)];

       sol = solvesos([c1 c4 c5 c8],-geomean(S),[],vars);
         
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