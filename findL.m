function [coeff1_,coeff2_,coeff3_,S_] = findL(P,Q,rho,rhodot,initRegion)

checkDependency('yalmip');
monomialOrder = 2;
N = length(P)-1;
coeff1_ = [];
coeff2_ = [];
coeff3_ = [];
S_ = [];

parfor i=1:N
   if i==1
       e = sdpvar(6,1);
       S = sdpvar(6,6);
       [L1,coeff1] = polynomial(e,monomialOrder);
       [L2,coeff2] = polynomial(e,monomialOrder);
       [L3,coeff3] = polynomial(e,monomialOrder);
       
       V = e'*P{i}*e;
       Vdot = e'*Q{i}*e;
       
       c1 = sos(rhodot(i) - Vdot - L1*(rho(i) - V));
       c2 = sos(L2);
       c3 = sos(rho(i) -V - L2*(1-e'*initRegion*e));
       c4 = sos(L3);
       c5 = sos(1 - e'*S*e - L3*(rho(i) - V)); 
       c6 = S >= 0;
       
       vars = [coeff1;coeff2;coeff3;S(:)]

       sol = solvesos([c1 c2 c3 c4 c5 c6],-geomean(S),[],vars);

       if sol.problem == 0
           coeff1_(:,i) = value(coeff1);
           coeff2_(:,i) = value(coeff2);
           coeff3_(:,i) = value(coeff3);
           S_(:,i) = value(S(:));
       else
           coeff1_(:,i) = zeros(size(coeff1));
           coeff2_(:,i) = zeros(size(coeff2));
           coeff3_(:,i) = zeros(size(coeff3));           
           S_(:,i) = zeros(size(S(:)));           
       end
   else
       e = sdpvar(6,1);
       S = sdpvar(6,6);
       [L1,coeff1] = polynomial(e,monomialOrder);
       [L3,coeff3] = polynomial(e,monomialOrder);

       V = e'*P{i}*e;
       Vdot = e'*Q{i}*e;

       c1 = sos(rhodot(i) - Vdot - L1*(rho(i) - V));
       c4 = sos(L3);
       c5 = sos(1 - e'*S*e - L3*(rho(i) - V));     
       c6 = S >= 0;
       
       vars = [coeff1;coeff3;S(:)];
       sol = solvesos([c1 c4 c5 c6],-geomean(S),[],vars);
       
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