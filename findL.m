function [coeff1_,coeff3_,coeff_Pe_,S_] = findL(dt,P,Q,rho,rhodot,Kp,Kd)

checkDependency('yalmip');
monomialOrder = 2;
N = length(rho)-1;
coeff1_ = [];
coeff3_ = [];
coeff_Pe_ = [];
S_ = [];
Er = 0.1;
max_ar = 0;

parfor i=1:N
       e = sdpvar(6,1);
       S = sdpvar(6,6);
%        Pe = sdpvar(1,1);
       Ke = sdpvar(1,1);
       
       [L1,coeff1] = polynomial([e;Ke],monomialOrder);
       [L3,coeff3] = polynomial([e;Ke],monomialOrder);
%        [L_Pe,coeff_Pe] = polynomial([e;Pe;Ke],monomialOrder);
       [L_Ke,coeff_Ke] = polynomial([e;Ke],monomialOrder);
                     
       V = e'*P{i}*e;
       Vdot = e'*Q{i}*e ...
              + Er*1*(Ke + max_ar) ...
              + e'*(P{i+1} - P{i})*e / dt;
       c1 = sos(rhodot(i) - Vdot ...
                - L1*(rho(i) - V) ...
                - L_Ke*(Ke - Kp*e(1:3) - Kd*e(4:6)) ...
                );
       c2 = sos(L3);
       c3 = sos(1 - e'*S*e ...
                - L3*(rho(i) - V));
       c4 = S >= 0;
       vars = [coeff1;coeff3;coeff_Ke;S(:)];

       sol = solvesos([c1 c2 c3 c4],-geomean(S),[],vars);
         
       if sol.problem == 0
           coeff1_(:,i) = value(coeff1);
           coeff3_(:,i) = value(coeff3);
%            coeff_Pe_(:,i) = value(coeff_Pe);
%            coeff_Pe_(:,i) = zeros(size(coeff_Pe));
           S_(:,i) = value(S(:));   
       else
           coeff1_(:,i) = zeros(size(coeff1));
           coeff3_(:,i) = zeros(size(coeff3)); 
%            coeff_Pe_(:,i) = zeros(size(coeff_Pe));
           S_(:,i) = zeros(size(S(:)));           
       end
end
coeff_Pe_ = [];
end

%                 - L_Pe*(Pe - P{i}(1:3,4:6)*e(1:3) - P{i}(4:6,4:6)*e(4:6)) ...