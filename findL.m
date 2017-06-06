function [coeff1_,coeff3_,S_] = findL(dt,P,Q,rho,rhodot,initRegion,Kp,Kd)

checkDependency('yalmip');
monomialOrder = 2;
N = length(rho)-1;
% N = 10;
coeff1_ = [];
coeff3_ = [];
S_ = [];
Er = 0.1;
max_ar = 5;

% maxPpv = max(max(P(1:3,4:6)));
% % maxPp = max(max(P(1:3,1:3)));
% maxPv = max(max(P(4:6,4:6)));
maxKp = max(max(Kp));
maxKd = max(max(Kd));

parfor i=1:N
   if i==1
       maxPpv = max(max(P{i}(1:3,4:6)));
       maxPv = max(max(P{i}(4:6,4:6)));
%        maxPpv = P{i}(3,6);
%        maxPv = P{i}(6,6);
       
       e = sdpvar(6,1);
       S = sdpvar(6,6);
       epbar = sdpvar(1,1);
       edbar = sdpvar(1,1);
       epep = e(1:3)'*e(1:3);
       eded = e(4:6)'*e(4:6); 
       
       [L1,coeff1] = polynomial([e;epbar;edbar],monomialOrder);
       [L2,coeff2] = polynomial([e;epbar;edbar],monomialOrder);
       [L3,coeff3] = polynomial([e;epbar;edbar],monomialOrder);

       [L4,coeff4] = polynomial([e;epbar;edbar],monomialOrder);
       [L5,coeff5] = polynomial([e;epbar;edbar],monomialOrder);
       [L6,coeff6] = polynomial([e;epbar;edbar],monomialOrder);
       [L7,coeff7] = polynomial([e;epbar;edbar],monomialOrder);
       
       [L8,coeff8] = polynomial([e;epbar;edbar],monomialOrder);
       [L9,coeff9] = polynomial([e;epbar;edbar],monomialOrder);
       [L10,coeff10] = polynomial([e;epbar;edbar],monomialOrder);
       [L11,coeff11] = polynomial([e;epbar;edbar],monomialOrder);

       [L12,coeff12] = polynomial([e;epbar;edbar],monomialOrder);
       [L13,coeff13] = polynomial([e;epbar;edbar],monomialOrder);
       [L14,coeff14] = polynomial([e;epbar;edbar],monomialOrder);
       [L15,coeff15] = polynomial([e;epbar;edbar],monomialOrder);
       
       V = e'*P{i}*e;
       Vdot = e'*Q{i}*e ...
              + Er*(maxPpv*epbar+maxPv*edbar)*(maxKp*epbar+maxKd*edbar + 9.8 + max_ar) ...
              + e'*(P{i+1} - P{i})*e / dt;
   
       c1 = sos(rhodot(i) - Vdot ...
                - L1*(rho(i) - V) ...
                - L4*(epbar^2 - epep) ...
                - L5*(edbar^2 - eded) ...
                - L6*(edbar) ...
                - L7*(epbar));    
       c2 = sos(L2);
       c3 = sos(rho(i) - V ...
                 - L2*(1-e'*initRegion*e) ...
                 - L8*(epbar^2 - epep) ...
                 - L9*(edbar^2 - eded) ...
                 - L10*(edbar) ...
                 - L11*(epbar));
       c4 = sos(L3);
       c5 = sos(1 - e'*S*e ...
                - L3*(rho(i) - V) ...
                - L12*(epbar^2 - epep) ...
                - L13*(edbar^2 - eded) ...
                - L14*(edbar) ...
                - L15*(epbar));
       c6 = sos(L6);
       c7 = sos(L7);
       c8 = S >= 0;
       c9 = sos(L10);
       c10 = sos(L11);
       c11 = sos(L14);
       c12 = sos(L15);
       
       vars = [coeff1;coeff2;coeff3;coeff4;coeff5;coeff6;coeff7;coeff8;coeff9;coeff10;coeff11;coeff12;coeff13;coeff14;coeff15;S(:)];

       sol = solvesos([c1 c2 c3 c4 c5 c6 c7 c8 c9 c10 c11 c12],-geomean(S),[],vars);

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
       maxPpv = max(max(P{i}(1:3,4:6)));
       maxPv = max(max(P{i}(4:6,4:6)));
%        maxPpv = P{i}(3,6);
%        maxPv = P{i}(6,6);
       
       e = sdpvar(6,1);
       S = sdpvar(6,6);
       epbar = sdpvar(1,1);
       edbar = sdpvar(1,1);
       epep = e(1:3)'*e(1:3);
       eded = e(4:6)'*e(4:6);
       
       [L1,coeff1] = polynomial([e;epbar;edbar],monomialOrder);
       [L3,coeff3] = polynomial([e;epbar;edbar],monomialOrder);

       [L4,coeff4] = polynomial([e;epbar;edbar],monomialOrder);
       [L5,coeff5] = polynomial([e;epbar;edbar],monomialOrder);
       [L6,coeff6] = polynomial([e;epbar;edbar],monomialOrder);
       [L7,coeff7] = polynomial([e;epbar;edbar],monomialOrder);
       
       [L8,coeff8] = polynomial([e;epbar;edbar],monomialOrder);
       [L9,coeff9] = polynomial([e;epbar;edbar],monomialOrder);
       [L10,coeff10] = polynomial([e;epbar;edbar],monomialOrder);
       [L11,coeff11] = polynomial([e;epbar;edbar],monomialOrder); 
       
       V = e'*P{i}*e;
       Vdot = e'*Q{i}*e ...
              + Er*(maxPpv*epbar+maxPv*edbar)*(maxKp*epbar+maxKd*edbar + 9.8 + max_ar) ...
              + e'*(P{i+1} - P{i})*e / dt;
   
       c1 = sos(rhodot(i) - Vdot ...
                - L1*(rho(i) - V) ...
                - L4*(epbar^2 - epep) ...
                - L5*(edbar^2 - eded) ...
                - L6*(edbar) ...
                - L7*(epbar));    
       c4 = sos(L3);
       c5 = sos(1 - e'*S*e ...
                - L3*(rho(i) - V) ... 
                - L8*(epbar^2 - epep) ...
                - L9*(edbar^2 - eded) ...
                - L10*(edbar) ...
                - L11*(epbar));
       c6 = sos(L6);
       c7 = sos(L7);
       c8 = S >= 0;
       c9 = sos(L10);
       c10 = sos(L11);
       
       vars = [coeff1;coeff3;coeff4;coeff5;coeff6;coeff7;coeff8;coeff9;coeff10;coeff11;S(:)];
       sol = solvesos([c1 c4 c5 c6 c7 c8 c9 c10],-geomean(S),[],vars);
         
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