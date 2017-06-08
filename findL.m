function [coeff1_,coeff3_,S_] = findL(dt,P,Q,rho,rhodot,Kp,Kd)

checkDependency('yalmip');
monomialOrder = 2;
N = length(rho)-1;
coeff1_ = [];
coeff3_ = [];
coeff_Pe_ = [];
S_ = [];
Er = 0.0;
ar = 15;
maxKp = max(max(Kp));
maxKd = max(max(Kd));

parfor i=1:N
       e = sdpvar(6,1);
       S = sdpvar(6,6);
       epbar = sdpvar(1,1);
       edbar = sdpvar(1,1);
       Ppv = sdpvar(1,1);
       Pv = sdpvar(1,1);
       x = sdpvar(1,1);

       [Lrho1,Crho1] =     polynomial([e;epbar;edbar],monomialOrder);
       [Lrho2,Crho2] =     polynomial([e;epbar;edbar],monomialOrder);       
       [Lep,Cep] =         polynomial([e;epbar;edbar],monomialOrder);
       [Led,Ced] =         polynomial([e;epbar;edbar],monomialOrder);
       [Lepsign,Cepsign] = polynomial([e;epbar;edbar],monomialOrder);
       [Ledsign,Cedsign] = polynomial([e;epbar;edbar],monomialOrder);
       [LPpv1,CPpv1] =     polynomial(x,monomialOrder);
       [LPpv2,CPpv2] =     polynomial(x,monomialOrder);
       [LPv1,CPv1] =       polynomial(x,monomialOrder);
       [LPv2,CPv2] =       polynomial(x,monomialOrder);
       
       V = e'*P{i}*e;
       Vdot = e'*Q{i}*e ...
              + Er*(maxKp*epbar + maxKd*edbar + ar)*(Ppv*epbar + Pv*edbar) ...
              + e'*(P{i+1} - P{i})*e / dt;
       
       c1 = sos(rhodot(i) - Vdot ...
                - Lrho1*(rho(i) - V) ...
                - Lep*(epbar^2 - e(1:3)'*e(1:3)) ...
                - Led*(edbar^2 - e(4:6)'*e(4:6)) ...
                - Lepsign*(epbar) ...
                - Ledsign*(edbar) ...
                );
       c2 = sos(Lepsign);
       c3 = sos(Ledsign);
       c4 = sos(1 - e'*S*e ...
                - Lrho2*(rho(i) - V) ...
                );
       c5 = sos(Lrho2);
       
       c6 = sos(Ppv - P{i}(1,4) - LPpv1);
       c7 = sos(Ppv + P{i}(1,4) - LPpv1);
       c8 = sos(Ppv - P{i}(3,6) - LPpv2);
       c9 = sos(Ppv + P{i}(3,6) - LPpv2);
       c10 = sos(LPpv1);
       c11 = sos(LPpv2);
    
       c12 = sos(Pv - P{i}(4,4) - LPv1);
       c13 = sos(Pv + P{i}(4,4) - LPv1);
       c14 = sos(Pv - P{i}(6,6) - LPv2);
       c15 = sos(Pv + P{i}(6,6) - LPv2);
       c16 = sos(LPv1);
       c17 = sos(LPv2);
       
       c18 = S >= 0;       
       
       constraints = [c1 c2 c3 c4 c5 c6 c7 c8 c9 c10 ...
                      c11 c12 c13 c14 c15 c16 c17 c18];

       vars = [Crho1;Crho2;Cep;Ced;Cepsign;Cedsign];
       vars = [vars;CPpv1;CPpv2;CPv1;CPv2];
       vars = [vars;Ppv;Pv;S(:)];           
           
       sol = solvesos(constraints,-geomean(S),[],vars);
         
       if sol.problem == 0
           coeff1_(:,i) = value(Crho1);
           coeff3_(:,i) = value(Crho2);
           S_(:,i) = value(S(:));   
       else
           coeff1_(:,i) = zeros(size(Crho1));
           coeff3_(:,i) = zeros(size(Crho2)); 
           S_(:,i) = zeros(size(S(:)));           
       end
end

end

