function [coeff1_,coeff3_,S_] = findL2(P,Q,rho,rhodot,unc)

monomialOrder = 2;
N = length(rho)-1;
coeff1_ = [];
coeff3_ = [];
S_ = [];

parfor i=1:N
       e = sdpvar(2,1);
       s = sdpvar(3,1);
       S = [s(1) s(2);
            s(2) s(3)];
       
       epbar = sdpvar(1,1);
       edbar = sdpvar(1,1);
       Ppv = P{i}(1,2);
       Pv =  P{i}(2,2);

       [Lrho1,Crho1] =     polynomial([e;epbar;edbar],monomialOrder);
       [Lrho2,Crho2] =     polynomial( e,monomialOrder);       
       [Lep,Cep] =         polynomial([e;epbar;edbar],monomialOrder);
       [Led,Ced] =         polynomial([e;epbar;edbar],monomialOrder);
       [Lepsign,Cepsign] = polynomial([e;epbar;edbar],monomialOrder);
       [Ledsign,Cedsign] = polynomial([e;epbar;edbar],monomialOrder);
       
       V = e'*P{i}*e;
       Vdot = e'*Q{i}*e + 2*unc*(Ppv*epbar + Pv*edbar);
          
       c1 = sos(rhodot(i) - Vdot ...
                - Lrho1*(rho(i) - V) ...
                - Lep*(epbar^2 - e(1)'*e(1)) ...
                - Led*(edbar^2 - e(2)'*e(2)) ...
                - Lepsign*(epbar) ...
                - Ledsign*(edbar) ...
                );
       c2 = sos(Lepsign);
       c3 = sos(Ledsign);
       c4 = sos(1 - e'*S*e ...
                - Lrho2*(rho(i) - V) ...
                );
       c5 = sos(Lrho2);
       
       c6 = S >= 0;       
       
       constraints = [c1 c2 c3 c4 c5 c6];

       vars = [Crho1;Crho2;Cep;Ced;Cepsign;Cedsign;s];
       
       sol = solvesos(constraints,-geomean(S),[],vars);
         
       if sol.problem == 0
           coeff1_(:,i) = value(Crho1);
           coeff3_(:,i) = value(Crho2);
           S_(:,i) = value(s);   
       else
           coeff1_(:,i) = zeros(size(Crho1));
           coeff3_(:,i) = zeros(size(Crho2)); 
           S_(:,i) = zeros(size(s));           
       end
end

end
