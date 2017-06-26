clear
close
clc

%%
K = diag([10 10 15 4 4 6]);
% P = randn(6,6); P = P'*P/100
P = diag([0.1 0.1 0.1 0.05 0.05 0.05]);
P = P*P;
det(P);
Q = inv(P)/(-2*log(30) - log((2*pi)^6*det(P)));

eval = eig(inv(Q)*K^2);
plz = sqrt(max(eval));

%%
minVal = 100;

for k=1:10000000
    e = 10*(rand(6,1)-0.5*ones(6,1));
    chad = e'*Q*e;
    e = e / sqrt(chad);
    
    chad = plz - norm(K*e);
    
    if(chad < minVal)
        minVal = chad;
        disp(num2str(minVal));
    end
end

