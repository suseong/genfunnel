clear
close
clc

%%
K = diag([10 10 10 4 4 4]);
% P = randn(6,6); P = P'*P/100
P = diag([0.05 0.05 0.05 0.05 0.05 0.05]);
P = P*P;

%%
m = 6.3; % 95 percent confidence
conf = (2 - 2*exp(-m) - 2*m*exp(-m) - m^2*exp(-m))/2

e = [1 0 0 0 0 0]';
chad = e'*inv(P)*e;
e = e / sqrt(chad) * sqrt(6.3*2);

e'*inv(P)*e

% d = 1/sqrt((2*pi)^6*det(P))*exp(-0.5*e'*inv(P)*e);
d = 1/sqrt((2*pi)^6*det(P))*exp(-0.5*12.6);

Q = inv(P)/(-2*log(d) - log((2*pi)^6*det(P)));

eval = eig(inv(Q)*K^2);
plz = sqrt(max(eval))

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
        e'
    end
end

%%
