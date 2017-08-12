% clear
clc

t = 0.05;
N = 150;

ang = linspace(-pi,pi,30);
tt = linspace(0,t*(N-1),N);

x = ones(30,N)*diag(tt);

for k=1:N
    p = SSS{1}(:,k);
    p1 = [p(1) p(3);p(3) p(5)];
    invp1 = inv(sqrtm(p1));
    
    for j=1:length(ang)
        xx(:,j) = invp1*[cos(ang(j));sin(ang(j))];
    end
    
    y(:,k) = xx(1,:);
    z(:,k) = xx(2,:);
end

%%
figure(100);clf
hold on

surf(x,y,z,'FaceColor',[0.3 0.9 0.9],'FaceAlpha',0.3,'LineStyle','none');
