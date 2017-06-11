function [cost,c_x,c_y] = calc_cost4(T,p,opts)

H = zeros(8*4);
for i=1:4
   H(8*(i-1)+1:8*i,8*(i-1)+1:8*i) = calc_H4([T(i) T(i+1)]); 
end
AA = calc_AA4(T);

% opts = optimoptions('quadprog','Display','off');

i=1; aa = [p(i,1);p(i,2);p(i,2);p(i,3);p(i,3);p(i,4);p(i,4);p(i,5);zeros(size(AA,1)-8,1)];
% [c_x,cost_x,exitVal] = quadprog(H,zeros(1,size(H,1)),zeros(size(H)),zeros(size(H,1),1),AA,aa,[],[],[],opts);
[c_x,cost_x] = qpOASES_mex(H,zeros(size(H,1),1),AA,[],[],aa,aa,opts);

i=2; aa = [p(i,1);p(i,2);p(i,2);p(i,3);p(i,3);p(i,4);p(i,4);p(i,5);zeros(size(AA,1)-8,1)];
[c_y,cost_y] = qpOASES_mex(H,zeros(size(H,1),1),AA,[],[],aa,aa,opts);
% [c_y,cost_y,exitVal] = quadprog(H,zeros(1,size(H,1)),zeros(size(H)),zeros(size(H,1),1),AA,aa,[],[],[],opts);

cost = cost_x + cost_y;

end