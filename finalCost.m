function [h,dh] = finalCost(t,x)

h = t;
dh = [1,zeros(1,size(x,1))];

end