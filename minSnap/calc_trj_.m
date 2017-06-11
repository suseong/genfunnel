function x = calc_trj_(a,t)

x = zeros(size(t));

for i=1:size(t,2)
    t_ = calc_t(t(i),8);
   x(i) = t_*a;
end

end