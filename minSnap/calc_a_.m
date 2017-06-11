function acc = calc_a_(a,t)

x = zeros(size(t));

chad = zeros(8,8);
chad(3,1) = 2;
chad(4,2) = 6;
chad(5,3) = 12;
chad(6,4) = 20;
chad(7,5) = 30;
chad(8,6) = 42;

for i=1:size(t,2)
    t_ = calc_t(t(i),8);
    acc(i) = a'*chad*t_';
end

end