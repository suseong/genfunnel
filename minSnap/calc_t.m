function output = calc_t(t,n)

chad = [];

for i=0:n-1
    chad = [chad t^i];
end

output = chad;

end