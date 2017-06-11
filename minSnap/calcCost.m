function cost = calcCost(A,x,coeffs)

chad1 = 0;
for i=1:size(x,2)-1
   chad1 = chad1 + 2*norm(x(:,i+1)-x(:,i))^2; 
end

chad2 = 0;
for i=2:size(x,2)-1
    B = x(:,i);
    [isNotFound,simplex] = isFound(A,B);
    depth = 0; 
    if(~isNotFound)
        depth = calcDepth(A,B,simplex);
    end
    chad2 = chad2 + depth;
end

chad3 = 0;
for i=1:size(x,2)-1
    B = [x(:,i+1) x(:,i)];
    [isNotFound,simplex] = isFound(A,B);
    depth = 0; 
    if(~isNotFound)
        depth = calcDepth(A,B,simplex);
    end
    chad3 = chad3 + depth;
end

cost = coeffs*[chad1;chad2;chad3];

end