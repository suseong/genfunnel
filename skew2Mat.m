function R = skew2Mat(w)

% if size(w,1)*size(w,2) == 3
    R = zeros(3,3);

    R = [0 -w(3) w(2);
         w(3) 0 -w(1);
        -w(2) w(1) 0];
    
% elseif size(w,1)*size(w,2) == 6
% 
%     R = [0 -w(6)   w(5) w(1);
%          w(6) 0   -w(4) w(2);
%         -w(5) w(4) 0    w(3);
%          0    0    0    0];
%     
% end

end