function [isSingular,acc] = singular_check(init,final,input)

isSingular = true;

pos = final(1) - init(1);
pSign = sign(pos);
pos = pSign*pos;

v0 = pSign*init(2); a0 = pSign*init(3);
vf = pSign*final(2); af = pSign*final(3);

checkAcc = 0; checkDec = 0;

if and(v0 > 0, vf > 0)
    % case 1 => acc / dec 
    checkAcc = 1; checkDec = 1;
elseif and(v0 < 0, vf < 0)
    % case 2 => safe
    isSingular = false;
elseif and(v0 > 0, vf < 0)
    checkDec = 1;
    % case 3 => dec
elseif and(v0 < 0, vf > 0)
    % case 4  => acc
    checkAcc = 1;
else
    disp('check init / final conditions')
    keyboard
end




end