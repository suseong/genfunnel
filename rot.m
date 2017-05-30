function R = rot(r,p,y)

Rr = [1 0 0;0 cos(r) -sin(r);0 sin(r) cos(r)];
Rp = [cos(p) 0 sin(p);0 1 0;-sin(p) 0 cos(p)];
Ry = [cos(y) -sin(y) 0;sin(y) cos(y) 0;0 0 1];

R = Ry*Rp*Rr;

end