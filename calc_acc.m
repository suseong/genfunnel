function acc = calc_acc(a0,tsq,u)

t1 = tsq(1);
t2 = tsq(2);
t3 = tsq(3);
t4 = tsq(4);
tf = tsq(5);

a1 = a0 + t1*u;
a2 = a1;
a3 = a2 - (t3-t2)*u;
a4 = a3;
a5 = a4 + (tf-t4)*u;

acc = [a0 a1 a2 a3 a4 a5];

end