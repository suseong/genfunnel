function AA = calc_AA4(in1)
%CALC_AA4
%    AA = CALC_AA4(IN1)

%    This function was generated by the Symbolic Math Toolbox version 6.3.
%    17-Jan-2017 10:26:00

t0 = in1(:,1);
t1 = in1(:,2);
t2 = in1(:,3);
t3 = in1(:,4);
t4 = in1(:,5);
t6 = t0.^2;
t7 = t6.^2;
t8 = t1.^2;
t9 = t8.^2;
t10 = t2.^2;
t11 = t10.^2;
t12 = t3.^2;
t13 = t12.^2;
t14 = t4.^2;
t15 = t14.^2;
t16 = t1.*2.0;
t17 = t8.*3.0;
t18 = t1.*t8.*4.0;
t19 = t9.*5.0;
t20 = t1.*t9.*6.0;
t21 = t8.*t9.*7.0;
t22 = t2.*2.0;
t23 = t10.*3.0;
t24 = t2.*t10.*4.0;
t25 = t11.*5.0;
t26 = t2.*t11.*6.0;
t27 = t10.*t11.*7.0;
t28 = t3.*2.0;
t29 = t12.*3.0;
t30 = t3.*t12.*4.0;
t31 = t13.*5.0;
t32 = t3.*t13.*6.0;
t33 = t12.*t13.*7.0;
t34 = t1.*6.0;
t35 = t8.*1.2e1;
t36 = t1.*t8.*2.0e1;
t37 = t9.*3.0e1;
t38 = t1.*t9.*4.2e1;
t39 = t2.*6.0;
t40 = t10.*1.2e1;
t41 = t2.*t10.*2.0e1;
t42 = t11.*3.0e1;
t43 = t2.*t11.*4.2e1;
t44 = t3.*6.0;
t45 = t12.*1.2e1;
t46 = t3.*t12.*2.0e1;
t47 = t13.*3.0e1;
t48 = t3.*t13.*4.2e1;
t49 = t1.*2.4e1;
t50 = t8.*6.0e1;
t51 = t1.*t8.*1.2e2;
t52 = t9.*2.1e2;
t53 = t2.*2.4e1;
t54 = t10.*6.0e1;
t55 = t2.*t10.*1.2e2;
t56 = t11.*2.1e2;
t57 = t3.*2.4e1;
t58 = t12.*6.0e1;
t59 = t3.*t12.*1.2e2;
t60 = t13.*2.1e2;
AA = reshape([1.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t0,t1,0.0,0.0,0.0,0.0,0.0,0.0,1.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t6,t8,0.0,0.0,0.0,0.0,0.0,0.0,t0.*2.0,t16,0.0,0.0,0.0,2.0,2.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t0.*t6,t1.*t8,0.0,0.0,0.0,0.0,0.0,0.0,t6.*3.0,t17,0.0,0.0,0.0,t0.*6.0,t34,0.0,0.0,0.0,6.0,6.0,0.0,0.0,0.0,t7,t9,0.0,0.0,0.0,0.0,0.0,0.0,t0.*t6.*4.0,t18,0.0,0.0,0.0,t6.*1.2e1,t35,0.0,0.0,0.0,t0.*2.4e1,t49,0.0,0.0,0.0,t0.*t7,t1.*t9,0.0,0.0,0.0,0.0,0.0,0.0,t7.*5.0,t19,0.0,0.0,0.0,t0.*t6.*2.0e1,t36,0.0,0.0,0.0,t6.*6.0e1,t50,0.0,0.0,0.0,t6.*t7,t8.*t9,0.0,0.0,0.0,0.0,0.0,0.0,t0.*t7.*6.0,t20,0.0,0.0,0.0,t7.*3.0e1,t37,0.0,0.0,0.0,t0.*t6.*1.2e2,t51,0.0,0.0,0.0,t0.*t6.*t7,t1.*t8.*t9,0.0,0.0,0.0,0.0,0.0,0.0,t6.*t7.*7.0,t21,0.0,0.0,0.0,t0.*t7.*4.2e1,t38,0.0,0.0,0.0,t7.*2.1e2,t52,0.0,0.0,0.0,0.0,0.0,1.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t1,t2,0.0,0.0,0.0,0.0,0.0,-1.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t8,t10,0.0,0.0,0.0,0.0,0.0,-t16,t22,0.0,0.0,0.0,-2.0,2.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t1.*t8,t2.*t10,0.0,0.0,0.0,0.0,0.0,-t17,t23,0.0,0.0,0.0,-t34,t39,0.0,0.0,0.0,-6.0,6.0,0.0,0.0,0.0,0.0,t9,t11,0.0,0.0,0.0,0.0,0.0,-t18,t24,0.0,0.0,0.0,-t35,t40,0.0,0.0,0.0,-t49,t53,0.0,0.0,0.0,0.0,t1.*t9,t2.*t11,0.0,0.0,0.0,0.0,0.0,-t19,t25,0.0,0.0,0.0,-t36,t41,0.0,0.0,0.0,-t50,t54,0.0,0.0,0.0,0.0,t8.*t9,t10.*t11,0.0,0.0,0.0,0.0,0.0,-t20,t26,0.0,0.0,0.0,-t37,t42,0.0,0.0,0.0,-t51,t55,0.0,0.0,0.0,0.0,t1.*t8.*t9,t2.*t10.*t11,0.0,0.0,0.0,0.0,0.0,-t21,t27,0.0,0.0,0.0,-t38,t43,0.0,0.0,0.0,-t52,t56,0.0,0.0,0.0,0.0,0.0,0.0,1.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t2,t3,0.0,0.0,0.0,0.0,-1.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t10,t12,0.0,0.0,0.0,0.0,-t22,t28,0.0,0.0,0.0,-2.0,2.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t2.*t10,t3.*t12,0.0,0.0,0.0,0.0,-t23,t29,0.0,0.0,0.0,-t39,t44,0.0,0.0,0.0,-6.0,6.0,0.0,0.0,0.0,0.0,0.0,t11,t13,0.0,0.0,0.0,0.0,-t24,t30,0.0,0.0,0.0,-t40,t45,0.0,0.0,0.0,-t53,t57,0.0,0.0,0.0,0.0,0.0,t2.*t11,t3.*t13,0.0,0.0,0.0,0.0,-t25,t31,0.0,0.0,0.0,-t41,t46,0.0,0.0,0.0,-t54,t58,0.0,0.0,0.0,0.0,0.0,t10.*t11,t12.*t13,0.0,0.0,0.0,0.0,-t26,t32,0.0,0.0,0.0,-t42,t47,0.0,0.0,0.0,-t55,t59,0.0,0.0,0.0,0.0,0.0,t2.*t10.*t11,t3.*t12.*t13,0.0,0.0,0.0,0.0,-t27,t33,0.0,0.0,0.0,-t43,t48,0.0,0.0,0.0,-t56,t60,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t3,t4,0.0,0.0,0.0,-1.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t12,t14,0.0,0.0,0.0,-t28,t4.*2.0,0.0,0.0,0.0,-2.0,2.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t3.*t12,t4.*t14,0.0,0.0,0.0,-t29,t14.*3.0,0.0,0.0,0.0,-t44,t4.*6.0,0.0,0.0,0.0,-6.0,6.0,0.0,0.0,0.0,0.0,0.0,0.0,t13,t15,0.0,0.0,0.0,-t30,t4.*t14.*4.0,0.0,0.0,0.0,-t45,t14.*1.2e1,0.0,0.0,0.0,-t57,t4.*2.4e1,0.0,0.0,0.0,0.0,0.0,0.0,t3.*t13,t4.*t15,0.0,0.0,0.0,-t31,t15.*5.0,0.0,0.0,0.0,-t46,t4.*t14.*2.0e1,0.0,0.0,0.0,-t58,t14.*6.0e1,0.0,0.0,0.0,0.0,0.0,0.0,t12.*t13,t14.*t15,0.0,0.0,0.0,-t32,t4.*t15.*6.0,0.0,0.0,0.0,-t47,t15.*3.0e1,0.0,0.0,0.0,-t59,t4.*t14.*1.2e2,0.0,0.0,0.0,0.0,0.0,0.0,t3.*t12.*t13,t4.*t14.*t15,0.0,0.0,0.0,-t33,t14.*t15.*7.0,0.0,0.0,0.0,-t48,t4.*t15.*4.2e1,0.0,0.0,0.0,-t60,t15.*2.1e2],[23,32]);
