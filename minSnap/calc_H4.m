function H_ = calc_H4(in1)
%CALC_H4
%    H_ = CALC_H4(IN1)

%    This function was generated by the Symbolic Math Toolbox version 6.3.
%    17-Jan-2017 10:35:29

tf = in1(:,2);
to = in1(:,1);
t2 = tf.^2;
t3 = to.^2;
t4 = t2.^2;
t5 = t3.^2;
t6 = t2.*7.2e1;
t7 = t3.*-7.2e1+t6;
t8 = t2.*tf.*1.2e2;
t9 = t8-t3.*to.*1.2e2;
t10 = t4.*3.6e2;
t11 = t5.*-3.6e2+t10;
t12 = t4.*1.8e2;
t13 = t5.*-1.8e2+t12;
t14 = t4.*tf.*5.76e2;
t15 = t14-t5.*to.*5.76e2;
t16 = t2.*t4.*1.2e3;
t17 = t16-t3.*t5.*1.2e3;
t18 = t4.*tf.*2.52e2;
t19 = t18-t5.*to.*2.52e2;
t20 = t2.*t4.*8.4e2;
t21 = t20-t3.*t5.*8.4e2;
t22 = t2.*t4.*tf.*1.8e3;
t23 = t22-t3.*t5.*to.*1.8e3;
t24 = t4.^2;
t25 = t24.*3.15e3;
t26 = t5.^2;
t27 = t25-t26.*3.15e3;
H_ = reshape([0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,tf.*3.6e1-to.*3.6e1,t7,t9,t13,t19,0.0,0.0,0.0,t7,t2.*tf.*1.92e2-t3.*to.*1.92e2,t11,t15,t21,0.0,0.0,0.0,t9,t11,t4.*tf.*7.2e2-t5.*to.*7.2e2,t17,t23,0.0,0.0,0.0,t13,t15,t17,t2.*t4.*tf.*2.057142857142857e3-t3.*t5.*to.*2.057142857142857e3,t27,0.0,0.0,0.0,t19,t21,t23,t27,t24.*tf.*4.9e3-t26.*to.*4.9e3],[8,8]);
