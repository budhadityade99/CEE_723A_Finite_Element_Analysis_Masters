clc
%HW2_P11(b)
%LAGRANGE MULTIPLIER METHOD
K=[100  -100   0    0    0    0    0   %stiffness matrix
   -100  200 -100   0    0    0    0
    0   -100  200 -100   0    0    0
    0     0  -100  200 -100   0    0
    0     0    0  -100  200 -100   0
    0     0    0    0  -100  200 -100
    0     0    0    0    0  -100  100];
d_ex=[0;0.27;0.275;0.25;0.185;0.07;0.14]; %exact nodal displacemenys given
f_ex=[-27;26.5;3;4;5;-18.5;7];%exact nodal forces
f=[1;2;3;4;5;6;7;];
A=[1 0 0 0 0 0 0
   0 1 0 0 0 -1 0];
b=[0;0.20;];
K_hat=zeros(9,9);
K_hat(1:7,1:7)=K;
K_hat(1:7,8:9)=A';
K_hat(8:9,1:7)=A;
f_hat=zeros(9,1);
f_hat(1:7,1)=f;
f_hat(8:9,1)=b;
d_lambda=K_hat\f_hat;
force=K*d_lambda(1:7);
fprintf("The combined displacement and lambda vector(last two entries are the lambdas)\n");
disp(d_lambda);
fprintf("The computed nodal forces(inc the reactions)\n");
disp(force);
R=force-f;
fprintf("The computed constrained forces:\n");
disp(R);
fprintf("The Lagrange Multipliers:\n");
disp(d_lambda(8:9,1));
% Thus we observe that the computed constrained forces are exactly same as the
% lagrange multipliers which again verifies that the lagrange multipliers
% are the the values of the constraint forces coming onto the system. Also
% the computed nodal displacement and nodal forces are exactly same as the
% actual exact results given.
