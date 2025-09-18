clc
%HW2_P11(a)
%PENALTY FUNCTION METHOD
K=[100  -100   0    0    0    0    0   %Stiffness matrix
   -100  200 -100   0    0    0    0
    0   -100  200 -100   0    0    0
    0     0  -100  200 -100   0    0
    0     0    0  -100  200 -100   0
    0     0    0    0  -100  200 -100
    0     0    0    0    0  -100  100];
dp=zeros(7,1);% stores the nodal displacements computed
f=[1;2;3;4;5;6;7;]; %nodal force vector
d_ex=[0;0.27;0.275;0.25;0.185;0.07;0.14]; %exact nodal displacement given
nc=2;%no of constraints 
Beta=zeros(nc,nc);%Initialising Beta matrix with zeros
A=[1 0 0 0 0 0 0    %A matrix containing info about the constraints
   0 1 0 0 0 -1 0];
b=[0;0.20;]; %stores the RHS vector of the constraint system of equations
k=[2 4 6 8 10 12 14 16 18 20]; %stores the power upto which Beta will be raised
d_err=zeros(7,length(k));% initialising error matrix for all ranges of k
for i=1:length(k)
    n=k(i);
    B=10^n;
    for j=1:nc
        for l=1:nc
            if l==j
                Beta(j,l)=B; %storing diagonal elements of Beta matrix with the Beta values
            end
        end
    end
    K_hat=K+A'*Beta*A; %modified stiffness matrix
    f_hat=f+A'*Beta*b; %modified force vector
    dp=K_hat\f_hat;
    d_err(:,i)=dp;%storing the computed nodal disp values for each value of k
    err=dp-d_ex;%error in calculations
    e_beta=norm(err);%L-2 norm of the error
    scatter(n,e_beta);
    set(gca,'yscale','log')
    xlabel('K=log(10)B');
    ylabel('Norm of the Errors(Log10 scale)');
    grid on
    hold on
end
%In the diagram we saw the norm of error is least for k=10, therefore
%invoking the displacement vector corresponding to k=10 for further
%computations:
f_ex=[-27;26.5;3;4;5;-18.5;7];%exact nodal force vector inc the constraint reactions
dp=d_err(:,5);%extracting the computed disp. closest to the exact solution
fprintf("The computed nodal displacements:\n");
disp(dp);
force=K*dp;%finding the nodal forces inc the reactions for this disp.
e_force=abs(force-f_ex);%finding the error in computed forces
fprintf("The computed nodal forces inc. the reactions:\n");
disp(force);
fprintf("The norm of error in nodal forces computed vs the exact solution(both inc. the reactions)= %5d",norm(e_force));

    