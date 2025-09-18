clc;
elem_node=[ 0 0;    % co-ordinates of the nodes of the structure
            5 0;
            10 0;
            15 0;
            20 0;];
       
elem_conn=[ 1 2;    % inter-connectivity of the nodes to determine the 
            2 3;    % members of the structure
            3 4;
            4 5;
            2 6;
            2 7;
            3 8;
            3 9;
            4 10;
            4 11;
            5 12;
            5 13;];  
        
elem_dof=[ 1 2 3 4;   % degrees of freedom for each element in the form
           3 4 5 6;   % U1y U1@ U2y U2@ where y-> vert. dof & @-> rotaional dof
           5 6 7 8;
           7 8 9 10;];
spr_dof=[ 3 11;
          4 12;
          5 13;
          6 14;
          7 15;
          8 16;
          9 17;
          10 18;];     
       
EI=[ 2*10^6 2*10^6 10^6 10^6;];   % Flexural rigidity of the members in N-m^2

l=[ 5 5 5 5;];       % length of members m
k_spr=zeros(size(spr_dof,1),1);   % vector for storing the stiffnesses of the springs
w= 1000;     % value of udl on all members in N/m
P= 5000;   % Value of Point load on Member 2 in N

% Initialising the nodal displacement vector
n_disp=[0;0;-0.01135;-0.0246;-0.1432;0.0126;-0.01;0.0236;0.0400;0.0028;
        0;0;0;0;0;0;0;0;];

% Equivalent Nodal Load Matrix for members in N
Force=[-(w*l(1))/2       -(w*l(2))/2       -(w*l(3))/2       -(w*l(4))/2;    
       -(w*(l(1)^2))/12  -(w*(l(2)^2))/12  -(w*(l(3)^2))/12  -(w*(l(4)^2))/12;   
       -(w*l(1))/2       -(w*l(2))/2       -(w*l(3))/2       -(w*l(4))/2;
       (w*(l(1)^2))/12   (w*(l(2)^2))/12   (w*(l(3)^2))/12   (w*(l(4)^2))/12;];
   
edof=[ 1 2 11 12 13 14 15 16 17 18];     % Array indicating known dofs
d=length(edof);
i_edof=[ 1 2 3 4 5 6 7 8 9 10];   % array to store the indices of edof in the rearranged nodal displacement vector
fdof=[ 3 4 5 6 7 8 9 10 ];    % array to indicate the unknown dofs
i_fdof=[ 11 12 13 14 15 16 17 18];  % similarly for fdof also


%======================== INPUT ENDS HERE ================================

num_nodes= size(elem_node,1);    % no. of nodes in the structure
num_dof= 2* num_nodes+size(spr_dof,1);      % no. of degrees of freedom
num_elem= size(elem_conn,1);        % no. of elements/members 
num_spr= size(spr_dof,1);    % no.of springs attached to the system
F_global=zeros(num_dof,1);   %initialising global nodal load vector
F_global(3,1)=-P;  % assigning the 2nd node with the given point load
a=zeros(length(fdof),num_elem);

flag1=2;       % rearrangement of the nodal disp. vector to account
for i=1:d      % for the boundary conditions
    if(edof(i)>flag1)
        for j=edof(i):-1:flag1
            temp=n_disp(j);
            n_disp(j)=n_disp(j-1);
            n_disp(j-1)=temp;
        end
    end
    flag1=flag1+1;
end
for e=1:num_elem
    if e<=size(elem_dof,1)
    k_elem= ((EI(e))/(l(e)^3))*[  12    (6*l(e))       -12       (6*l(e));   
                               (6*l(e)) (4*(l(e)^2)) -(6*l(e))  (2*(l(e)^2));
                                 -12    -(6*l(e))       12       -(6*l(e));
                               (6*l(e)) (2*(l(e)^2)) -(6*l(e))  (4*(l(e)^2));];
                
    m= elem_dof(e,:);    % array to store the degrees of freedom of the concerned element
    L=zeros(4,num_dof);  % defining gather operator for beam elements
    for j=1:4
        for k=1:num_dof
            if(m(j)==k)
                L(j,k)=1;
            end
        end
    end
    K=L'*k_elem*L;
    F_global(m,1)= F_global(m,1) + Force(:,e);
    flag3=2;     % rearrangement of the rows of K
    for i=1:d
        if(edof(i)>flag3)
            for j=edof(i):-1:flag3
            temp=K(j,:);
            K(j,:)=K(j-1,:);
            K(j-1,:)=temp;
            end
        end
        flag3=flag3+1;
    end
    flag4=2;      % rearrangement of the columns of K
    for i=1:d
        if(edof(i)>flag4)
            for j=edof(i):-1:flag4
            temp=K(:,j);
            K(:,j)=K(:,j-1);
            K(:,j-1)=temp;
            end
        end
        flag4=flag4+1;
    end
    a(:,e)=K(i_fdof,i_fdof)*n_disp(i_fdof,1)+K(i_fdof,i_edof)*n_disp(i_edof,1);
    else
        k_elem=[ 1  -1;
                -1   1;];
        m=spr_dof(e-size(elem_dof,1),:);
        L=zeros(2,num_dof);   % defining gather operator for springs
        for j=1:2
            for k=1:num_dof
                if(m(j)==k)
                L(j,k)=1;
                end
            end
        end
        K=L'*k_elem*L;
        flag3=2;     % rearrangement of the rows of K
        for i=1:d
            if(edof(i)>flag3)
                for j=edof(i):-1:flag3
                    temp=K(j,:);
                    K(j,:)=K(j-1,:);
                    K(j-1,:)=temp;
                end
            end
            flag3=flag3+1;
        end
        flag4=2;      % rearrangement of the columns of K
        for i=1:d
            if(edof(i)>flag4)
                for j=edof(i):-1:flag4
                    temp=K(:,j);
                    K(:,j)=K(:,j-1);
                    K(:,j-1)=temp;
                end
            end
            flag4=flag4+1;
        end
        a(:,e)=K(i_fdof,i_fdof)*n_disp(i_fdof,1)+K(i_fdof,i_edof)*n_disp(i_edof,1);
    end
end
a_dof=a(:,5:12);
a1=zeros(size(fdof,1),1);
for i=1:4
    a1=a1+a(:,i);
end
% Rearrangement of the nodal load vector
flag2=2;
for i=1:d
    if(edof(i)>flag2)
    for j=edof(i):-1:flag2
        temp=F_global(j);
        F_global(j)=F_global(j-1);
        F_global(j-1)=temp;
    end
    end
    flag2=flag2+1;
end
k_spr=a_dof\(F_global(i_fdof,1)-a1);
fprintf(' SL.No        SPRING        STIFFNESS\n')
fprintf('   1           Kt1          %8.3e\n',k_spr(1));
fprintf('   2           Kr1          %8.3e\n',k_spr(2));
fprintf('   3           Kt2          %8.3e\n',k_spr(3));
fprintf('   4           Kr2          %8.3e\n',k_spr(4));
fprintf('   5           Kt3          %8.3e\n',k_spr(5));
fprintf('   6           Kr3          %8.3e\n',k_spr(6));
fprintf('   7           Kt4          %8.3e\n',k_spr(7));
fprintf('   8           Kr4          %8.3e\n',k_spr(8));
%========================== END OF PROGRAM ================================
    
    

