clc;
elem_node=[ 0 0 0;    % co-ordinates of the nodes of the structure in
            1 0 0;   % (x,y,z) format all in 'm'
            2 0 0;
            3 0 0;
            4 0 0;
            5 0 0;
            1 1 0;
            2 1 0;
            3 1 0;
            4 1 0];
       
elem_conn=[ 1 2;    % inter-connectivity of the nodes to determine the 
            2 3;    % members of the structure
            3 4;
            4 5;
            5 6;
            7 8;
            8 9;
            9 10;
            2 7;
            3 8;
            4 9;
            5 10;
            1 7;
            3 7;
            2 8;
            4 8;
            3 9;
            5 9;
            4 10;
            6 10];
                                                    
elem_dof=[ 1 2 3 4;     % degrees of freedom for each element in the form
           3 4 5 6;     % u1x,u1y followed by u2x,u2y
           5 6 7 8;
           7 8 9 10;
           9 10 11 12;
           13 14 15 16;
           15 16 17 18;
           17 18 19 20;
           3 4 13 14;
           5 6 15 16;
           7 8 17 18;
           9 10 19 20;
           1 2 13 14;
           5 6 13 14;
           3 4 15 16;
           7 8 15 16;
           5 6 17 18;
           9 10 17 18;
           7 8 19 20;
           11 12 19 20;];

% Given nodal displacements in the question in m
n_disp=[0;0;0.005;-0.0478;0.011;-0.0683;0.0183;-0.0683;0.0244;-0.0478;0.0294;
        0;0.0303;-0.0454;0.0204;-0.0672;0.009;-0.0672;-0.0009;-0.0454;];
A=2.5*10^(-3);      % Cross- sectional area of the members in m^2    
Force=zeros(20,1);    % Initialising the nodal load vector in N
Force([4 6 8 10],1)=-100;
Force(11,1)=50;
mem_type=[5 3 4 8];  % array to hold the no.of member in each group
E=zeros(length(mem_type),1);  % array to store the unknown Elastic moduli of the members

% array to indicate the degrees of freedom
edof=[ 1 2 12;];
i_edof=[ 1 2 3];
d=length(edof);
fdof=[ 3 4 5 6 7 8 9 10 11 13 14 15 16 17 18 19 20 ];   
i_fdof=[ 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20];

%======================== INPUT ENDS HERE ================================

num_nodes= size(elem_node,1);    % no. of nodes in the structure
num_dof= 2* num_nodes;           % no. of degrees of freedom
num_elem= size(elem_conn,1);        % no. of elements/members 
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
    
    n1= elem_conn(e,1);   % 1st node of the concerned element
    n2= elem_conn(e,2);   % 2nd node of the concerned element
    
    x1= elem_node(n1,1);  % co-ordinates of the 1st node
    y1= elem_node(n1,2);
    z1= elem_node(n1,3);
    
    x2= elem_node(n2,1);  % co-ordinates of the 2nd node
    y2= elem_node(n2,2);
    z2= elem_node(n2,3);
    
    l=sqrt((x2-x1)^2 + (y2-y1)^2 + (z2-z1)^2);  % length of the element
    C= (x2-x1)/l;    % calculations of the direction cosines of the element 
    S= (y2-y1)/l;    % w.r.t the cartesian system
        
    T=[ C  S  0  0;      % Transformation matrix to map the element level
        0  0  C  S;];   % parameters to the global equivalents
    
    k_elem=(A/l)*[ 1  -1;     % stiffness matrix for the member is local cordinate system
                   -1  1;];                  
    m= elem_dof(e,:);    % variables to store the degrees of freedom of
                         % the concerned element
    L=zeros(4,num_dof);
    for j=1:4
        for k=1:num_dof
            if(m(j)==k)
                L(j,k)=1;
            end
        end
    end
    K=L'*T'*k_elem*T*L;
    flag3=2;     %rearrangement of the rows of K
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
    a(:,e)=K(i_fdof,i_fdof)*n_disp(i_fdof,1);
end
a_dof=zeros(length(fdof),4);   
ctr1=0;ctr2=1;
for i=1:4
    for j=ctr2:mem_type(i)+ctr1
        a_dof(:,i)=a_dof(:,i)+a(:,j);
    end
    ctr1=ctr1+mem_type(i);
    ctr2=ctr1+1;
end
flag2=2;        % rearrangement of nodal load vector
for i=1:d
    if(edof(i)>flag2)
    for j=edof(i):-1:flag2
        temp=Force(j);
        Force(j)=Force(j-1);
        Force(j-1)=temp;
    end
    end
    flag2=flag2+1;
end
% calculation of the unknown elastic moduli of the members
E=E+(a_dof'*a_dof)\a_dof'*Force(i_fdof,1);
fprintf(' SL.No        MEMBER TYPE       ELASTIC MODULUS\n')
fprintf('   1         BOTTOM CHORD         %8.3e\n',E(1));
fprintf('   2          TOP CHORD           %8.3e\n',E(2));
fprintf('   3           VERTICAL           %8.3e\n',E(3));
fprintf('   4           INCLINED           %8.3e\n',E(4));
%========================== END OF PROGRAM ================================
    
    



