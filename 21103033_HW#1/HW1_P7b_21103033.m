clc;
elem_node=[ 0 0 0;    % co-ordinates of the nodes of the structure in
            0 1 0;    % (x,y,z) format all in 'm'
            1 1 0;];
       
elem_conn=[ 1 2;    % inter-connectivity of the nodes to determine the 
            2 3;    % members of the structure
            1 3;
            2 4;];
        
node_dof=[1 2;     %degrees of freedom associated with each node
          3 4;
          5 6;
          7 8;];
                                                    
elem_dof=[ 1 2 3 4;     % degrees of freedom for each element in the form
           3 4 5 6;     % u1x,u1y followed by u2x,u2y
           1 2 5 6;
           3 4 7 8;];
       
A=[6*10^(-4) 6*10^(-4) 8.485*10^(-4)];   % Cross- sectional area of the members
E=2.1*10^11;     % Young's Modulus of the Material in N/m^2

Force=zeros(8,1);    % Initialising the nodal load vector
Force(4,1)=1000000;

ang_node=[0 0 45 0;];   % array to store angles at which the nodes are inclined

% arrays to indicate the degrees of freedom
edof=[ 1 2 6 7 8]; %restrained dofs
d=length(edof);
i_edof=[ 1 2 3 4 5];
fdof=[ 3 4 5 ]; %unrestrained dofs
i_fdof=[6 7 8];
%======================== INPUT ENDS HERE ================================
num_nodes= size(elem_conn,1);    % no. of nodes in the structure
num_dof= 2* num_nodes;           % no. of degrees of freedom
num_elem= size(elem_conn,1);        % no. of elements/members 
k_global=sym(zeros(num_dof, num_dof));   % initialising global stiffness matrix
syms u1 u2 u3 u4 u5 u6 u7 u8;
U=sym('u%d',[num_dof 1]); %symbolic nodal disp. vector
G=sym(zeros(length(fdof),1)); %G vector in the form KU-F=0
J=sym(zeros(length(fdof),length(fdof))); %Jacobian matrix
U_star=[0;0;0;-.3;0;-7.126*10^(-4);0.0079;-3.3593*10^(-4);]; %starting guess nodal disp. vector
for e=1:num_elem
    if(e<4)   
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
    
    k_elem= ((A(e)*E)/l)*[ 1  -1;     % stiffness matrix for the member is local cordinate system
                           -1  1;];
    
    K=T'*k_elem*T;     % transforming the stiffness matrix from local to global system
                        
    m= elem_dof(e,:);    % variables to store the degrees of freedom of
                         % the concerned element
    
    k_global(m,m)= k_global(m,m) + K; % arranging the element stiffness 
                                      % to corresponding positions in the
                                      % global stiffness matrix
    else
    T=[ cosd(180) sind(180) 0 0;
         0 0 cosd(180) sind(180);];
    k_elem=(2*10^5+5*10^6*(u3-0.3)^2)*[ 1  -1;
                                       -1   1;];
    K=T'*k_elem*T;
    m= elem_dof(e,:);
    k_global(m,m)= k_global(m,m) + (K);
    end
end
Ts=zeros(num_dof,num_dof); %transformation matrix for support inclination
for k=1:size(node_dof,1)
    ts=[cosd(ang_node(k)) -sind(ang_node(k));
        sind(ang_node(k)) cosd(ang_node(k));];
    Ts(node_dof(k,:),node_dof(k,:))=ts;
end
k_global=Ts'*k_global*Ts;
Force=Ts*Force;
%rearranging the nodal disp. vector in the new global co-ordinate system
flag1=2;
for i=1:d
    if(edof(i)>flag1)
    for j=edof(i):-1:flag1
        temp=U(j);
        U(j)=U(j-1);
        U(j-1)=temp;
    end
    end
    flag1=flag1+1;
end
% Rearrangement of the nodal load vector in the new global system
flag2=2;
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
% Rearrangemnt of rows of the new global stiffness matrix 
flag3=2;
for i=1:d
    if(edof(i)>flag3)
    for j=edof(i):-1:flag3
        temp=k_global(j,:);
        k_global(j,:)=k_global(j-1,:);
        k_global(j-1,:)=temp;
    end
    end
    flag3=flag3+1;
end
% Rearrangement of columns of the new global stiffness matrix
flag4=2;
for i=1:d
    if(edof(i)>flag4)
    for j=edof(i):-1:flag4
        temp=k_global(:,j);
        k_global(:,j)=k_global(:,j-1);
        k_global(:,j-1)=temp;
    end
    end
    flag4=flag4+1;
end
%Initialising The G vector as KU-F=0
G=k_global(i_fdof,i_edof)*U(i_edof,1)+k_global(i_fdof,i_fdof)*U(i_fdof)-Force(i_fdof);
for i=1:length(fdof)
    for j=1:length(fdof)
        J(i,j)=diff(G(i),U(j+5)); %initialising the jacobian matrix
    end
end
err=ones(size(fdof,1),1);
cpy_J=J; %creating a copy of J as it will go on updating itself
cpy_G=G; %creating a copy of G as it will go on updating itself
for i=1:100
    U=U_star;
    J=cpy_J;G=cpy_G;
    J=double(vpa(subs(J,[u1 u2 u6 u7 u8 u3 u4 u5],U')));
    G=double(vpa(subs(G,[u1 u2 u6 u7 u8 u3 u4 u5],U'))); 
    del_U=-(J\G); %increment in the nodal disp.
    U_star(i_fdof)=U_star(i_fdof)+del_U; %updating the nodal disp vector
    err=abs(U-U_star);
end
%disp(k_global);
k_global=subs(k_global,u3,U_star(5));
%Reaction forces the changed global co-ordinate system
R=k_global(i_edof,i_edof)*U_star(i_edof,1)+k_global(i_edof,i_fdof)*U_star(i_fdof,1)-Force(i_edof);
fprintf('\nSUPPORT REACTIONS(in KN):\n');
fprintf('R1x=%8.3e\n',R(1)/1000);
fprintf('R1y=%8.3e\n',R(2)/1000);
fprintf('R3y(in the inclined axes)=%8.3e\n',R(3)/1000);
fprintf('R4x=%8.3e\n',R(4)/1000);
fprintf('R4y=%8.3e\n',R(5)/1000);
% Again rearragning the nodal disp. vector back to its initial positions
flag5=d;
for i=d:-1:1
    if edof(i)>flag5
    for j= flag5:edof(i)-1
        temp=U_star(j);
        U_star(j)=U_star(j+1);
        U_star(j+1)=temp;
    end
    end
    flag5=flag5-1;
end
% reconverting the nodal displacement vector from new global system to the original global system
U_star=Ts*U_star; 
fprintf('\n NODAL DISPLACEMENTS(in mm):')
fprintf('\n   NODE       HORIZONTAL      VERTICAL\n');
for i=1:length(elem_conn)
   fprintf(' % 5d       % 8.3e      % 8.3e\n',i,U_star(2*i-1)*1000,U_star(2*i)*1000);
end
%========================== END OF PROGRAM ================================