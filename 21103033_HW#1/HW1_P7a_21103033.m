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
k_global= zeros(num_dof, num_dof);   % initialising global stiffness matrix
n_disp=zeros(num_dof,1);   % initialising nodal displacement vector
n_disp(7,1)=-0.3;
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
    k_elem=2*10^5*[ 1  -1;
                   -1   1;];
    K=T'*k_elem*T;
    m= elem_dof(e,:);
    k_global(m,m)= k_global(m,m) + K;
    end
end
Ts=zeros(num_dof,num_dof);
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
        temp=n_disp(j);
        n_disp(j)=n_disp(j-1);
        n_disp(j-1)=temp;
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
% calculation of the nodal displacements in changed global system
n_disp(i_fdof)=(k_global(i_fdof,i_fdof))\(Force(i_fdof)-(k_global(i_fdof,i_edof)*n_disp(i_edof)));
% Again rearragning the nodal disp. vector back to its initial positions
flag5=d;
for i=d:-1:1
    if edof(i)>flag5
    for j= flag5:edof(i)-1
        temp=n_disp(j);
        n_disp(j)=n_disp(j+1);
        n_disp(j+1)=temp;
    end
    end
    flag5=flag5-1;
end
n_disp=Ts*n_disp;
fprintf('\n NODAL DISPLACEMENTS:\n')
fprintf('  NODE         X-AXIS          Y-AXIS\n')
for i=1:num_nodes
    fprintf('% 5d       % 8.3e      % 8.3e\n',i,n_disp(2*i-1),n_disp(2*i))
end
%calculation of member forces
fprintf('\nMEMBER END FORCES:\n')
fprintf(' MEMBER        AXIAL FORCE @ NODE 1       AXIAL FORCE @ NODE 2\n')
for i=1:num_elem
    if(i<4)
    m=elem_dof(i,:);      % getting the global dofs of the concerned member
    n1= elem_conn(i,1);   % 1st node of the concerned element
    n2= elem_conn(i,2);   % 2nd node of the concerned element
    
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
    
    k_elem= ((A(i)*E)/l)*[ 1  -1;     % stiffness matrix for the member is local cordinate system
                           -1  1;];
    
    u_elem= T*n_disp(m,1);    % nodal displacements in the element local coordinates
    
    f_elem=k_elem*u_elem;     % internal axial forces along the length of the member
    fprintf('%5d               % 8.3e                % 8.3e\n',i,f_elem(1),f_elem(2))
    else
    m=elem_dof(i,:);
    T=[ cosd(180) sind(180) 0 0;
         0 0 cosd(180) sind(180);];
    k_elem=2*10^5*[ 1  -1;
                   -1   1;];
    u_elem=T*n_disp(m,1);
    f_elem=k_elem*u_elem;
    fprintf('%5d               % 8.3e                % 8.3e\n',i,f_elem(1),f_elem(2))
    end
end
%========================== END OF PROGRAM ================================
    
    



