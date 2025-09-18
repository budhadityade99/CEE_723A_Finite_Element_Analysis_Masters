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
       
A=2.5*10^(-3);      % Cross- sectional area of the member
E=2*10^7;     % Young's Modulus of the Material in N/m^2

Force=zeros(20,1);    % Initialising the nodal load vector
Force([4 6 8 10],1)=-100;
Force(11,1)=50;

% array to indicate the free degrees of freedom
fdof=[ 3 4 5 6 7 8 9 10 11 13 14 15 16 17 18 19 20 ];   

%======================== INPUT ENDS HERE ================================

num_nodes= size(elem_node,1);    % no. of nodes in the structure
num_dof= 2* num_nodes;           % no. of degrees of freedom
num_elem= size(elem_conn,1);        % no. of elements/members 
k_global= zeros(num_dof, num_dof);   % initialising global stiffness matrix
n_disp=zeros(num_dof,1);   % initialising nodal displacement vector

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
    
    k_elem= ((A*E)/l)*[ 1  -1;     % stiffness matrix for the member is local cordinate system
                        -1  1;];
    
    K=T'*k_elem*T;     % transforming the stiffness matrix from local to global system
                        
    m= elem_dof(e,:);    % variables to store the degrees of freedom of
                         % the concerned element
    
    k_global(m,m)= k_global(m,m) + K; % arranging the element stiffness 
                                      % to corresponding positions in the
                                      % global stiffness matrix
end

% calculation of the nodal displacements
n_disp(fdof)=(k_global(fdof,fdof))\Force(fdof);

fprintf('\n NODAL DISPLACEMENTS:\n')
fprintf('  NODE         X-AXIS          Y-AXIS\n')
for i=1:num_nodes
    fprintf('% 5d       % 8.3e      % 8.3e\n',i,n_disp(2*i-1),n_disp(2*i))
end

%calculation of member forces
fprintf('\nMEMBER END FORCES:\n')
fprintf(' MEMBER        AXIAL FORCE @ NODE 1       AXIAL FORCE @ NODE 2\n')
for i=1:num_elem
    
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
    
    k_elem= ((A*E)/l)*[ 1  -1;     % stiffness matrix for the member is local cordinate system
                        -1  1;];
    
    u_elem= T*n_disp(m,1);    % nodal displacements in the element local coordinates
    
    f_elem=k_elem*u_elem;     % internal axial forces along the length of the member
    fprintf('%5d               % 8.3e                % 8.3e\n',i,f_elem(1),f_elem(2))
end


%========================== END OF PROGRAM ================================
    
    



