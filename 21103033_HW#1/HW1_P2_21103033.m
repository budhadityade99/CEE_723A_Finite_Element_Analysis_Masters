clc;
elem_node=[ 1 1 1;    % co-ordinates of the nodes of the structure in
            1 -1 1;   % (x,y,z) format all in 'm'
           -1 0 1;
            0 0 0;];
       
elem_conn=[ 1 4;    % inter-connectivity of the nodes to determine the 
            2 4;    % members of the structure
            3 4;];
        
elem_dof=[ 1 2 3 10 11 12;     % degrees of freedom for each element in the form
           4 5 6 10 11 12;     % u1x,u1y,u1z followed by u2x,u2y,u2z
           7 8 9 10 11 12;];
       
A=1.0;      % Cross- sectional area of the member
E=200;     % Young's Modulus of the Material in N/m^2

Force=zeros(12,1);    % Initialising the nodal load vector
Force(12)=-10;

fdof=[ 10 11 12 ];    % array to indicate the free degrees of freedom

%======================== INPUT ENDS HERE ================================

num_nodes= size(elem_node,1);    % no. of nodes in the structure
num_dof= 3* num_nodes;           % no. of degrees of freedom
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
    C1= (x2-x1)/l;    % calculations of the direction cosines of the element 
    C2= (y2-y1)/l;    % w.r.t the cartesian system
    C3= (z2-z1)/l;
    
    T=[ C1 C2 C3  0  0  0;      % Transformation matrix to map the element level
         0  0  0  C1 C2 C3;];   % parameters to the global equivalents
    
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
fprintf('  NODE         X-AXIS          Y-AXIS         Z-AXIS\n')
for i=1:num_nodes
    fprintf('% 5d       % 8.3e      % 8.3e      % 8.3e\n',i,n_disp(3*i-2),n_disp(3*i-1),n_disp(3*i))
end

%calculation of member forces
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
    C1= (x2-x1)/l;    % calculations of the direction cosines of the element 
    C2= (y2-y1)/l;    % w.r.t the cartesian system
    C3= (z2-z1)/l;
    
    T=[ C1 C2 C3  0  0  0;      % Transformation matrix to map the element level
         0  0  0  C1 C2 C3;];   % parameters to the global equivalents
    
    k_elem= ((A*E)/l)*[ 1  -1;     % stiffness matrix for the member is local cordinate system
                        -1  1;];
    
    u_elem= T*n_disp(m,1);    % nodal displacements in the element local coordinates
    
    f_elem=k_elem*u_elem;     % internal axial forces along the length of the member
    
    fprintf('MEMBER END FORCES:\n')
    fprintf('  MEMBER        AXIAL FORCE @ NODE 1       AXIAL FORCE @ NODE 2\n')
    fprintf('%5d               % 8.3e                % 8.3e\n',i,f_elem(1),f_elem(2))
end


%========================== END OF PROGRAM ================================
    
    



