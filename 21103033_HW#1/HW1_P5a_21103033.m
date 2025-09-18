clc;
elem_node=[ 0 0;    % co-ordinates of the nodes of the structure
            5 0;
            10 0;
            15 0;
            20 0;];
       
elem_conn=[ 1 2;    % inter-connectivity of the nodes to determine the 
            2 3;    % members of the structure
            3 4;
            4 5;];  
        
elem_dof=[ 1 2 3 4;   % degrees of freedom for each element in the form
           3 4 5 6;   % U1y U1@ U2y U2@ where y-> vert. dof & @-> rotaional dof
           5 6 7 8;
           7 8 9 10;];    
       
EI=[ 2*10^6 2*10^6 10^6 10^6;];   % Flexural rigidity of the members in N-m^2

l=[ 5 5 5 5;];       % length of members m

w= 1000;     % value of udl on all members in N/m
P= 5000;   % Value of Point load on Member 2 in N

% Equivalent Nodal Load Matrix for members in N
Force=[-(w*l(1))/2       -(w*l(2))/2       -(w*l(3))/2       -(w*l(4))/2;    
       -(w*(l(1)^2))/12  -(w*(l(2)^2))/12  -(w*(l(3)^2))/12  -(w*(l(4)^2))/12;   
       -(w*l(1))/2       -(w*l(2))/2       -(w*l(3))/2       -(w*l(4))/2;
       (w*(l(1)^2))/12   (w*(l(2)^2))/12   (w*(l(3)^2))/12   (w*(l(4)^2))/12;];
   
edof=[ 1 2 7 ];     % Array indicating known dofs
d=length(edof);
i_edof=[ 1 2 3 ];   % array to store the indices of edof in the rearranged nodal displacement vector
fdof=[ 3 4 5 6 8 9 10 ];    % array to indicate the unknown dofs
i_fdof=[ 4 5 6 7 8 9 10 ];  % similarly for fdof also


%======================== INPUT ENDS HERE ================================

num_nodes= size(elem_node,1);    % no. of nodes in the structure
num_dof= 2* num_nodes;          % no. of degrees of freedom
num_elem= size(elem_conn,1);        % no. of elements/members 
k_global= zeros(num_dof, num_dof);   % initialising global stiffness matrix
n_disp=zeros(num_dof,1);   % initialising nodal displacement vector
%n_disp(7,1)=-0.01;   % assigning the 7th dof with the given support settlement
F_global=zeros(num_dof,1);   %initialising global nodal load vector
F_global(3,1)=-P;  % assigning the 2nd node with the given point load
for e=1:num_elem
    
    k_elem= ((EI(e))/(l(e)^3))*[  12    (6*l(e))       -12       (6*l(e));   
                               (6*l(e)) (4*(l(e)^2)) -(6*l(e))  (2*(l(e)^2));
                                 -12    -(6*l(e))       12       -(6*l(e));
                               (6*l(e)) (2*(l(e)^2)) -(6*l(e))  (4*(l(e)^2));];
                
    m= elem_dof(e,:);    % variables to store the degrees of freedom of the concerned element
    k_global(m,m)= k_global(m,m) + k_elem;
    F_global(m,1)= F_global(m,1) + Force(:,e);
end

% Rearrangement of the nodal displacement vector
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
% Rearrangemnt of rows of the global stiffness matrix 
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
% Rearrangement pf columns of the global stiffness matrix
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
% calculation of the nodal displacements
n_disp(i_fdof)=(k_global(i_fdof,i_fdof))\(F_global(i_fdof)-(k_global(i_fdof,i_edof)*n_disp(i_edof)));

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
fprintf('\n NODAL DISPLACEMENTS:\n')
fprintf('  NODE      VERTICAL(mm)    ROTATION(rads)\n')
for i=1:num_nodes
    fprintf('% 5d       % 8.3e      % 8.3e\n',i,n_disp(2*i-1)*1000,n_disp(2*i))
end

%calculation of member end forces
fprintf('MEMBER END FORCES:\n')
fprintf('  NODE       Y-FORCE(N)     Z-MOMENT(Nm)\n')

for e=1:num_elem
    
  k_elem= ((EI(e))/(l(e)^3))*[    12    (6*l(e))       -12       (6*l(e));   
                               (6*l(e)) (4*(l(e)^2)) -(6*l(e))  (2*(l(e)^2));
                                 -12    -(6*l(e))       12       -(6*l(e));
                               (6*l(e)) (2*(l(e)^2)) -(6*l(e))  (4*(l(e)^2));];
 
  m= elem_dof(e,:);
  u=zeros(4,1);
  for i=1:4
      u(i)=n_disp(m(i));
  end
  
  f= (k_elem*u) - Force(:,e);
  
  fprintf('  MEMBER:%i\n',e);
  for i=1:2
    fprintf('%5d       % 8.3e     % 8.3e\n',i,f(2*i-1)/1000,f(2*i)/1000)
  end
end
%========================== END OF PROGRAM ================================
    
    



