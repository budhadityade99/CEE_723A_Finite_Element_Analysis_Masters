clc;
elem_node=[ 0      0;    % co-ordinates of the nodes of the structure
            0      6;
           3.464   8;
           6.928   6;
           6.928   0;];
       
elem_conn=[ 1 2;    % inter-connectivity of the nodes to determine the 
            2 3;    % members of the structure
            3 4;
            4 5;];  
        
elem_dof=[ 1  2  3  4  5  6;   % degrees of freedom for each element in the form
           4  5  6  7  8  9;   % U1x U1y U1@ U2x U2y U2@ where x-> hor. dof y-> vert. dof 
           7  8  9 10 11 12;   % & @-> rotaional dof
          10 11 12 13 14 15;];    
      
hinge=[ 0 0 0 4 0;];  % array to indicate the location of internal hinges
joint_inc=[ 0 0 0 0 45;];  % array to store the inclination angle of the joints    
I=[ 5.7*10^(-5) 5.5*10^(-5) 5.5*10^(-5) 5.7*10^(-5);];   % MOI of members in m^4
E= 210*10^6;   %Elastic Midulus for all members in KN/m^2
l1=[ 6 4 4 6;];    % length of members in m
A=[ 4.5*10^(-3) 4.0*10^(-3) 4.0*10^(-3) 4.5*10^(-3);];   % C/S area of the members

w= 1;     % value of udl on members 1 and 2 in KN/m

% Equivalent Nodal Load Matrix for members in N
Force=[     0                  0             0      0;
       -(w*l1(1))/2       -(w*l1(2))/2       0      0;    
       -(w*(l1(1)^2))/12  -(w*(l1(2)^2))/12  0      0;   
            0                  0             0      0;
       -(w*l1(1))/2       -(w*l1(2))/2       0      0;
       (w*(l1(1)^2))/12   (w*(l1(2)^2))/12   0      0;];
   
edof=[ 1 2 3 13 ];     % Array indicating known dofs
d=length(edof);
i_edof=[ 1 2 3 4 ];   % array to store the indices of edof in the rearranged nodal displacement vector
fdof=[ 4 5 6 7 8 9 10 11 12 14 ];    % array to indicate the unknown dofs
i_fdof=[ 5 6 7 8 9 10 11 12 13 14 ];  % similarly for fdof also

%======================== INPUT ENDS HERE ================================

num_nodes= size(elem_node,1);    % no. of nodes in the structure
num_dof= 3* num_nodes;          % no. of degrees of freedom
num_elem= size(elem_conn,1);        % no. of elements/members 
k_global= zeros(num_dof-1, num_dof-1);   % initialising global stiffness matrix
n_disp=zeros(num_dof-1,1);   % initialising nodal displacement vector
F_global=zeros(num_dof-1,1);   %initialising global nodal load vector
for e=1:num_elem
    n1=elem_conn(e,1);
    n2=elem_conn(e,2);
    l=l1(e);
    a= (A(e)*E)/l;
    b= (E*I(e))/(l^3);
    
    x1= elem_node(n1,1);
    y1= elem_node(n1,2);
    
    x2= elem_node(n2,1);
    y2= elem_node(n2,2);
    
    C= (x2-x1)/l;
    S= (y2-y1)/l;
    if(n2==hinge(e+1))
        T=[ C    S    0    0    0;
           -S    C    0    0    0;
            0    0    1    0    0;
            0    0    0    C    S;
            0    0    0   -S    C;];            
        k_elem=[ a       0       0       -a       0;
                 0      3*b    3*b*l      0     -3*b;
                 0     3*b*l  3*b*(l^2)   0    -3*b*l;
                -a       0       0        a       0;
                 0     -3*b   -3*b*l      0      3*b;];
        m= elem_dof(e,1:5);
        k_global(m,m)= k_global(m,m) + T'*k_elem*T;
        F_global(m,1)= F_global(m,1) + T'*Force(1:5,e);
    elseif(n1==hinge(e))
        T=[ C     S      0      0      0;
           -S     C      0      0      0;
            0     0      C      S      0;
            0     0     -S      C      0;
            0     0      0      0      1;];
         k_elem=[ a      0      -a       0       0;
                  0     3*b      0     -3*b    3*b*l;
                 -a      0       a       0       0;
                  0    -3*b      0      3*b   -3*b*l;
                  0    3*b*l     0    -3*b*l   3*b*(l^2);];                 
         m= elem_dof(e,1:5);
         k_global(m,m)= k_global(m,m) + T'*k_elem*T;
         F_global(m,1)= F_global(m,1) + T'*Force([1 2 4 5 6],e);
     else
         T=[  C    S    0    0    0    0;
             -S    C    0    0    0    0;
              0    0    1    0    0    0;
              0    0    0    C    S    0;
              0    0    0   -S    C    0;
              0    0    0    0    0    1;];        
         k_elem= [ a       0       0      -a        0        0;
                   0     12*b   6*b*l      0      -12*b    6*b*l;
                   0     6*b*l  4*b*(l^2)  0      -6*b*l   2*b*(l^2);
                  -a       0       0       a        0        0;
                   0    -12*b   -6*b*l     0       12*b    -6*b*l;
                   0     6*b*l  2*b*(l^2)  0      -6*b*l   4*b*(l^2);];          
         m= elem_dof(e,:);    % array to store the degrees of freedom of the current element
         k_global(m,m)= k_global(m,m) + T'*k_elem*T;
         F_global(m,1)= F_global(m,1) + T'*Force(:,e);
    end 
end
% Conversion of the usual global system to a new global system to
% incorporate the inclined supports
Ts=zeros(num_dof-1,num_dof-1);
for i=1:num_nodes
    t=[ cosd(joint_inc(i))  -sind(joint_inc(i))  0;
        sind(joint_inc(i))   cosd(joint_inc(i))  0;
                0                   0            1;];
    Ts([3*i-2 3*i-1 3*i],[3*i-2 3*i-1 3*i])=t;
end
Ts(12,:)=[];
Ts(:,12)=[];
k_global=Ts'*k_global*Ts;
F_global=Ts'*F_global;
% Rearrangement of the nodal displacement vector in the changed global system
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
        temp=F_global(j);
        F_global(j)=F_global(j-1);
        F_global(j-1)=temp;
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
% reconverting the nodal displacement vector from new global system to the original global system
n_disp=Ts*n_disp; 
fprintf('\n NODAL DISPLACEMENTS:')
fprintf('\n   NODE      HORIZONTAL       VERTICAL        ROTATION\n')
f=zeros(6,num_elem);
for e=1:num_elem
    u=zeros(6,1);
    u_elem=zeros(6,1);
    n1=elem_conn(e,1);
    n2=elem_conn(e,2);
    l=l1(e);
    a= (A(e)*E)/l;
    b= (E*I(e))/(l^3);
       
    x1= elem_node(n1,1);
    y1= elem_node(n1,2);
    
    x2= elem_node(n2,1);
    y2= elem_node(n2,2);
    
    C= (x2-x1)/l;
    S= (y2-y1)/l;    
    if(n2==hinge(e+1))
        m=elem_dof(e,1:5);
        u(1:5,1)=n_disp(m,1);
        T=[ C    S    0    0    0;
           -S    C    0    0    0;
            0    0    1    0    0;
            0    0    0    C    S;
            0    0    0   -S    C;];
        k_elem=[ a       0       0       -a       0;
                 0      3*b    3*b*l      0     -3*b;
                 0     3*b*l  3*b*(l^2)   0    -3*b*l;
                -a       0       0        a       0;
                 0     -3*b   -3*b*l      0      3*b;];
        u_elem(1:5,1)=T*u(1:5,1);
        u(6,1)=u(6,1) - ([0 6*b*l 2*b*(l^2) 0 -6*b*l]*u_elem(1:5,1))/(4*b*(l^2));
        f(1:5,e)= f(1:5,e) + (k_elem*u_elem(1:5,1))-Force(1:5,e);
        fprintf('\n MEMBER NO. %i:\n',e);
        for i=1:2
            fprintf('% 5d       % 8.3e      % 8.3e      % 8.3e\n',n1,u(3*i-2)*1000,u(3*i-1)*1000,u(3*i))
            n1=n1+1;
        end
    elseif(n1==hinge(e))
        T=[ C     S      0      0      0;
           -S     C      0      0      0;
            0     0      C      S      0;
            0     0     -S      C      0;
            0     0      0      0      1;];
         k_elem=[ a      0      -a       0       0;
                  0     3*b      0     -3*b    3*b*l;
                 -a      0       a       0       0;
                  0    -3*b      0      3*b   -3*b*l;
                  0    3*b*l     0    -3*b*l   3*b*(l^2);];           
        m= elem_dof(e,1:5);
        u([1 2 4 5 6])=n_disp(m,1);
        u_elem(1:5,1)=u_elem(1:5,1)+T*u([1 2 4 5 6],1);
        u(3,1)=u(3,1)-([0 6*b*l 0 -6*b*l 2*b*(l^2)]*u_elem(1:5,1))/(4*b*(l^2));
        f([1 2 4 5 6],e)=f([1 2 4 5 6],e) + (k_elem*u_elem(1:5,1))-Force([1 2 4 5 6],e);
        fprintf('\n MEMBER NO. %i:\n',e);
        for i=1:2
           fprintf('% 5d       % 8.3e      % 8.3e      % 8.3e\n',n1,u(3*i-2)*1000,u(3*i-1)*1000,u(3*i))
           n1=n1+1;
        end  
    else
         T=[  C    S    0    0    0    0;
             -S    C    0    0    0    0;
              0    0    1    0    0    0;
              0    0    0    C    S    0;
              0    0    0   -S    C    0;
              0    0    0    0    0    1;];
         k_elem= [ a       0       0      -a        0        0;
                   0     12*b   6*b*l      0      -12*b    6*b*l;
                   0     6*b*l  4*b*(l^2)  0      -6*b*l   2*b*(l^2);
                  -a       0       0       a        0        0;
                   0    -12*b   -6*b*l     0       12*b    -6*b*l;
                   0     6*b*l  2*b*(l^2)  0      -6*b*l   4*b*(l^2);];          
         m= elem_dof(e,:);
         u= n_disp(m,1);
         u_elem(:,1)=T*u;
         f(:,e)=f(:,e)+(k_elem*u_elem)-Force(:,e);
         fprintf('\n MEMBER NO. %i:\n',e);
         for i=1:2
           fprintf('% 5d       % 8.3e      % 8.3e      % 8.3e\n',n1,u(3*i-2)*1000,u(3*i-1)*1000,u(3*i))
           n1=n1+1;
         end
    end      
end
fprintf('\n MEMBER END FORCES:\n')
fprintf('   NODE         AXIAL          SHEAR           MOMENT\n')
for e=1:num_elem
  fprintf('\n MEMBER NO.%i:\n',e);
  for i=1:2
    fprintf('%5d       % 8.3e      % 8.3e       % 8.3e\n',i,f(3*i-2,e),f(3*i-1,e),f(3*i,e))
  end
end
%========================== END OF PROGRAM ================================
    
    



