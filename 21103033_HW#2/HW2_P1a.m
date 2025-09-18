clc;
%HW2_P1a QUADRATIC SHAPE FUNCTIONS USING NGP=2
nel=[2 4 10;];    %% Array to hold No. of Elements for the system
A=zeros(1,2);% To store the C/S area at gauss pts
Q=zeros(1,2);% To store the UDL value at gauss pts
zeta=[(-1/sqrt(3)) (1/sqrt(3)) 1;];% Gauss pts value and traction point array
zeta2=[-1 -0.5 0 0.5 1];%stores the points where displacements are interpolated
for i=1:length(nel)
    n= nel(i);
    E= 1.3e6;   % Young's Modulus of bar material
    l=40/n;     % Length of each element
    elem_conn=zeros(n,3);   %Initialising a matrix to store the connectivity of each element w.r.t global nodes
    K_global=zeros(2*n+1,2*n+1);
    for j=1:n
        elem_conn(j,1)=2*j-1;
        elem_conn(j,2)=2*j;
        elem_conn(j,3)=2*j+1;
    end
    u=zeros(2*n+1,1); % Initiating Global Nodal Disp. vector  
    fdof=zeros(2*n,1);
    for j=1:2*n
        fdof(j)=fdof(j)+(j+1);
    end
    f=zeros(2*n+1,1); % Initialising Global external nodal load vector
    q=20;   % Value of peak of distributed load
    for k=1:n
        y1= l*(k-1);
        q1= q*(1-(y1/40));
        y_e1=0.5*l*(1+zeta(1));
        y_e2=0.5*l*(1+zeta(2));
        r1= 2.5*10^(-3)*(((y1+y_e1)^2)-80*(y1+y_e1)+2000);
        r2= 2.5*10^(-3)*(((y1+y_e2)^2)-80*(y1+y_e2)+2000);
        A1=pi*r1^2; A2=pi*r2^2;% C/S areas at the gauss pts
        A(1,1)=A1; A(1,2)=A2;
        K_elem=zeros(3,3);
        for j=1:2 %iterating over no.of gauss points(here ngp=2)
            B_zeta=[zeta(1,j)-0.5  -2*zeta(1,j)  zeta(1,j)+0.5;];
            K_hat=(B_zeta)'*B_zeta;
            K_elem=K_elem+(((A(1,j)*E*2)/l)*K_hat);
        end
        p= elem_conn(k,:);    %invoking the global nodal connection of the respective element      
        K_global(p,p)= K_global(p,p)+ K_elem;   %Arranging into global stiffness matrix
        q1_e1=q1-0.5*y_e1; q1_e2=q1-0.5*y_e2;% UDL values at gauss pts
        Q(1,1)=q1_e1; Q(1,2)=q1_e2;
        f_elem=zeros(3,1);
        for j=1:2
            N_zeta=[0.5*(zeta(1,j)^2-zeta(1,j)) 1-zeta(1,j)^2 0.5*(zeta(1,j)^2+zeta(1,j));];
            f_elem=f_elem+(Q(1,j)*(N_zeta)'*0.5*l);
        end
        if(k==n)
            N_zeta=[0.5*(zeta(1,3)^2-zeta(1,3)) 1-zeta(1,3)^2 0.5*(zeta(1,3)^2+zeta(1,3));];
            f_elem=f_elem+N_zeta'*500;% including traction BC for last element
        end
        f(p,1)=f(p,1)+f_elem;
    end
    u(fdof)=(K_global(fdof,fdof))\f(fdof);  %Finding the nodal displacements
    u_e=zeros(5,1);%stores the interpolated displacements
    S_g=zeros(4*n+1,1);%stores the interpolated stress values
    y_g=zeros(4*n+1,1);%stores the interpolated global y values for the entire bar
    y=zeros(5,1);%stores the global y values corresponding to interpolated disp. on lement basis
    cpy=0;%stores a copy of the stress at the end node of each element to be used for averaging
    for k=1:n
      p=elem_conn(k,:);
      y1= l*(k-1);
      d_e=u(p,1);
      for j=1:length(zeta2)
          N_zeta=[0.5*(zeta2(1,j)^2-zeta2(1,j)) 1-zeta2(1,j)^2 0.5*(zeta2(1,j)^2+zeta2(1,j));];
          u_e(j,1)=N_zeta*d_e;
      end
      for j=1:length(zeta2)
          y_e=0.5*l*(1+zeta2(j));
          y(j,1)=y1+y_e;
      end
      if cpy==0
          y_g(1:5,1)=y;
      else
          y_g(cpy:cpy+5-1,1)=y;
      end
      subplot(2,1,1)
      plot(y,u_e)
      xlabel('Nodal Height(m)');
      ylabel('Axial Displacement(m)');
      hold on
      if k==1
          for j=1:length(zeta2)
              B_zeta=[zeta2(1,j)-0.5  -2*zeta2(1,j)  zeta2(1,j)+0.5;];
              S_g(j,1)=(E*B_zeta*d_e*2)/l;
          end
          cpy=j;
      else
          for j=1:length(zeta2)
              if j==1
                  B_zeta=[zeta2(1,j)-0.5  -2*zeta2(1,j)  zeta2(1,j)+0.5;];
                  S_g(cpy,1)=0.5*(((E*B_zeta*d_e*2)/l)+S_g(cpy));
              else
                  B_zeta=[zeta2(1,j)-0.5  -2*zeta2(1,j)  zeta2(1,j)+0.5;];
                  S_g(cpy+j-1,1)=(E*B_zeta*d_e*2)/l;
              end
          end
          cpy=cpy+j-1;
      end
    end
    subplot(2,1,2)
    plot(y_g,S_g)
    xlabel('Nodal Height(m)');
    ylabel('Axial Stresses(N/m)');
    hold on
end
%========================== END OF PROGRAM ================================
