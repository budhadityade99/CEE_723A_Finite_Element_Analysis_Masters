clc;
%HW2_P1c LINEAR SHAPE FUNCTIONS WITH NGP=1
nel=[2 4 10;];    % Array to hold No. of Elements for the system
A=zeros(1,1);
ngp=1;%stores the no.of gauss points
Q=zeros(1,1);% To store the UDL value at gauss pts
zeta=[0 1;];% Gauss pts value and traction point array
zeta2=[-1 -0.5 0 0.5 1];%stores the points where displacements are interpolated
for i=1:length(nel)
    n= nel(i);
    E= 1.3e6;   % Young's Modulus of bar material
    l=40/n;     % Length of each element
    elem_conn=zeros(n,2);   %Initialising a matrix to store the connectivity of each element w.r.t global nodes
    K_global= zeros(n+1,n+1);
    for j=1:n
        elem_conn(j,1)=j;
        elem_conn(j,2)=j+1;
    end
    u=zeros(n+1,1); % Initiating Global Nodal Disp. vector  
    fdof=zeros(n,1);
    for j=1:n
        fdof(j)=fdof(j)+(j+1);
    end
    f=zeros(n+1,1); % Initialising Global external nodal load vector
    q=20;   % Value of peak of distributed load
    for k=1:n
        y1= l*(k-1);
        q1= q*(1-(y1/40));
        y_e1=0.5*l*(1+zeta(1));
        r1= 2.5*10^(-3)*(((y1+y_e1)^2)-80*(y1+y_e1)+2000);
        A1=pi*r1^2;
        A(1,1)=A1;
        K_elem=zeros(2,2);
        B_zeta=[-0.5 0.5];
        for j=1:ngp %iterating over no.of gauss points(here ngp=1)
            K_hat=(B_zeta)'*B_zeta;
            K_elem=K_elem+(((A(1,j)*E*2)/l)*K_hat*2);%weights in this case =2
        end
        p= elem_conn(k,:);      %invoking the global nodal connection of the respective element      
        K_global(p,p)= K_global(p,p)+ K_elem;   %Arranging into global stiffness matrix
        q1_e1=q1-0.5*y_e1;
        Q(1,1)=q1_e1;
        f_elem=zeros(2,1);
        for j=1:ngp
            N_zeta=[0.5*(1-zeta(1,j)) 0.5*(1+zeta(1,j));];
            f_elem=f_elem+(Q(1,j)*(N_zeta)'*0.5*l);
        end
        if(k==n)
            N_zeta=[0.5*(1-zeta(1,2)) 0.5*(1+zeta(1,2));];
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
          N_zeta=[0.5*(1-zeta2(1,j)) 0.5*(1+zeta2(1,j));];
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
              S_g(j,1)=(E*B_zeta*d_e*2)/l;
          end
          cpy=j;
      else
          for j=1:length(zeta2)
              if j==1
                  S_g(cpy,1)=0.5*(((E*B_zeta*d_e*2)/l)+S_g(cpy));
              else
                  S_g(cpy+j-1,1)=(E*B_zeta*d_e*2)/l;
              end
          end
          cpy=cpy+j-1;
      end
    end
    subplot(2,1,2)
    plot(y_g,S_g)
    xlabel('Nodal Height(m)');
    ylabel('Stresses(N/m^2)');
    hold on
end
%========================== END OF PROGRAM ================================