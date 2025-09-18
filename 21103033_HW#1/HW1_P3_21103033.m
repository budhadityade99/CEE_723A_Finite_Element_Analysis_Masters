clc;
nel=[2 4 10 50 100 200 500;];    %% Array to hold No. of Elements for the system
u_20=zeros(length(nel),1);  %vector to store disp at y=20 for differenyt nel
u_40=zeros(length(nel),1);  %vector to store disp at y=40 for differenyt nel
s_0=zeros(length(nel),1);  %vector to store stress at y=0 for differenyt nel
s_20=zeros(length(nel),1);  %vector to store stress at y=20 for differenyt nel
u_analytic20=zeros(length(nel),1);
u_analytic40=zeros(length(nel),1);
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
    f(n+1,1)=500;
    q=20;   % Value of peak of distributed load
    f_eq=zeros(n+1,1);  %Initialising equivalent nodal load vector due to dist. load
    for k=1:n
        y1= l*(k-1);
        y2= l*k;
        q1= q*(1-(y1/40));
        q2= q*(1-(y2/40));
        f_eq(k)= f_eq(k)+ 0.5*(0.5*(q1+q2)*l);      % Adding the equivalent point load at respective node from
        f_eq(k+1)= f_eq(k+1)+ 0.5*(0.5*(q1+q2)*l);  % Distributed load
        r1= 2.5*10^(-3)*((y1^2)-80*y1+2000);
        r2= 2.5*10^(-3)*((y2^2)-80*y2+2000);
        A1= pi*(r1^2);
        A2= pi*(r2^2);
        A= 0.5*(A1+A2);
        p= elem_conn(k,:);      %invoking the global nodal connection of the respective element
        K_elem= ((A*E)/l)*[ 1   -1;     % Initialising element stiffness matrix;
                           -1    1;];       
        K_global(p,p)= K_global(p,p)+ K_elem;   %Arranging into global stiffness matrix
    end

    F= f+f_eq; % Global nodal load vector including equivalent nodal load terms

    u(fdof)=(K_global(fdof,fdof))\F(fdof);  %Finding the nodal displacements
    fprintf('\nFOR NO.OF ELEMENTS= %2d',n)
    fprintf('\nNODAL DISPLACEMENT AT y=20m (in mm): %5d',u((n+2)/2)*1000)
    fprintf('\nNODAL DISPLACEMENT AT y=40m (in mm): %5d\n',u(n+1)*1000)
    str_1=((u(2)-u(1))/l)*E;
    str_2=(u(((n+2)/2)+1)-u((n+2)/2))*E/l;
    str_3=(u((n+2)/2)-u((n+2)/2-1))*E/l;
    str_4=0.5*(str_2+str_3);
    fprintf('NODAL STRESS AT y=0m (in N/m^2): %5d',str_1)
    s_0(i)=str_1;
    fprintf('\nNODAL STRESS AT y=20m (in N/m^2): %5d\n',str_4)
    s_20(i)=str_4;
    u_20(i)=u_20(i)+u((n+2)/2);
    u_40(i)=u_40(i)+u(n+1);
    u_analytic20(i)=0.38;
    u_analytic40(i)=2.018;
end
subplot(2,2,1)
plot(nel,u_20*1000)
title('Displacement at y=20m')
xlabel('No.of elements');
ylabel('Displacement(in mm)');
subplot(2,2,2)
plot(nel,u_40*1000)
title('Displacement at y=40m')
xlabel('No.of elements');
ylabel('Displacement(in mm)');
subplot(2,2,3)
plot(nel,s_0)
title('Stress at y=0m')
xlabel('No.of elements');
ylabel('Stress(in N/m^2)');
subplot(2,2,4)
plot(nel,s_20)
title('Stress at y=20m')
xlabel('No.of elements');
ylabel('Stress(in N/m^2)');
%========================== END OF PROGRAM ================================






    
    
    
    