clc
nodes=[ 0 0;
        1 0;
        1 1;
        0 2;
        1 3;
        0 3;];
elem_conn=[ 1 2 3 4
            4 3 5 6];
elem_dof=[1 2 3 4 5 6 7 8;
          7 8 5 6 9 10 11 12;];
E=3*10^7;
neu=0.3;
zeta=[-1/sqrt(3) 1/sqrt(3)];
eeta=[-1/sqrt(3) 1/sqrt(3)];
f_dof=[5 6 7 8 9 10 11 12];
e_dof=[1 2 3 4];
%==========================INPUT ENDS HERE=================================
num_elem=size(elem_conn,1);
num_nodes=size(nodes,1);
K_global=zeros(2*num_nodes,2*num_nodes);
f_global=zeros(2*num_nodes,1);
n_disp=zeros(2*num_nodes,1);
De=(E/(1-neu^2))*[1    neu    0
                 neu    1     0
                  0     0  (1-neu)/2];
Xe=zeros(4,2);
for e=1:num_elem
    p=elem_conn(e,:);
    for i=1:length(p)
        m=p(i);
        Xe(i,1)=nodes(m,1);
        Xe(i,2)=nodes(m,2);
    end
%     plot(Xe(:,1),Xe(:,2),'b')
%     hold on
    k_elem=zeros(2*length(p),2*length(p));
    for j=1:2
        for k=1:2
            z=zeta(k);
            n=eeta(j);
            N1=0.25*(1-z)*(1-n);
            N2=0.25*(1+z)*(1-n);
            N3=0.25*(1+z)*(1+n);
            N4=0.25*(1-z)*(1+n);
            H=0.25*[(n-1)   (1-n)   (1+n)  -1*(1+n)
                    (z-1)  -1*(1+z)  (1+z)   (1-z)];
            J=H*Xe;
            H_hat=J\H;
            Be=[H_hat(1,1) 0 H_hat(1,2) 0 H_hat(1,3) 0 H_hat(1,4) 0
                  0 H_hat(2,1) 0 H_hat(2,2) 0 H_hat(2,3) 0 H_hat(2,4)
                H_hat(2,1) H_hat(1,1) H_hat(2,2) H_hat(1,2) H_hat(2,3) H_hat(1,3) H_hat(2,4) H_hat(1,4)];
            k_elem=k_elem+Be'*De*Be*det(J);
        end
    end
    K_global(elem_dof(e,:),elem_dof(e,:))=K_global(elem_dof(e,:),elem_dof(e,:))+k_elem;
    le=abs(Xe(1,2)-Xe(4,2));
    f_elem=zeros(2*length(p),1);
    for j=1:2
        n=eeta(j);
        Ny=[0.5*(1-n) 0 0 0 0 0 0.5*(1+n) 0];
        y=Ny([1,3,5,7])*Xe(:,2);
        t=30*(1-(y/3));
        f_elem=f_elem+(Ny'*t*le*0.5);
    end
    f_global(elem_dof(e,:),1)=f_global(elem_dof(e,:),1)+f_elem;
end
n_disp(f_dof,1)=K_global(f_dof,f_dof)\f_global(f_dof,1);
fprintf('Nodal Displacements(mm)\n');
fprintf('NODE NO.       X-DISP         Y-DISP\n');
for j=1:length(nodes)
    fprintf('%5d        %8.3e       %8.3e\n',j,n_disp(2*j-1)*1000,n_disp(2*j)*1000);
end
node_disp=zeros(num_nodes,2);
for j=1:num_nodes
    node_disp(j,1)=n_disp(2*j-1,1);
    node_disp(j,2)=n_disp(2*j);
end
node_disp=node_disp+nodes;
for j=1:num_elem
    p=elem_conn(j,:);
    for i=1:length(p)
        m=p(i);
        Xe(i,1)=node_disp(m,1);
        Xe(i,2)=node_disp(m,2);
    end
%     plot(Xe(:,1),Xe(:,2),'r');
%     hold on
end            
    
            