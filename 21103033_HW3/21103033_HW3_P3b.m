clc
fid = fopen('Nodes_Circular_Hole.txt') ;  % open the text file
S = textscan(fid,'%s');   % text scan the data
fclose(fid) ;      % close the file
S = S{1} ;
nodes = cellfun(@(x)str2double(x), S);  % convert the cell array to double 
% Remove NaN's which were strings earlier
nodes(isnan(nodes))=[];
col = 3;
count = 0;
temp_arr =[];
temp_row = [];
for i = 1:length(nodes)
    if count == col
        temp_arr = [temp_arr;
                    temp_row];
        count = 0;
        temp_row = [];
    end   
    temp_row = [temp_row,nodes(i)];
    count = count +1;
end
temp_arr = [temp_arr;
            temp_row];
nodes = temp_arr(:,2:end);
nodes=nodes.*0.001;
clear temp_arr temp_row S;
%Creating Element Connectivity Matrix
fid = fopen('Elements_Circular_Hole.txt') ;  % open the text file
S = textscan(fid,'%s');   % text scan the data
fclose(fid) ;      % close the file
S = S{1} ;
elem_conn = cellfun(@(x)str2double(x), S);  % convert the cell array to double 
% Remove NaN's which were strings earlier
elem_conn(isnan(elem_conn))=[];
col = 5;
count = 0;
temp_arr1 =[];
temp_row1 = [];
for i = 1:length(elem_conn)
    if count == col
        temp_arr1 = [temp_arr1;
                    temp_row1];
        count = 0;
        temp_row1 = [];
    end
    temp_row1 = [temp_row1,elem_conn(i)];
    count = count +1;
end
temp_arr1 = [temp_arr1;
            temp_row1];
elem_conn = temp_arr1(:,2:end);
clear temp_arr1 temp_row1 S i;
E=210*10^9;
neu=0.3;
zeta=[-1/sqrt(3) 1/sqrt(3)];
eeta=[-1/sqrt(3) 1/sqrt(3)];
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
            H=0.25*[(n-1)   (1-n)    (1+n)  -1*(1+n)
                    (z-1)  -1*(1+z)  (1+z)   (1-z)];
            J=H*Xe;
            H_hat=J\H;
            Be=[H_hat(1,1) 0 H_hat(1,2) 0 H_hat(1,3) 0 H_hat(1,4) 0
                  0 H_hat(2,1) 0 H_hat(2,2) 0 H_hat(2,3) 0 H_hat(2,4)
                  H_hat(2,1) H_hat(1,1) H_hat(2,2) H_hat(1,2) H_hat(2,3) H_hat(1,3) H_hat(2,4) H_hat(1,4)];
            k_elem=k_elem+Be'*De*Be*det(J);
        end
    end
    elem_dof=[(2*p(1)-1) (2*p(1)) (2*p(2)-1) (2*p(2)) (2*p(3)-1) (2*p(3)) (2*p(4)-1) (2*p(4))];
    K_global(elem_dof,elem_dof)=K_global(elem_dof,elem_dof)+k_elem;
end
dof_id1=zeros(2000,1);
ct1=1;
for j=1:num_nodes
    if(nodes(j,1)==0)
        dof_id1(ct1)=2*j-1;
        ct1=ct1+1;
    end
end
dof_id2=zeros(2000,1);
ct2=1;
for j=1:num_nodes
    if(nodes(j,2)==0)
        dof_id2(ct2)=2*j;
        ct2=ct2+1;
    end
end
dof_id=zeros(1,ct1+ct2-2);
for j=1:ct1-1
    dof_id(1,2*j-1)=dof_id1(j);
    dof_id(1,2*j)=dof_id2(j);
end
A=zeros(length(dof_id),length(n_disp));
A(:,dof_id)=eye(length(dof_id));
ct3=1;
for j=1:num_nodes
    if(nodes(j,1)==0.5)
        trac_id1(ct3)=j;
        ct3=ct3+1;
    end
end
ct4=1;
for j=1:num_nodes
    if(nodes(j,1)==-0.5)
        trac_id2(ct4)=j;
        ct4=ct4+1;
    end
end
% trac_id1(trac_id1==0)=[];
% trac_id2(trac_id2==0)=[];
l=1/(length(trac_id1)-1);
for j=1:length(trac_id1)
    if trac_id1(j)==11 || trac_id1(j)==12
        f_global(2*trac_id1(j)-1,1)=0.5*10^6*l*0.01;
    else
        f_global(2*trac_id1(j)-1,1)=10^6*l*0.01;
    end
end
for j=1:length(trac_id2)
    if (trac_id2(j)==3 || trac_id2(j)==8)
        f_global(2*trac_id2(j)-1,1)=-(0.5*10^6*l*0.01);
    else
        f_global(2*trac_id2(j)-1,1)=-(10^6*l*0.01);
    end
end
aug_K=zeros(size(K_global,1)+length(dof_id),size(K_global,1)+length(dof_id));
aug_K(1:length(n_disp),1:length(n_disp))=K_global;
aug_K([length(n_disp)+1:length(n_disp)+length(dof_id)],1:length(n_disp))=A;
aug_K(1:length(n_disp),[length(n_disp)+1:length(n_disp)+length(dof_id)])=A';
aug_F=zeros(size(K_global,1)+length(dof_id),1);
aug_F(1:length(n_disp),1)=f_global;
n_disp=aug_K\aug_F;
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
u_e=zeros(8,1);
x_e=zeros(length(num_elem)*4,1);
y_e=zeros(length(num_elem)*4,1);
s_e=zeros(length(num_elem)*4,1);
for e=1:num_elem
    p=elem_conn(e,:);
    for i=1:length(p)
        m=p(i);
        u_e(2*i-1)=n_disp(2*m-1);
        u_e(2*i)=n_disp(2*m);
        Xe(i,1)=nodes(m,1);
        Xe(i,2)=nodes(m,2);
    end
    k_elem=zeros(2*length(p),2*length(p));
    r=3;
    for j=1:2
        for k=1:2
            z=zeta(k);
            n=eeta(j);
            N1=0.25*(1-z)*(1-n);
            N2=0.25*(1+z)*(1-n);
            N3=0.25*(1+z)*(1+n);
            N4=0.25*(1-z)*(1+n);
            H=0.25*[(n-1)   (1-n)    (1+n)  -1*(1+n)
                    (z-1)  -1*(1+z)  (1+z)   (1-z)];
            J=H*Xe;
            H_hat=J\H;
            Be=[H_hat(1,1) 0 H_hat(1,2) 0 H_hat(1,3) 0 H_hat(1,4) 0
                  0 H_hat(2,1) 0 H_hat(2,2) 0 H_hat(2,3) 0 H_hat(2,4)
                  H_hat(2,1) H_hat(1,1) H_hat(2,2) H_hat(1,2) H_hat(2,3) H_hat(1,3) H_hat(2,4) H_hat(1,4)];
            S=E*Be*u_e;
            s_e(4*e-r)=S(1,1);
            x_e(4*e-r)=[N1 N2 N3 N4]*Xe(:,1);
            y_e(4*e-r)=[N1 N2 N3 N4]*Xe(:,2);
            r=r-1;
        end
    end
end
scatter3(x_e,y_e,s_e,'+');
xlabel('x');
ylabel('y');
zlabel('Stress(N/m^2)');
disp(max(s_e));