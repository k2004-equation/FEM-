function [A,b] = dirichlet_boundary_process(Dirichlet_boundary_function_name1,Dirichlet_boundary_function_name2,Dirichlet_boundary_function_name3,A,b,t,boundary_nodes,P_trial,P)

nbn=size(boundary_nodes,2);
Nb=length(P_trial);

for k=1:nbn

    if boundary_nodes(1,k)==-1 
        i=boundary_nodes(2,k);
        A(i,:)=0;
        A(i,i)=1;
        b(i,1)=feval(Dirichlet_boundary_function_name1,P_trial(1,i),P_trial(2,i),t);
        A(Nb+i,:)=0;
        A(Nb+i,Nb+i)=1;
        b(Nb+i,1)=feval(Dirichlet_boundary_function_name2,P_trial(1,i),P_trial(2,i),t);
        
    end
end
s=2;
A(2*Nb+s,:)=0;
A(2*Nb+s,2*Nb+s)=1;
b(2*Nb+s,1)=feval(Dirichlet_boundary_function_name3,P(1,s),P(2,s),t);


