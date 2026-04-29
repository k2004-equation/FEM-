function [A,b] = dirichlet_boundary_process(Dirichlet_boundary_function_name,A,b,t,boundary_nodes,P_trial)

nbn=size(boundary_nodes,2);

for k=1:nbn

    if boundary_nodes(1,k)==-1 
        i=boundary_nodes(2,k);
        A(i,:)=0;
        A(i,i)=1;
        b(i,1)=feval(Dirichlet_boundary_function_name,P_trial(1,i),P_trial(2,i),t);
    end

end
