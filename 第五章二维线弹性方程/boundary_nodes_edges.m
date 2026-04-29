function [boundary_nodes,boundary_edges]=boundary_nodes_edges(basis_type,num_of_element_r,num_of_element_c)
    
if basis_type==202
        basis_element_r=num_of_element_r*2;
        basis_element_c=num_of_element_c*2;
elseif basis_type==201
        basis_element_r=num_of_element_r;
        basis_element_c=num_of_element_c;
end



    nbn=2*(basis_element_r+basis_element_c);
    boundary_nodes=zeros(2,nbn);
    
    %The following boundary condition may change for different problems.
    %All Dirichlet boundary nodes.
    boundary_nodes(1,:)=-1;
    
    %The index in the following is associated with the index in "generate_M_T_triangle.m".
    %bottom boundary nodes.
    for k=1:basis_element_r
        boundary_nodes(2,k)=(k-1)*(basis_element_c+1)+1;
    end
    
    %right boundary nodes.
    for k=basis_element_r+1:basis_element_r+basis_element_c
        boundary_nodes(2,k)=basis_element_r*(basis_element_c+1)+k-basis_element_r;
    end
    
    %top boundary nodes.
    for k=basis_element_r+basis_element_c+1:2*basis_element_r+basis_element_c
        boundary_nodes(2,k)=(2*basis_element_r+basis_element_c+2-k)*(basis_element_c+1);
    end
    
    %left boundary nodes.
    for k=2*basis_element_r+basis_element_c+1:nbn
        boundary_nodes(2,k)=2*basis_element_r+2*basis_element_c+2-k;
    end



        
    nbe=2*(num_of_element_r+num_of_element_c);
    boundary_edges=zeros(4,nbe);
    
    boundary_edges(1,:)=-1;
    
    %The index in the following is associated with the index in "generate_M_T_triangle.m".
    %bottom boundary edges.
    for k=1:num_of_element_r
        boundary_edges(2,k)=(k-1)*2*num_of_element_c+1;
        boundary_edges(3,k)=(k-1)*(num_of_element_c+1)+1;
        boundary_edges(4,k)=k*(num_of_element_c+1)+1;
    end
    
    %right boundary edges.
    for k=basis_element_r+1:num_of_element_r+num_of_element_c
        boundary_edges(2,k)=(num_of_element_r-1)*2*num_of_element_c+2*(k-num_of_element_r);
        boundary_edges(3,k)=num_of_element_r*(num_of_element_c+1)+k-num_of_element_r;
        boundary_edges(4,k)=num_of_element_r*(num_of_element_c+1)+k-num_of_element_r+1;
    end
    
    %top boundary edges.
    for k=num_of_element_r+num_of_element_c+1:2*num_of_element_r+num_of_element_c
        boundary_edges(2,k)=(2*num_of_element_r+num_of_element_c+1-k)*2*num_of_element_c;
        boundary_edges(3,k)=(2*num_of_element_r+num_of_element_c+2-k)*(num_of_element_c+1);
        boundary_edges(4,k)=(2*num_of_element_r+num_of_element_c+1-k)*(num_of_element_c+1);
    end
    
    %left boundary edges.
    for k=2*num_of_element_r+num_of_element_c+1:nbe
        boundary_edges(2,k)=2*(2*num_of_element_r+2*num_of_element_c+1-k)-1;
        boundary_edges(3,k)=2*num_of_element_r+2*num_of_element_c+2-k;
        boundary_edges(4,k)=2*num_of_element_r+2*num_of_element_c+1-k;
    end

end