function uh = assemble_vector_2D_u(s,vector_size,  function_name,ul,trial_basis_type,trial_basis_der_b_x ,trial_basis_der_b_y,P, T,P_trial,T_trial,P_Bi, T_Bi,P_Bitrial,T_Bitrial)

uh = zeros(vector_size,1);

num_local_trial=size(T_Bitrial,1);
if trial_basis_type==202
    s=s*2;
end
for ni = 1:size(T,2)

    vertices=P(:,T(:,ni));
    uh_local=ul(T_trial(:,ni));
    if mod(ni,2)~=0
        hr=vertices(1,2)-vertices(1,1);
        hc=vertices(2,3)-vertices(2,2);
        k=0;
        for j=1:s+1
            for i=1:s+1-(j-1)
                k=k+1;
                vertices_h(1,k)=vertices(1,1)+(j-1)*hr/s;
                vertices_h(2,k)=vertices(2,2)+(i-1)*hc/s;
            end
        end
    else
        hr=vertices(1,2)-vertices(1,1);
        hc=vertices(2,3)-vertices(2,2);
        k=0;
        for j=1:s+1
            for i=s+1-(j-1):s+1
                k=k+1;
                vertices_h(1,k)=vertices(1,1)+(j-1)*hr/s;
                vertices_h(2,k)=vertices(2,2)+(i-1)*hc/s;
            end
        end
    end

    for g=1:size(vertices_h,2)
        for n=1:size(P_Bitrial,2)
            if vertices_h(:,g)==P_Bitrial(:,n)
                
                T_h(g)=n;break;
            end
        end
    end

    

    vertices_s =vertices_h';

    feval_c = feval(function_name,num_local_trial,vertices_s(:,1),vertices_s(:,2),vertices,uh_local,trial_basis_type,trial_basis_der_b_x,trial_basis_der_b_y);          
    for g=1:size(vertices_h,2)
    if uh(T_h(g))==0
    uh(T_h(g)) = uh(T_h(g)) + feval_c(g);
    end
    end

     % uh(T_h) = uh(T_h) + feval_c;
    

    
end

end
