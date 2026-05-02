function uh = assemble_vector_2D_u(s,vector_size,  function_name,ul,trial_basis_type,trial_basis_der_b_x ,trial_basis_der_b_y,P, T,P_trial,T_trial,P_Bi, T_Bi,P_Bitrial,T_Bitrial)

uh = zeros(vector_size,1);

num_local_trial=size(T_Bitrial,1);
ni=1;
for n = 1:size(T_Bi,2)

    if ni*s^2<n
        ni=ni+1;
    end
    vertices =P_Bitrial(:,T_Bitrial(:,n));
    vertices_u =P_trial(:,T_trial(:,ni))';
    uh_local=ul(T_trial(:,ni));
    
    feval_c = feval(function_name,num_local_trial,vertices_u(:,1),vertices_u(:,2),vertices,uh_local,trial_basis_type,trial_basis_der_b_x,trial_basis_der_b_y);          
    
    uh(T_Bitrial(:,n)) = uh(T_Bitrial(:,n)) + feval_c;

end

end
