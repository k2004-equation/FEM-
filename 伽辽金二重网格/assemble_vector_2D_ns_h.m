function b = assemble_vector_2D_ns_h(vector_size, gauss_point, gauss_weight, function_name,function_name_c,ul,ul1, trial_basis_type,trial_basis_der_b_x ,trial_basis_der_b_y, P, T,P_trial,T_trial, T_test, test_basis_type, test_basis_der_x,test_basis_der_y)

b = sparse(vector_size,1);

num_local_test = size(T_test,1);
num_local_trial=size(T_trial,1);
for n = 1:size(T,2)

    vertices =P(:,T(:,n));
    
    [gauss_weight_local,gauss_point_local]=generate_Gauss_local_triangle(gauss_weight,gauss_point,vertices);
 
        uh_local=ul(T_trial(:,n));
        uh1_local=ul1(T_trial(:,n));
    for beta = 1:num_local_test

        val = Gauss_quad_vector_2D_ns_h(gauss_point_local,gauss_weight_local,vertices,function_name,function_name_c,num_local_trial,uh_local, uh1_local,trial_basis_type,trial_basis_der_b_x ,trial_basis_der_b_y,test_basis_type,beta,test_basis_der_x,test_basis_der_y);           
        
        b(T_test(beta,n)) = b(T_test(beta,n)) + val;
    end
end

end
