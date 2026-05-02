function b = assemble_vector_2D_ns(vector_size, gauss_point, gauss_weight, function_name,ul,trial_type_1,trial_type_2,trial_basis_type_c1, trial_basis_der_b_x_c1 ,trial_basis_der_b_y_c1 ,trial_basis_type_c2,trial_basis_der_b_x_c2 ,trial_basis_der_b_y_c2, P, T,P_trial,T_trial, T_test, test_basis_type, test_basis_der_x,test_basis_der_y)

b = sparse(vector_size,1);


num_local_test = size(T_test,1);

for n = 1:size(T,2)

    vertices =P(:,T(:,n));
    
    [gauss_weight_local,gauss_point_local]=generate_Gauss_local_rectangle(gauss_weight,gauss_point,vertices);
 
        uh_local_u1=ul(T_trial(:,n));
   
        uh_local_u2=ul(length(P_trial)+T_trial(:,n));
    
    for beta = 1:num_local_test

        val = Gauss_quad_vector_2D_ns(gauss_point_local,gauss_weight_local,vertices,function_name,uh_local_u1,uh_local_u2, trial_type_1,trial_type_2,trial_basis_type_c1, trial_basis_der_b_x_c1 ,trial_basis_der_b_y_c1 ,trial_basis_type_c2,trial_basis_der_b_x_c2 ,trial_basis_der_b_y_c2,test_basis_type,beta,test_basis_der_x,test_basis_der_y);           
        
        b(T_test(beta,n)) = b(T_test(beta,n)) + val;
    end
end

end
