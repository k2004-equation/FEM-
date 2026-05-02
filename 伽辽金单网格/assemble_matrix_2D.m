function A = assemble_matrix_2D(matrix_size, gauss_point, gauss_weight, function_name, P, T, T_trial, T_test, trial_basis_type, trial_basis_der_x,trial_basis_der_y, test_basis_type, test_basis_der_x,test_basis_der_y)

A = sparse(matrix_size(1),matrix_size(2));

num_local_trial = size(T_trial,1); 
num_local_test  = size(T_test ,1);

for n = 1:size(T,2)

    vertices=P(:,T(:,n));  
    
    [gauss_weight_local,gauss_point_local]=generate_Gauss_local_triangle(gauss_weight,gauss_point,vertices);
    
    for alpha = 1:num_local_trial 

        for beta = 1:num_local_test 
            
            val = Gauss_quad_matrix_2D(gauss_point_local,gauss_weight_local,vertices,function_name,trial_basis_type,alpha,trial_basis_der_x,trial_basis_der_y, test_basis_type,beta,test_basis_der_x,test_basis_der_y);

            A(T_test(beta,n),T_trial(alpha,n)) = A(T_test(beta,n),T_trial(alpha,n)) + val;
            %组装每个单元基函数到全局刚度矩阵
        end
    end
end

end