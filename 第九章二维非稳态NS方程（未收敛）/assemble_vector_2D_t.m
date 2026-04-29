function b = assemble_vector_2D_t(nu,vector_size, gauss_point, gauss_weight, function_name,t, P, T, T_test, test_basis_type, test_basis_der_x,test_basis_der_y)

b = sparse(vector_size,1);

num_local_test = size(T_test,1);

for n = 1:size(T,2)

    vertices =P(:,T(:,n));
    
    [gauss_weight_local,gauss_point_local]=generate_Gauss_local_triangle(gauss_weight,gauss_point,vertices);
    
    for beta = 1:num_local_test

        val = Gauss_quad_vector_2D_t(nu,gauss_point_local,gauss_weight_local,vertices,function_name, t,test_basis_type,beta,test_basis_der_x,test_basis_der_y);           
        
        b(T_test(beta,n)) = b(T_test(beta,n)) + val;
    end
end

end
