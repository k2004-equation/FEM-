function A = assemble_matrix_2D(nu,matrix_size, gauss_point, gauss_weight, function_name, P, T, T_trial, T_test, trial_basis_type, trial_basis_der_x,trial_basis_der_y, test_basis_type, test_basis_der_x,test_basis_der_y)

A = sparse(matrix_size(1),matrix_size(2));%稀疏矩阵，节省内存

num_local_trial = size(T_trial,1); %试探函数的基函数个数（一个单元上节点个数）
num_local_test  = size(T_test ,1);%测试函数的基函数个数（一个单元上节点个数）

for n = 1:size(T,2)

    vertices=P(:,T(:,n));  %提取第n个单元上所有网格节点（线性元形式）坐标  
    
    [gauss_weight_local,gauss_point_local]=generate_Gauss_local_rectangle(gauss_weight,gauss_point,vertices);
    
    for alpha = 1:num_local_trial  %Nlb_trial

        for beta = 1:num_local_test %Nlb_test
            
            val = Gauss_quad_matrix_2D(nu,gauss_point_local,gauss_weight_local,vertices,function_name,trial_basis_type,alpha,trial_basis_der_x,trial_basis_der_y, test_basis_type,beta,test_basis_der_x,test_basis_der_y);

            A(T_test(beta,n),T_trial(alpha,n)) = A(T_test(beta,n),T_trial(alpha,n)) + val;
            %组装每个单元基函数到全局刚度矩阵
        end
    end
end

end