function result = Gauss_quad_matrix_1d(gauss_point_local,gauss_weight_local,vertices,function_name,trial_basis_type,trial_basis_index,trial_basis_der, test_basis_type,test_basis_index,test_basis_der)


%利用高斯节点和权重，经变量替换和雅克比h/2，计算局部单元上的高斯积分值
trial_val = FE_basis_local_fun_1D(gauss_point_local,vertices,trial_basis_type,trial_basis_index,trial_basis_der);
test_val  = FE_basis_local_fun_1D(gauss_point_local,vertices, test_basis_type, test_basis_index, test_basis_der);
 
func_val  = feval(function_name,gauss_point_local);


result = gauss_weight_local*(trial_val.*test_val.*func_val)';
end

