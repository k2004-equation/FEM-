function result = Gauss_quad_matrix_2D(lambda,mu,gauss_point_local,gauss_weight_local,vertices,function_name,trial_basis_type,trial_basis_index,trial_basis_der_x,trial_basis_der_y, test_basis_type,test_basis_index,test_basis_der_x,test_basis_der_y)


trial_val=local_basis_2D(gauss_point_local(:,1),gauss_point_local(:,2),vertices,trial_basis_type,trial_basis_index,trial_basis_der_x,trial_basis_der_y);

test_val=local_basis_2D(gauss_point_local(:,1),gauss_point_local(:,2),vertices,test_basis_type,test_basis_index,test_basis_der_x,test_basis_der_y);

func_val  = feval(function_name,gauss_point_local(:,1),gauss_point_local(:,2),lambda,mu);

result = gauss_weight_local*(trial_val.*test_val*func_val);




