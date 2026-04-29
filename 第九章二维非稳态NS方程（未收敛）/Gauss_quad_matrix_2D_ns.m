function result = Gauss_quad_matrix_2D_ns(gauss_point_local,gauss_weight_local,vertices,function_name,uh_local,trial_basis_type_c,trial_basis_der_x_c,trial_basis_der_y_c,trial_basis_type,trial_basis_index,trial_basis_der_x,trial_basis_der_y, test_basis_type,test_basis_index,test_basis_der_x,test_basis_der_y)


trial_val=local_basis_2D(gauss_point_local(:,1),gauss_point_local(:,2),vertices,trial_basis_type,trial_basis_index,trial_basis_der_x,trial_basis_der_y);

test_val=local_basis_2D(gauss_point_local(:,1),gauss_point_local(:,2),vertices,test_basis_type,test_basis_index,test_basis_der_x,test_basis_der_y);

func_val=0;

for alpha=1:6

func_val  =func_val+ feval(function_name,gauss_point_local(:,1),gauss_point_local(:,2),vertices,uh_local,trial_basis_type_c,alpha, trial_basis_der_x_c,trial_basis_der_y_c);

end

result = gauss_weight_local*(trial_val.*test_val.*func_val);




