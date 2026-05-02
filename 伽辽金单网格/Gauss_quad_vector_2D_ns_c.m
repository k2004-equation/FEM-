function result = Gauss_quad_vector_2D_ns_c(gauss_point_local,gauss_weight_local,vertices,function_name,function_name_c,num_local_trial,uh_local ,trial_basis_type,trial_basis_der_b_x ,trial_basis_der_b_y,test_basis_type,test_basis_index,test_basis_der_x,test_basis_der_y)  


test_val=local_basis_2D(gauss_point_local(:,1),gauss_point_local(:,2),vertices,test_basis_type,test_basis_index,test_basis_der_x,test_basis_der_y);

feval_c= - feval(function_name,num_local_trial,gauss_point_local(:,1),gauss_point_local(:,2),vertices,uh_local,trial_basis_type,trial_basis_der_b_x,trial_basis_der_b_y);

feval_c0=feval(function_name_c,gauss_point_local(:,1),gauss_point_local(:,2));

result = gauss_weight_local*(test_val.*feval_c.*feval_c0);

