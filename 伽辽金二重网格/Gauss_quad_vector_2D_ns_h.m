function result = Gauss_quad_vector_2D_ns_h(gauss_point_local,gauss_weight_local,vertices,function_name,function_name_gu,num_local_trial,uh_local ,uh1_local ,trial_basis_type,trial_basis_der_b_x ,trial_basis_der_b_y,test_basis_type,test_basis_index,test_basis_der_x,test_basis_der_y)  


test_val=local_basis_2D(gauss_point_local(:,1),gauss_point_local(:,2),vertices,test_basis_type,test_basis_index,test_basis_der_x,test_basis_der_y);

feval_c=feval(function_name,num_local_trial,function_name_gu,gauss_point_local(:,1),gauss_point_local(:,2),vertices,uh_local,uh1_local ,trial_basis_type,trial_basis_der_b_x,trial_basis_der_b_y);

result = gauss_weight_local*(test_val.*feval_c);

