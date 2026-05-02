function result = Gauss_quad_vector_2D_ns(gauss_point_local,gauss_weight_local,vertices,function_name,uh_local_u1,uh_local_u2, trial_type_1,trial_type_2,trial_basis_type_c1, trial_basis_der_b_x_c1 ,trial_basis_der_b_y_c1 ,trial_basis_type_c2,trial_basis_der_b_x_c2 ,trial_basis_der_b_y_c2,test_basis_type,test_basis_index,test_basis_der_x,test_basis_der_y)  


test_val=local_basis_2D(gauss_point_local(:,1),gauss_point_local(:,2),vertices,test_basis_type,test_basis_index,test_basis_der_x,test_basis_der_y);

func_val1=0;func_val2=0;


for alpha1=1:6

func_val1  = func_val1 + feval(function_name,gauss_point_local(:,1),gauss_point_local(:,2),uh_local_u1,uh_local_u2, vertices,alpha1,trial_type_1,trial_basis_type_c1, trial_basis_der_b_x_c1 ,trial_basis_der_b_y_c1 );

end

for alpha2=1:6

func_val2  = func_val2 + feval(function_name,gauss_point_local(:,1),gauss_point_local(:,2),uh_local_u1,uh_local_u2, vertices,alpha2,trial_type_2,trial_basis_type_c2,trial_basis_der_b_x_c2 ,trial_basis_der_b_y_c2);

end

result = gauss_weight_local*(test_val.*func_val1.*func_val2);




