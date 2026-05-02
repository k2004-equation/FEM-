function result = Gauss_quad_matrix_2D_ns_ns(gauss_point_local,gauss_weight_local,vertices,num_local_trial,function_name,function_name_gu,uh_local,um_local,trial_basis_type_p, trial_basis_der_x_p,trial_basis_der_y_p,trial_basis_type_q, trial_basis_der_x_q,trial_basis_der_y_q,trial_basis_type,trial_basis_index,trial_basis_der_x,trial_basis_der_y, test_basis_type,test_basis_index,test_basis_der_x,test_basis_der_y)


trial_val=local_basis_2D(gauss_point_local(:,1),gauss_point_local(:,2),vertices,trial_basis_type,trial_basis_index,trial_basis_der_x,trial_basis_der_y);

test_val=local_basis_2D(gauss_point_local(:,1),gauss_point_local(:,2),vertices,test_basis_type,test_basis_index,test_basis_der_x,test_basis_der_y);


feval_c= - feval(function_name,num_local_trial,function_name_gu,gauss_point_local(:,1),gauss_point_local(:,2),vertices,uh_local,trial_basis_type_p, trial_basis_der_x_p,trial_basis_der_y_p);

feval_t=0;
for alpha=1:num_local_trial

    feval_t=feval_t+um_local(alpha).*local_basis_2D(gauss_point_local(:,1),gauss_point_local(:,2),vertices,trial_basis_type_q,alpha, trial_basis_der_x_q,trial_basis_der_y_q);
end


result = gauss_weight_local*(trial_val.*test_val.*feval_c.*feval_t);




