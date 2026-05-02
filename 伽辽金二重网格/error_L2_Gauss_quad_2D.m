function result=error_L2_Gauss_quad_2D(gauss_point_local,gauss_weight_local,vertices,u_local,ul_local,trial_basis_type,alpha,trial_basis_der_x,trial_basis_der_y);



trial_val=local_basis_2D(gauss_point_local(:,1),gauss_point_local(:,2),vertices,trial_basis_type,alpha,trial_basis_der_x,trial_basis_der_y);


result = gauss_weight_local*((u_local(alpha)-ul_local(alpha))*trial_val).^2;


