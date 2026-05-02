function a =error_L2_assemble_2D(u,ul,gauss_point, gauss_weight,P, T, T_trial, T_test, trial_basis_type, trial_basis_der_x,trial_basis_der_y)


a=0;

num_local_trial = size(T_trial,1); 


for n = 1:size(T,2)

    vertices=P(:,T(:,n));  
    
    [gauss_weight_local,gauss_point_local]=generate_Gauss_local_triangle(gauss_weight,gauss_point,vertices);
    
    u_local=u(T_trial(:,n));

    ul_local=ul(T_trial(:,n));
    for alpha = 1:num_local_trial 

        
            
            val = error_L2_Gauss_quad_2D(gauss_point_local,gauss_weight_local,vertices,u_local,ul_local,trial_basis_type,alpha,trial_basis_der_x,trial_basis_der_y);

            a=a+val;
        
    end
end

end