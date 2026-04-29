function r = function_uh(x,y,uh_local_u1,uh_local_u2, vertices,basis_index_c,trial_type,trial_basis_type_c, trial_basis_der_x_c,trial_basis_der_y_c)
k=length(uh_local_u1);
r=zeros(9,1);

if trial_type==1
    
        
    r=uh_local_u1(basis_index_c)*local_basis_2D(x,y,vertices,trial_basis_type_c,basis_index_c,trial_basis_der_x_c,trial_basis_der_y_c);
    
    

elseif trial_type==2
    
    
        
    r=uh_local_u2(basis_index_c)*local_basis_2D(x,y,vertices,trial_basis_type_c,basis_index_c,trial_basis_der_x_c ,trial_basis_der_y_c);
    
    
    



end
end

