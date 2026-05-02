tic

[P_Bi,T_Bi] = triangular_generate_PbTb(num_of_element_r^2,num_of_element_c^2,interval_r,interval_c, basis_type);

[P_Bitrial,T_Bitrial] = triangular_generate_PbTb(num_of_element_r^2,num_of_element_c^2,interval_r,interval_c, trial_basis_type);

vector_size_u=length(P_Bitrial);


for n=1:Nt
    for i=1:Nt
    uL(:,(n-1)*Nt+i)=ul(:,n)+(ul(:,n+1)-ul(:,n))/Nt*(i-1);
    end
end
uL(:,Nt^2+1)=ul(:,Nt+1);

for n=1:Nt^2+1

uh(:,n) = assemble_vector_2D_u(s,vector_size_u, 'function_Hh',uL(:,n),trial_basis_type,trial_basis_der_b_x ,trial_basis_der_b_y, P,T,P_trial,T_trial,P_Bi, T_Bi,P_Bitrial,T_Bitrial);

end

toc