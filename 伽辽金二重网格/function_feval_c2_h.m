function feval_c2 = function_feval_c2_h(num_local_trial,function_name_gudu,x,y,vertices,uh,uh1,basis_type_c,basis_der_x_c,basis_der_y_c)

r1=0;
for alpha1=1:num_local_trial

        r1=r1+uh(alpha1).*local_basis_2D(x,y,vertices,basis_type_c,alpha1,basis_der_x_c,basis_der_y_c);
end
r2=0;
for alpha2=1:num_local_trial

        r2=r2+uh1(alpha2).*local_basis_2D(x,y,vertices,basis_type_c,alpha2,basis_der_x_c,basis_der_y_c);
end

feval_g=feval(function_name_gudu,r1);

feval_c2=feval_g.*(r1-r2);

end