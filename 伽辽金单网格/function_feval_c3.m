function feval_c3 = function_feval_c3(num_local_trial,x,y,vertices,uh,basis_type_c,basis_der_x_c,basis_der_y_c)

r=0;
for alpha=1:num_local_trial

        r=r+uh(alpha).*local_basis_2D(x,y,vertices,basis_type_c,alpha,basis_der_x_c,basis_der_y_c);
end

feval_c3=r;

end