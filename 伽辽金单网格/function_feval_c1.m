function feval_c1 = function_feval_c1(num_local_trial,function_name_gu,x,y,vertices,uh,basis_type_c,basis_der_x_c,basis_der_y_c)

r=0;
for alpha=1:num_local_trial

        r=r+uh(alpha).*local_basis_2D(x,y,vertices,basis_type_c,alpha,basis_der_x_c,basis_der_y_c);
end

feval_g=feval(function_name_gu,r);

feval_c1=feval_g;

end