function result = Gauss_quad_vector_1d(gauss_point_local, gauss_weight_local, vertices, function_name, test_basis_type, test_basis_index, test_basis_der)  



test_val  = FE_basis_local_fun_1D(gauss_point_local, vertices, test_basis_type, test_basis_index, test_basis_der);
func_val  = feval(function_name, gauss_point_local);

result = gauss_weight_local*(test_val.*func_val)';
end
