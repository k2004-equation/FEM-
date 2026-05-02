function result=local_basis_2D(x,y,vertices,basis_type,basis_index,basis_der_x,basis_der_y)

xn1=vertices(1,1);
xn2=vertices(1,2);
xn3=vertices(1,3);
xn4=vertices(1,4);
yn1=vertices(2,1);
yn2=vertices(2,2);
yn3=vertices(2,3);
yn4=vertices(2,4);

h1=xn3-xn1;h2=yn2-yn1;

xh=(2*x-2*xn1-h1)/h1;
yh=(2*y-2*yn1-h2)/h2;

if basis_der_x==0 && basis_der_y==0

    result=reference_basis_2D(xh,yh,basis_type,basis_index,0,0);

elseif basis_der_x==1 && basis_der_y==0

    der_x_result=reference_basis_2D(xh,yh,basis_type,basis_index,1,0);
    

    result=der_x_result*2/h1;
  

elseif basis_der_x==0 && basis_der_y==1

    der_y_result=reference_basis_2D(xh,yh,basis_type,basis_index,0,1);

    result=der_y_result*2/h2;
 

elseif basis_der_x==1 && basis_der_y==1

    der_xy_result=reference_basis_2D(xh,yh,basis_type,basis_index,1,1);

    result=der_xy_result*4/h1/h2;

elseif basis_der_x==2 && basis_der_y==0

    der_xx_result=reference_basis_2D(xh,yh,basis_type,basis_index,2,0);

    result=der_xx_result*4/h1/h1;

elseif basis_der_x==0 && basis_der_y==2

    der_yy_result=reference_basis_2D(xh,yh,basis_type,basis_index,0,2);

    result=der_yy_result*4/h2/h2;

end
