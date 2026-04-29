function result = reference_basis_2D(x,y,basis_type,basis_index,basis_der_x,basis_der_y)
%三角形剖分

if basis_type == 201 % 线性元（分片二元一次多项式三点插值）
    switch basis_index
        case 1
            if basis_der_x==0&&basis_der_y==0

                result=1-x-y;

            elseif basis_der_x==1&&basis_der_y==0

                result=-1*ones(length(x),1);

            elseif basis_der_x==0 && basis_der_y==1

                result=-1*ones(length(x),1);

            else

                result=0*ones(length(x),1);

            end
            
        case 2
            if basis_der_x==0&&basis_der_y==0

                result=x;

            elseif basis_der_x==1&&basis_der_y==0

                result=1*ones(length(x),1);

            elseif basis_der_x==0&&basis_der_y==1

                result=0*ones(length(x),1);

            else 
                result=0*ones(length(x),1);
            end
        case 3
            if basis_der_x==0&&basis_der_y==0

                result=y;

            elseif basis_der_x==1&&basis_der_y==0

                result=0*ones(length(x),1);

            elseif basis_der_x==0&&basis_der_y==1

                result=1*ones(length(x),1);

            else 
                result=0*ones(length(x),1);
            end
    end

elseif basis_type == 202 %二次元（二维二次多项式六点插值）
    
    switch basis_index
        case 1  
            if basis_der_x == 0 && basis_der_y == 0

                result=2*x.^2+2*y.^2+4*x.*y-3*y-3*x+1;

            elseif basis_der_x == 1 && basis_der_y == 0

                result=4*x+4*y-3;

            elseif basis_der_x == 0 && basis_der_y == 1

                result=4*y+4*x-3;

            elseif basis_der_x == 2 && basis_der_y == 0

                result=4*ones(length(x),1);

            elseif basis_der_x == 0 && basis_der_y == 2

                result=4*ones(length(x),1);

            elseif basis_der_x == 1 && basis_der_y == 1

                result=4*ones(length(x),1);

            end
        case 2 
            if basis_der_x == 0 && basis_der_y == 0

                result=2*x.^2-x;

            elseif basis_der_x == 1 && basis_der_y == 0

                result=4*x-1;

            elseif basis_der_x == 0 && basis_der_y == 1

                result=0*ones(length(x),1);

            elseif basis_der_x == 2 && basis_der_y == 0

                result=4*ones(length(x),1);

            elseif basis_der_x == 0 && basis_der_y == 2

                result=0*ones(length(x),1);

            elseif basis_der_x == 1 && basis_der_y == 1

                result=0*ones(length(x),1);

            end
        case 3
            if basis_der_x == 0 && basis_der_y == 0

                result=2.*y.^2-y;

            elseif basis_der_x == 1 && basis_der_y == 0

                result=0*ones(length(x),1);

            elseif basis_der_x == 0 && basis_der_y == 1

                result=4.*y-1;

            elseif basis_der_x == 2 && basis_der_y == 0

                result=0*ones(length(x),1);

            elseif basis_der_x == 0 && basis_der_y == 2

                result=4*ones(length(x),1);

            elseif basis_der_x == 1 && basis_der_y == 1

                result=0*ones(length(x),1);

            end
        case 4
            if basis_der_x == 0 && basis_der_y == 0

                result=-4*x.^2-4*x.*y+4*x;

            elseif basis_der_x == 1 && basis_der_y == 0

                result=-8*x-4*y+4;

            elseif basis_der_x == 0 && basis_der_y == 1

                result=-4*x;
                
            elseif basis_der_x == 2 && basis_der_y == 0

                result=-8*ones(length(x),1);

            elseif basis_der_x == 0 && basis_der_y == 2

                result=0*ones(length(x),1);

            elseif basis_der_x == 1 && basis_der_y == 1

                result=-4*ones(length(x),1);

            end
        case 5
            if basis_der_x == 0 && basis_der_y == 0

                result=4*x.*y;

            elseif basis_der_x == 1 && basis_der_y == 0

                result=4*y;

            elseif basis_der_x == 0 && basis_der_y == 1

                result=4*x;

            elseif basis_der_x == 2 && basis_der_y == 0

                result=0*ones(length(x),1);

            elseif basis_der_x == 0 && basis_der_y == 2

                result=0*ones(length(x),1);

            elseif basis_der_x == 1 && basis_der_y == 1

                result=4*ones(length(x),1);

            end
        case 6
            if basis_der_x == 0 && basis_der_y == 0

                result=-4*y.^2-4*x.*y+4*y;

            elseif basis_der_x == 1 && basis_der_y == 0

                result=-4*y;

            elseif basis_der_x == 0 && basis_der_y == 1

                result=-8*y-4*x+4;

            elseif basis_der_x == 2 && basis_der_y == 0

                result=0*ones(length(x),1);

            elseif basis_der_x == 0 && basis_der_y == 2

                result=-8*ones(length(x),1);

            elseif basis_der_x == 1 && basis_der_y == 1

                result=-4*ones(length(x),1);

            end
    end

end