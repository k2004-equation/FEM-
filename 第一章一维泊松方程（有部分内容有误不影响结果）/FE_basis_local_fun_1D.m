function result = FE_basis_local_fun_1D(x,vertices,basis_type,basis_index,basis_der)

xl = vertices(1);     
xr = vertices(2);   

h = xr - xl;

if basis_type == 101 % 线性元（一次多项式两点插值）
    switch basis_index
        case 1%左端点对应的拉格朗日插值基函数
            switch basis_der
                case 0
                    result = (x - xr)/(xl - xr); 
                case 1
                    result = -1/h*ones(1,length(x));
                case 2
                    result = 0;
                otherwise
                    warning('All derivative >=2 is 0, confirm your input again. Use 2 for higher order.');
            end
        case 2%右端点对应的拉格朗日插值基函数
            switch basis_der
                case 0
                    result = (x - xl)/(xr - xl);
                case 1
                    result = 1/h*ones(1,length(x));
                case 2
                    result = 0;
                otherwise
                    warning('All derivative >=2 is 0, confirm your input again. Use 2 for higher order.');
            end
    end

elseif basis_type == 102 %二次元（二次多项式三点插值）
    xm = (xl+xr)/2;  %中点
    switch basis_index
        case 1  % 左端点的拉格朗日差值基函数
            switch basis_der
                case 0
                    result = (x-xm).*(x-xr)/(xl-xm)/(xl-xr);
                case 1
                    result = (x-xr + x-xm)/(xl-xm)/(xl-xr);
                case 2
                    result = 2/(xl-xm)/(xl-xr)*ones(1,length(x));
                case 3
                    result = 0;
                otherwise
                    warning('All derivative >=3 is 0, confirm your input again. Use 3 for higher order.');
            end
        case 2  %右端点的拉格朗日差值基函数
            switch basis_der
                case 0
                    result = (x-xl).*(x-xm)/(xr-xl)/(xr-xm);
                case 1
                    result = (x-xl + x-xm)/(xr-xl)/(xr-xm);
                case 2
                    result = 2/(xr-xl)/(xr-xm)*ones(1,length(x));
                case 3
                    result = 0;
                otherwise
                    warning('All derivative >=3 is 0, confirm your input again. Use 3 for higher order.');
            end
        case 3  %中点的拉格朗日差值基函数
            switch basis_der
                case 0
                    result = (x-xl).*(x-xr)/(xm-xl)/(xm-xr);
                case 1
                    result = (x-xl + x-xr)/(xm-xl)/(xm-xr);
                case 2
                    result = 2/(xm-xl)/(xm-xr)*ones(1,length(x));
                case 3
                    result = 0;
                otherwise
                    warning('All derivative >=3 is 0, confirm your input again. Use 3 for higher order.');
            end
    end

elseif basis_type == 103  % 三次元（三次多项式的四点插值）
    xb = 2/3*xl + 1/3*xr; 
    xt = 2/3*xr + 1/3*xl;  
    switch basis_index
        case 1  % 左一
            switch basis_der
                case 0
                    result = (x-xb).*(x-xt).*(x-xr)/(xl-xb)/(xl-xt)/(xl-xr);
                case 1
                    result = ((x-xt).*(x-xr) + (x-xb).*(x-xr) + (x-xb).*(x-xt))/(xl-xb)/(xl-xt)/(xl-xr);
                case 2
                    result = ((x-xt)+(x-xr) + (x-xb)+(x-xr) + (x-xb)+(x-xt))/(xl-xb)/(xl-xt)/(xl-xr);
                case 3
                    result = 6/(xl-xb)/(xl-xt)/(xl-xr);
                case 4
                    result = 0;
                otherwise
                    warning('All derivative >=4 is 0, confirm your input again. Use 4 for higher order.');
            end
        case 2  % 右一
            switch basis_der
                case 0
                    result = (x-xl).*(x-xb).*(x-xt)/(xr-xl)/(xr-xb)/(xr-xt);
                case 1
                    result = ((x-xl).*(x-xb) + (x-xb).*(x-xt) + (x-xl).*(x-xt))/(xr-xl)/(xr-xb)/(xr-xt);
                case 2
                    result = ((x-xl)+(x-xb) + (x-xb)+(x-xt) + (x-xl)+(x-xt))/(xr-xl)/(xr-xb)/(xr-xt);
                case 3
                    result = 6/(xr-xl)/(xr-xb)/(xr-xt);
                case 4
                    result = 0;
                otherwise
                    warning('All derivative >=4 is 0, confirm your input again. Use 4 for higher order.');
            end
        case 3  % 左二
            switch basis_der
                case 0
                    result = (x-xl).*(x-xt).*(x-xr)/(xb-xl)/(xb-xt)/(xb-xr);
                case 1
                    result = ((x-xl).*(x-xt) + (x-xt).*(x-xr) + (x-xl).*(x-xr))/(xb-xl)/(xb-xt)/(xb-xr);
                case 2
                    result = ((x-xl)+(x-xt) + (x-xt)+(x-xr) + (x-xl)+(x-xr))/(xb-xl)/(xb-xt)/(xb-xr);
                case 3
                    result = 6/(xb-xl)/(xb-xt)/(xb-xr);
                case 4
                    result = 0;
                otherwise
                    warning('All derivative >=4 is 0, confirm your input again. Use 4 for higher order.');
            end
        case 4  % 右二
            switch basis_der
                case 0
                    result = (x-xl).*(x-xb).*(x-xr)/(xt-xl)/(xt-xb)/(xt-xr);
                case 1
                    result = ((x-xl).*(x-xb) + (x-xb).*(x-xr) + (x-xl).*(x-xr))/(xt-xl)/(xt-xb)/(xt-xr);
                case 2
                    result = ((x-xl)+(x-xb) + (x-xb)+(x-xr) + (x-xl)+(x-xr))/(xt-xl)/(xt-xb)/(xt-xr);
                case 3
                    result = 6/(xt-xl)/(xt-xb)/(xt-xr);
                case 4
                    result = 0;
                otherwise
                    warning('All derivative >=4 is 0, confirm your input again. Use 4 for higher order.');
            end
    end
else
    warining('I dont construct this kind of element.');
end 
end

