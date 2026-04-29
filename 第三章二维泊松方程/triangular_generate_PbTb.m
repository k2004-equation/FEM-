function [Pb,Tb] = triangular_generate_PbTb(number_of_element_r,number_of_element_c,interval_r,interval_c,basis_type)

Nr = number_of_element_r;
Nc = number_of_element_c;
xl=interval_r(1);xr=interval_r(2);
yl=interval_c(1);yr=interval_c(2);
switch basis_type
    case 201   % 线性元
        for j=1:Nr+1
            for i=1:Nc+1
                k=(j-1)*(Nc+1)+i;
                Pb(1,k)=xl+(j-1)*(xr-xl)/Nr;
                Pb(2,k)=yl+(i-1)*(yr-yl)/Nc;
            end
        end
        for j=1:Nr
            for i=1:Nc
                n1=(j-1)*2*Nc+2*i-1;
                Tb(1,n1)=(j-1)*(Nc+1)+i;
                Tb(2,n1)=(j-1+1)*(Nc+1)+i;
                Tb(3,n1)=(j-1)*(Nc+1)+i+1;
                n2=(j-1)*2*Nc+2*i;
                Tb(1,n2)=(j-1)*(Nc+1)+i+1;
                Tb(2,n2)=(j-1+1)*(Nc+1)+i;
                Tb(3,n2)=(j-1+1)*(Nc+1)+i+1;
            end
        end
    case 202   % 二次元
        for j=1:2*Nr+1
            for i=1:2*Nc+1
                k=(j-1)*(2*Nc+1)+i;
                Pb(1,k)=xl+(j-1)*(xr-xl)/Nr/2;
                Pb(2,k)=yl+(i-1)*(yr-yl)/Nc/2;
            end
        end
        for j=1:Nr
            for i=1:Nc
                n1=(j-1)*2*Nc+2*i-1;
                Tb(1,n1)=(2*j-2)*(2*Nc+1)+2*i-1;
                Tb(2,n1)=(2*j)*(2*Nc+1)+2*i-1;
                Tb(3,n1)=(2*j-2)*(2*Nc+1)+2*i+1;
                Tb(4,n1)=(2*j-1)*(2*Nc+1)+2*i-1;
                Tb(5,n1)=(2*j-1)*(2*Nc+1)+2*i;
                Tb(6,n1)=(2*j-2)*(2*Nc+1)+2*i;
                n2=(j-1)*2*Nc+2*i;
                Tb(1,n2)=(2*j-2)*(2*Nc+1)+2*i+1;
                Tb(2,n2)=(2*j)*(2*Nc+1)+2*i-1;
                Tb(3,n2)=(2*j)*(2*Nc+1)+2*i+1;
                Tb(4,n2)=(2*j-1)*(2*Nc+1)+2*i;
                Tb(5,n2)=(2*j)*(2*Nc+1)+2*i;
                Tb(6,n2)=(2*j-1)*(2*Nc+1)+2*i+1;
            end
        end
end
end

