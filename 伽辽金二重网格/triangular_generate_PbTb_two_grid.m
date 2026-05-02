function [Pb,Tb] = triangular_generate_PbTb_two_grid(number_of_element_r,number_of_element_c,interval_r,interval_c,basis_type)

Nr = number_of_element_r;
Nc = number_of_element_c;
Nhr = number_of_element_r^2;
Nhc = number_of_element_c^2;

xl=interval_r(1);xr=interval_r(2);
yl=interval_c(1);yr=interval_c(2);

switch basis_type
    case 201   % 线性元
        for j=1:Nr
            for i=1:Nc
                k=(j-1)*(Nc)+i;
                tx(1,k)=xl+(j-1)*(xr-xl)/Nr;
                ty(2,k)=yl+(i-1)*(yr-yl)/Nc;
                for jh=1:Nr+1
                    for ih=1:Nc+1
                        kh=(k-1)*(Nr+1)*(Nc+1)+(jh-1)*(Nc+1)+ih;
                        Pb(1,kh)=tx(1,k)+(jh-1)*(xr-xl)/Nr/Nr;
                        Pb(2,kh)=ty(2,k)+(ih-1)*(yr-yl)/Nc/Nc;
                    end
                end
            end
        end
       

        for j=1:Nr
            for i=1:Nc

                kh=(Nr+1)*(Nc+1)*Nr*(j-1)+(Nr+1)*(Nc+1)*(i-1);
                nh=2*Nc*(Nc*Nr*(j-1)+Nc*(i-1));
                for jh=1:Nr
                    for ih=1:Nc
                        if jh==1&&ih==1
                        n11=nh+(jh-1)*2*Nc+2*ih-1;
                        n12=nh+(jh-1)*2*Nc+2*ih;
                        
                        Tb(1,n11)=kh+(jh-1)*(Nc+1)+ih;
                        Tb(2,n11)=kh+(jh-1+1)*(Nc+1)+ih;
                        Tb(3,n11)=kh+(jh-1)*(Nc+1)+ih+1;
                                
                        Tb(1,n12)=kh+(jh-1)*(Nc+1)+ih+1;
                        Tb(2,n12)=kh+(jh-1+1)*(Nc+1)+ih;
                        Tb(3,n12)=kh+(jh-1+1)*(Nc+1)+ih+1;
                        end

                        if jh==1&&ih==2
                        n21=nh+(jh-1)*2*Nc+2*ih-1;
                        n22=nh+(jh-1)*2*Nc+2*ih;
                        
                        Tb(1,n21)=kh+(jh-1)*(Nc+1)+ih;
                        Tb(2,n21)=kh+(jh-1+1)*(Nc+1)+ih;
                        Tb(3,n21)=kh+(jh-1)*(Nc+1)+ih+1;
                                
                        

                        Tb(1,n22)=kh+(jh-1+1)*(Nc+1)+ih-1;
                        Tb(2,n22)=kh+(jh-1+1+1)*(Nc+1)+ih-1;
                        Tb(3,n22)=kh+(jh-1+1)*(Nc+1)+ih+1-1;
                        end

                        if jh==2&&ih==1
                        n31=nh+(jh-1)*2*Nc+2*ih-1;
                        n32=nh+(jh-1)*2*Nc+2*ih;
                        
                        Tb(1,n31)=kh+(jh-1-1)*(Nc+1)+ih+1+1;
                        Tb(2,n31)=kh+(jh-1+1-1)*(Nc+1)+ih+1;
                        Tb(3,n31)=kh+(jh-1+1-1)*(Nc+1)+ih+1+1;
                                
                        Tb(1,n32)=kh+(jh-1)*(Nc+1)+ih+1;
                        Tb(2,n32)=kh+(jh-1+1)*(Nc+1)+ih;
                        Tb(3,n32)=kh+(jh-1+1)*(Nc+1)+ih+1;
                        end

                        if jh==2&&ih==2
                        n41=nh+(jh-1)*2*Nc+2*ih-1;
                        n42=nh+(jh-1)*2*Nc+2*ih;
                        
                        Tb(1,n41)=kh+(jh-1)*(Nc+1)+ih;
                        Tb(2,n41)=kh+(jh-1+1)*(Nc+1)+ih;
                        Tb(3,n41)=kh+(jh-1)*(Nc+1)+ih+1;
                                
                        Tb(1,n42)=kh+(jh-1)*(Nc+1)+ih+1;
                        Tb(2,n42)=kh+(jh-1+1)*(Nc+1)+ih;
                        Tb(3,n42)=kh+(jh-1+1)*(Nc+1)+ih+1;
                        end
                    end
                end
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

