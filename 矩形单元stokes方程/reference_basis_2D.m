function result = reference_basis_2D(x,y,basis_type,basis_index,basis_der_x,basis_der_y)
%矩形剖分

if basis_type == 201 % 线性元
    switch basis_index
        case 1
            if basis_der_x==0&&basis_der_y==0
                result=(1-x-y+x.*y)/4;
            elseif basis_der_x==1&&basis_der_y==0
                result=(-1+y)/4;
            elseif basis_der_x==0 && basis_der_y==1
                result=(-1+x)/4;
            elseif basis_der_x==1&&basis_der_y==1
                result=1/4*ones(length(x),1);
            else
                result=0*ones(length(x),1);
            end
            
        case 2
            if basis_der_x==0&&basis_der_y==0
                result=(1-x+y-x.*y)/4;
            elseif basis_der_x==1&&basis_der_y==0
                result=(-1-y)/4;
            elseif basis_der_x==0&&basis_der_y==1
                result=(1-x)/4;
            elseif basis_der_x==1&&basis_der_y==1
                result=-1/4*ones(length(x),1);
            else 
                result=0*ones(length(x),1);
            end
        
        case 3
            if basis_der_x==0&&basis_der_y==0
                result=(1+x-y-x.*y)/4;
            elseif basis_der_x==1&&basis_der_y==0
                result=(1-y)/4;
            elseif basis_der_x==0&&basis_der_y==1
                result=(-1-x)/4;
            elseif basis_der_x==1&&basis_der_y==1
                result=-1/4*ones(length(x),1);
            else 
                result=0*ones(length(x),1);
            end
            
        case 4
            if basis_der_x==0&&basis_der_y==0
                result=(1+x+y+x.*y)/4;
            elseif basis_der_x==1&&basis_der_y==0
                result=(1+y)/4;
            elseif basis_der_x==0&&basis_der_y==1
                result=(1+x)/4;
            elseif basis_der_x==1&&basis_der_y==1
                result=1/4*ones(length(x),1);
            else 
                result=0*ones(length(x),1);
            end
    end

elseif basis_type == 202 % 二次元
    
    switch basis_index
        case 1  % 顶点基函数1
            if basis_der_x == 0 && basis_der_y == 0
                result=(1-x).*(1-y).*x.*y/4;  
            elseif basis_der_x == 1 && basis_der_y == 0
                result=(1-y).*y.*(1-2*x)/4;  
            elseif basis_der_x == 0 && basis_der_y == 1
                result=(1-x).*x.*(1-2*y)/4;  
            elseif basis_der_x == 2 && basis_der_y == 0
                result=-(1-y).*y/2;       
            elseif basis_der_x == 0 && basis_der_y == 2
                result=-(1-x).*x/2;         
            elseif basis_der_x == 1 && basis_der_y == 1
                result=(1-2*x).*(1-2*y)/4;   
            else
                result=0*ones(length(x),1);
            end
            
        case 2  % 边中点基函数2
            if basis_der_x == 0 && basis_der_y == 0
                result=-(1-x).*x.*(1-y.^2)/2;
            elseif basis_der_x == 1 && basis_der_y == 0
                result=-(1-y.^2).*(1-2*x)/2;  
            elseif basis_der_x == 0 && basis_der_y == 1
                result=(1-x).*x.*y;       
            elseif basis_der_x == 2 && basis_der_y == 0
                result=(1-y.^2);           
            elseif basis_der_x == 0 && basis_der_y == 2
                result=(1-x).*x;              
            elseif basis_der_x == 1 && basis_der_y == 1
                result=(1-2*x).*y;          
            else
                result=0*ones(length(x),1);
            end
            
        case 3  % 顶点基函数3
            if basis_der_x == 0 && basis_der_y == 0
                result=-(1-x).*(1+y).*x.*y/4;  
            elseif basis_der_x == 1 && basis_der_y == 0
                result=-(1+y).*y.*(1-2*x)/4;  
            elseif basis_der_x == 0 && basis_der_y == 1
                result=-(1-x).*x.*(1+2*y)/4;  
            elseif basis_der_x == 2 && basis_der_y == 0
                result=(1+y).*y/2;             
            elseif basis_der_x == 0 && basis_der_y == 2
                result=-(1-x).*x/2;            
            elseif basis_der_x == 1 && basis_der_y == 1
                result=-(1-2*x).*(1+2*y)/4;    
            else
                result=0*ones(length(x),1);
            end
            
        case 4  % 边中点基函数4
            if basis_der_x == 0 && basis_der_y == 0
                result=-(1-x.^2).*(1-y).*y/2;  
            elseif basis_der_x == 1 && basis_der_y == 0
                result=x.*(1-y).*y;            
            elseif basis_der_x == 0 && basis_der_y == 1
                result=-(1-x.^2).*(1-2*y)/2;   
            elseif basis_der_x == 2 && basis_der_y == 0
                result=(1-y).*y;              
            elseif basis_der_x == 0 && basis_der_y == 2
                result=(1-x.^2);            
            elseif basis_der_x == 1 && basis_der_y == 1
                result=x.*(1-2*y);       
            else
                result=0*ones(length(x),1);
            end
            
        case 5  % 中心基函数
            if basis_der_x == 0 && basis_der_y == 0
                result=(1-x.^2).*(1-y.^2);
            elseif basis_der_x == 1 && basis_der_y == 0
                result=2*x.*(y.^2-1);
            elseif basis_der_x == 0 && basis_der_y == 1
                result=2*y.*(x.^2-1);
            elseif basis_der_x == 2 && basis_der_y == 0
                result=2*(y.^2-1);
            elseif basis_der_x == 0 && basis_der_y == 2
                result=2*(x.^2-1);
            elseif basis_der_x == 1 && basis_der_y == 1
                result=4*x.*y;
            else
                result=0*ones(length(x),1);
            end
            
        case 6  % 边中点基函数6
            if basis_der_x == 0 && basis_der_y == 0
                result=(1-x.^2).*(1+y).*y/2;   
            elseif basis_der_x == 1 && basis_der_y == 0
                result=-x.*(1+y).*y;           
            elseif basis_der_x == 0 && basis_der_y == 1
                result=(1-x.^2).*(1+2*y)/2;    
            elseif basis_der_x == 2 && basis_der_y == 0
                result=-(1+y).*y;              
            elseif basis_der_x == 0 && basis_der_y == 2
                result=(1-x.^2);               
            elseif basis_der_x == 1 && basis_der_y == 1
                result=-x.*(1+2*y);            
            else
                result=0*ones(length(x),1);
            end
            
        case 7  % 顶点基函数7
            if basis_der_x == 0 && basis_der_y == 0
                result=-(1+x).*(1-y).*x.*y/4;  
            elseif basis_der_x == 1 && basis_der_y == 0
                result=-(1-y).*y.*(1+2*x)/4;   
            elseif basis_der_x == 0 && basis_der_y == 1
                result=-(1+x).*x.*(1-2*y)/4;   
            elseif basis_der_x == 2 && basis_der_y == 0
                result=-(1-y).*y/2;            
            elseif basis_der_x == 0 && basis_der_y == 2
                result=(1+x).*x/2;             
            elseif basis_der_x == 1 && basis_der_y == 1
                result=-(1+2*x).*(1-2*y)/4;    
            else
                result=0*ones(length(x),1);
            end
            
        case 8  % 边中点基函数8
            if basis_der_x == 0 && basis_der_y == 0
                result=(1+x).*x.*(1-y.^2)/2;   
            elseif basis_der_x == 1 && basis_der_y == 0
                result=(1-y.^2).*(1+2*x)/2;    
            elseif basis_der_x == 0 && basis_der_y == 1
                result=-(1+x).*x.*y;           
            elseif basis_der_x == 2 && basis_der_y == 0
                result=(1-y.^2);               
            elseif basis_der_x == 0 && basis_der_y == 2
                result=-(1+x).*x;              
            elseif basis_der_x == 1 && basis_der_y == 1
                result=-y.*(1+2*x);            
            else
                result=0*ones(length(x),1);
            end
            
        case 9  % 顶点基函数9
            if basis_der_x == 0 && basis_der_y == 0
                result=(1+x).*(1+y).*x.*y/4;   
            elseif basis_der_x == 1 && basis_der_y == 0
                result=(1+y).*y.*(1+2*x)/4;    
            elseif basis_der_x == 0 && basis_der_y == 1
                result=(1+x).*x.*(1+2*y)/4;    
            elseif basis_der_x == 2 && basis_der_y == 0
                result=(1+y).*y/2;             
            elseif basis_der_x == 0 && basis_der_y == 2
                result=(1+x).*x/2;             
            elseif basis_der_x == 1 && basis_der_y == 1
                result=(1+2*x).*(1+2*y)/4;     
            else
                result=0*ones(length(x),1);
            end
    end

end
end
