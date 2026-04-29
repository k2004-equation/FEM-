clear;
s=4;
num_of_element_r=s;%x轴行单元数
num_of_element_c=s;%y轴列单元数
h=1/s;
interval_r=[0,1];interval_c=[-0.25,0];%左右区间，上下区间
interval_t=[0,1];
nu=1;L=5;dt=1/64;
%theta=1;
gauss_type=9;%三角形单元n点高斯公式，可供使用3,4,9点（至少满足二次多项式精度要求）

trial_basis_type=202;%试探函数u的基函数类型
test_basis_type=202;%测试函数v的基函数类型
%202表示二维二次triangular,只能用这个
%试探函数p和测试函数q都是201，直接使用P、T

%单元数4×4,  二维maxerror_u=0.0054-----maxerror_p=1.8869
%单元数8×8,  二维maxerror_u=0.0047------maxerror_p=2.3626
%单元数16×16,二维maxerror_u=0.0046------maxerror_p=2.6706





%速度u使用二次元，压力p使用线性元

%以下是二维Navier-Stokes方程问题的求解器（采用三角形单元）
% -------------------------------------
%区分网格概念和有限元概念，在线性元下网格剖分跟有限元是相同的，故默认P，T是201的
basis_type=201;
[P,T] = triangular_generate_PbTb(num_of_element_r,num_of_element_c,interval_r,interval_c,basis_type);
%P是所有单元网格节点的坐标矩阵，T是所有单元上网格节点的全局编号矩阵

[P_trial,T_trial] = triangular_generate_PbTb(num_of_element_r,num_of_element_c,interval_r,interval_c, trial_basis_type);
[P_test,  T_test] = triangular_generate_PbTb(num_of_element_r,num_of_element_c,interval_r,interval_c,  test_basis_type);
%两种函数根据有限元类型分别设置相应的Pb和Tb

matrix_size_123 = [length(P_test),length(P_trial)];%A123刚度矩阵
matrix_size_56 = [length(P_test),length(P)];%A56刚度矩阵
vector_size = length(P_trial);%b12载荷向量

[boundary_nodes,boundary_edges]=boundary_nodes_edges(trial_basis_type,num_of_element_r,num_of_element_c);
%存储边界边和顶点的代数信息（顶点信息用于处理固定边界条件）

[gauss_weight,gauss_point]=generate_Gauss_reference_triangle(gauss_type);
%根据确定的高斯积分公式提取相应的高斯节点和权重（2D使用三角单元高斯节点和权重）


trial_basis_der_M_x = 0;trial_basis_der_M_y = 0;%矩阵M的试探函数的导数阶数
test_basis_der_M_x  = 0;test_basis_der_M_y  = 0;%矩阵M的测试函数的导数阶数
Me=assemble_matrix_2D(nu,matrix_size_123, gauss_point, gauss_weight, 'function_one', P, T, T_trial, T_test, trial_basis_type, trial_basis_der_M_x,trial_basis_der_M_y, test_basis_type, test_basis_der_M_x,test_basis_der_M_y);


trial_basis_der_A_x = 1;trial_basis_der_A_y = 0;%矩阵A中试探函数的导数阶数
test_basis_der_A_x  = 1;test_basis_der_A_y  = 0;%矩阵A中测试函数的导数阶数
A1 = assemble_matrix_2D(nu,matrix_size_123, gauss_point, gauss_weight, 'function_c1', P, T, T_trial, T_test, trial_basis_type, trial_basis_der_A_x,trial_basis_der_A_y, test_basis_type, test_basis_der_A_x,test_basis_der_A_y);    

trial_basis_der_A_x = 0;trial_basis_der_A_y = 1;
test_basis_der_A_x  = 0;test_basis_der_A_y  = 1;
A2 = assemble_matrix_2D(nu,matrix_size_123, gauss_point, gauss_weight, 'function_c1', P, T, T_trial, T_test, trial_basis_type, trial_basis_der_A_x,trial_basis_der_A_y, test_basis_type, test_basis_der_A_x,test_basis_der_A_y);    

trial_basis_der_A_x = 1;trial_basis_der_A_y = 0;
test_basis_der_A_x  = 0;test_basis_der_A_y  = 1;
A3 = assemble_matrix_2D(nu,matrix_size_123, gauss_point, gauss_weight, 'function_c1', P, T, T_trial, T_test, trial_basis_type, trial_basis_der_A_x,trial_basis_der_A_y, test_basis_type, test_basis_der_A_x,test_basis_der_A_y);    

basis_der_A_x = 0;basis_der_A_y = 0;
test_basis_der_A_x  = 1;test_basis_der_A_y  = 0;
A5 = assemble_matrix_2D(nu,matrix_size_56, gauss_point, gauss_weight, 'function_c2', P, T, T, T_test, basis_type, basis_der_A_x,basis_der_A_y, test_basis_type, test_basis_der_A_x,test_basis_der_A_y);    

basis_der_A_x = 0;basis_der_A_y = 0;
test_basis_der_A_x  = 0;test_basis_der_A_y  = 1;
A6 = assemble_matrix_2D(nu,matrix_size_56, gauss_point, gauss_weight, 'function_c2', P, T, T, T_test, basis_type, basis_der_A_x,basis_der_A_y, test_basis_type, test_basis_der_A_x,test_basis_der_A_y);    


O1=zeros(length(P),length(P));
O2=zeros(length(P_test),length(P));
O3=zeros(length(P_trial),length(P_trial));

M=[Me,O3,O2;
    O3',Me,O2;
    O2',O2',O1];%组装质量矩阵M

A=[2*A1+A2,A3,A5;
    A3',2*A2+A1,A6;
    A5',A6',O1];%组装刚度矩阵A


Xm_u1=inital_vector_u1(P_trial(1,:)',P_trial(2,:)');
Xm_u2=inital_vector_u2(P_trial(1,:)',P_trial(2,:)');
Xm_p=inital_vector_p(P(1,:)',P(2,:)');
Xm=[Xm_u1;Xm_u2;Xm_p];
%初值条件（当前时间解向量）

Nt=(interval_t(2)-interval_t(1))/dt;


%遍历时间迭代解向量，（单步法）且下一个时间解向量只需要当前时间解向量参与计算
for n=1:Nt

    tmp1=n*dt;
    
    test_basis_der_b_x = 0;test_basis_der_b_y=0;%向量b中测试函数的导数阶数
    b1 = assemble_vector_2D_t(nu,vector_size, gauss_point, gauss_weight, 'function_f1',tmp1, P, T, T_test, test_basis_type, test_basis_der_b_x,test_basis_der_b_y);

    test_basis_der_b_x = 0;test_basis_der_b_y=0;%向量b中测试函数的导数阶数
    b2 = assemble_vector_2D_t(nu,vector_size, gauss_point, gauss_weight, 'function_f2',tmp1, P, T, T_test, test_basis_type, test_basis_der_b_x,test_basis_der_b_y);
    
    O4=zeros(length(P),1);
    
    bmp1=[b1;b2;O4];%组装载荷向量b
    XL=Xm;
    for i=1:L

        ul=XL;
        

        trial_basis_type_c=202;
    
        trial_basis_der_A_x_c = 1;trial_basis_der_A_y_c = 0;
        trial_basis_der_A_x = 0;trial_basis_der_A_y = 0;%矩阵A中试探函数的导数阶数
        test_basis_der_A_x  = 0;test_basis_der_A_y  = 0;%矩阵A中测试函数的导数阶数
        AN1 = assemble_matrix_2D_ns(matrix_size_123, gauss_point, gauss_weight, 'function_c1h',ul,trial_basis_type_c, trial_basis_der_A_x_c ,trial_basis_der_A_y_c , P, T,P_trial, T_trial, T_test, trial_basis_type, trial_basis_der_A_x,trial_basis_der_A_y, test_basis_type, test_basis_der_A_x,test_basis_der_A_y);    
        
        trial_basis_der_A_x_c = 0;trial_basis_der_A_y_c = 0;
        trial_basis_der_A_x = 1;trial_basis_der_A_y = 0;
        test_basis_der_A_x  = 0;test_basis_der_A_y  = 0;
        AN2 = assemble_matrix_2D_ns(matrix_size_123, gauss_point, gauss_weight, 'function_c1h', ul,trial_basis_type_c, trial_basis_der_A_x_c ,trial_basis_der_A_y_c ,P, T,P_trial, T_trial, T_test, trial_basis_type, trial_basis_der_A_x,trial_basis_der_A_y, test_basis_type, test_basis_der_A_x,test_basis_der_A_y);    
        
        trial_basis_der_A_x_c = 0;trial_basis_der_A_y_c = 0;
        trial_basis_der_A_x = 0;trial_basis_der_A_y = 1;
        test_basis_der_A_x  = 0;test_basis_der_A_y  = 0;
        AN3 = assemble_matrix_2D_ns(matrix_size_123, gauss_point, gauss_weight, 'function_c2h', ul,trial_basis_type_c, trial_basis_der_A_x_c ,trial_basis_der_A_y_c ,P, T,P_trial, T_trial, T_test, trial_basis_type, trial_basis_der_A_x,trial_basis_der_A_y, test_basis_type, test_basis_der_A_x,test_basis_der_A_y);    
        
        trial_basis_der_A_x_c = 0;trial_basis_der_A_y_c = 1;
        trial_basis_der_A_x = 0;trial_basis_der_A_y = 0;
        test_basis_der_A_x  = 0;test_basis_der_A_y  = 0;
        AN4 = assemble_matrix_2D_ns(matrix_size_123, gauss_point, gauss_weight, 'function_c1h', ul,trial_basis_type_c, trial_basis_der_A_x_c ,trial_basis_der_A_y_c ,P, T,P_trial, T_trial, T_test, trial_basis_type, trial_basis_der_A_x,trial_basis_der_A_y, test_basis_type, test_basis_der_A_x,test_basis_der_A_y);    
        
        trial_basis_der_A_x_c = 1;trial_basis_der_A_y_c = 0;
        trial_basis_der_A_x = 0;trial_basis_der_A_y = 0;
        test_basis_der_A_x  = 0;test_basis_der_A_y  = 0;
        AN5 = assemble_matrix_2D_ns(matrix_size_123, gauss_point, gauss_weight, 'function_c2h', ul,trial_basis_type_c, trial_basis_der_A_x_c ,trial_basis_der_A_y_c ,P, T,P_trial, T_trial, T_test, trial_basis_type, trial_basis_der_A_x,trial_basis_der_A_y, test_basis_type, test_basis_der_A_x,test_basis_der_A_y);    
        
        trial_basis_der_A_x_c = 0;trial_basis_der_A_y_c = 1;
        trial_basis_der_A_x = 0;trial_basis_der_A_y = 0;
        test_basis_der_A_x  = 0;test_basis_der_A_y  = 0;
        AN6 = assemble_matrix_2D_ns(matrix_size_123, gauss_point, gauss_weight, 'function_c2h', ul,trial_basis_type_c, trial_basis_der_A_x_c ,trial_basis_der_A_y_c ,P, T, P_trial,T_trial, T_test, trial_basis_type, trial_basis_der_A_x,trial_basis_der_A_y, test_basis_type, test_basis_der_A_x,test_basis_der_A_y);    
        
        
        AN=[AN1+AN2+AN3,AN4,O2;
            AN5,AN6+AN2+AN3,O2;
            O2',O2',O1];
    
    
        trial_basis_type_c1=202;trial_basis_type_c2=202;
    
        trial_type_1=1;trial_type_2=1;
        trial_basis_der_b_x_c1 = 0;trial_basis_der_b_y_c1 = 0;
        trial_basis_der_b_x_c2 = 1;trial_basis_der_b_y_c2 = 0;
        test_basis_der_b_x = 0;test_basis_der_b_y=0;%向量b中测试函数的导数阶数
        bN1 = assemble_vector_2D_ns(vector_size, gauss_point, gauss_weight, 'function_uh', ul,trial_type_1,trial_type_2,trial_basis_type_c1, trial_basis_der_b_x_c1 ,trial_basis_der_b_y_c1 ,trial_basis_type_c2, trial_basis_der_b_x_c2 ,trial_basis_der_b_y_c2 ,P, T,P_trial,T_trial, T_test, test_basis_type, test_basis_der_b_x,test_basis_der_b_y);
        
        trial_type_1=2;trial_type_2=1;
        trial_basis_der_b_x_c1 = 0;trial_basis_der_b_y_c1 = 0;
        trial_basis_der_b_x_c2 = 0;trial_basis_der_b_y_c2 = 1;
        test_basis_der_b_x = 0;test_basis_der_b_y=0;%向量b中测试函数的导数阶数
        bN2 = assemble_vector_2D_ns(vector_size, gauss_point, gauss_weight, 'function_uh', ul,trial_type_1,trial_type_2,trial_basis_type_c1, trial_basis_der_b_x_c1 ,trial_basis_der_b_y_c1 ,trial_basis_type_c2, trial_basis_der_b_x_c2 ,trial_basis_der_b_y_c2 ,P, T,P_trial,T_trial, T_test, test_basis_type, test_basis_der_b_x,test_basis_der_b_y);
        
        trial_type_1=1;trial_type_2=2;
        trial_basis_der_b_x_c1 = 0;trial_basis_der_b_y_c1 = 0;
        trial_basis_der_b_x_c2 = 1;trial_basis_der_b_y_c2 = 0;
        test_basis_der_b_x = 0;test_basis_der_b_y=0;%向量b中测试函数的导数阶数
        bN3 = assemble_vector_2D_ns(vector_size, gauss_point, gauss_weight, 'function_uh', ul,trial_type_1,trial_type_2,trial_basis_type_c1, trial_basis_der_b_x_c1 ,trial_basis_der_b_y_c1 ,trial_basis_type_c2, trial_basis_der_b_x_c2 ,trial_basis_der_b_y_c2 ,P, T, P_trial,T_trial,T_test, test_basis_type, test_basis_der_b_x,test_basis_der_b_y);
        
        trial_type_1=2;trial_type_2=2;
        trial_basis_der_b_x_c1 = 0;trial_basis_der_b_y_c1 = 0;
        trial_basis_der_b_x_c2 = 0;trial_basis_der_b_y_c2 = 1;
        test_basis_der_b_x = 0;test_basis_der_b_y=0;%向量b中测试函数的导数阶数
        bN4 = assemble_vector_2D_ns(vector_size, gauss_point, gauss_weight, 'function_uh', ul,trial_type_1,trial_type_2,trial_basis_type_c1, trial_basis_der_b_x_c1 ,trial_basis_der_b_y_c1 ,trial_basis_type_c2, trial_basis_der_b_x_c2 ,trial_basis_der_b_y_c2 ,P, T,P_trial,T_trial, T_test, test_basis_type, test_basis_der_b_x,test_basis_der_b_y);
        
        bN=[bN1+bN2;bN3+bN4;O4];

        AL=M/dt+A+AN;
        bL=bmp1+M/dt*Xm+bN;
        
        [AL,bL] = dirichlet_boundary_process('function_g1','function_g2','function_p',AL,bL,tmp1,boundary_nodes,P_trial,P);
        
        XL=AL\bL;%下一时间的解向量

        
    end

    Xm=XL;
end


ul=Xm;

u=[function_u1(P_trial(1,:)',P_trial(2,:)',tmp1);function_u2(P_trial(1,:)',P_trial(2,:)',tmp1);function_p(P(1,:)',P(2,:)',tmp1)];
%离散解析解

%最大模误差

error_u = max(abs(full(ul(1:2*length(P_trial)))-u(1:2*length(P_trial))));   
error_p = max(abs(full(ul(2*length(P_trial)+1:end))-u(2*length(P_trial)+1:end)));   

u_comparison_u=[u,full(ul)];


