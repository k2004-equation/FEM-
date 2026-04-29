s=4;
num_of_element_r=s;%x轴行单元数
num_of_element_c=s;%y轴列单元数

interval_r=[0,1];interval_c=[-0.25,0];%左右区间，上下区间

nu=1;

gauss_type=9;%三角形单元n点高斯公式，可供使用3,4,9点（至少满足二次多项式精度要求）

trial_basis_type=202;%试探函数u的基函数类型
test_basis_type=202;%测试函数v的基函数类型
%202表示二维二次triangular,只能用这个
%试探函数p和测试函数q都是201，直接使用P、T

%单元数2×2,  二维线性maxerror_u=0.0046--------maxerror_p=1.1051
%单元数4×4,  二维线性maxerror_u=0.0011--------maxerror_p=0.2040
%单元数8×8,  二维线性maxerror_u=9.3090e-5-----maxerror_p=0.0436
%单元数16×16,二维线性maxerror_u=6.5267e-6-----maxerror_p=0.0103
%单元数32×32,二维线性maxerror_u=4.2574e-7-----maxerror_p=0.0026
%单元数64×64,二维线性maxerror_u=2.6614e-8-----maxerror_p=6.4363e-4
%误差阶--------------4次（十六倍）---------------2次（四倍）
%理论上速度u的误差阶应该是3次，此处出现超收敛现象？

%速度u使用二次元，压力p使用线性元

%以下是二维Stokes向量方程问题的求解器（采用三角形单元）
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

A=[2*A1+A2,A3,A5;
    A3',2*A2+A1,A6;
    A5',A6',O1];%组装刚度矩阵A



test_basis_der_b_x = 0;test_basis_der_b_y=0;%向量b中测试函数的导数阶数
b1 = assemble_vector_2D(nu,vector_size, gauss_point, gauss_weight, 'function_f1', P, T, T_test, test_basis_type, test_basis_der_b_x,test_basis_der_b_y);

test_basis_der_b_x = 0;test_basis_der_b_y=0;%向量b中测试函数的导数阶数
b2 = assemble_vector_2D(nu,vector_size, gauss_point, gauss_weight, 'function_f2', P, T, T_test, test_basis_type, test_basis_der_b_x,test_basis_der_b_y);

O2=zeros(length(P),1);

b=[b1;b2;O2];%组装载荷向量b

[A,b] = dirichlet_boundary_process('function_g1','function_g2','function_p',A,b,boundary_nodes,P_trial,P);
%处理狄利克雷边界条件
%处理诺尔曼和罗宾等混合边界条件，根据变分形式相应的添加到对应的分量

sol = A\b;

u=[function_u1(P_trial(1,:)',P_trial(2,:)');function_u2(P_trial(1,:)',P_trial(2,:)');function_p(P(1,:)',P(2,:)')];
%离散解析解

%最大模误差

error_u = max(abs(full(sol(1:2*length(P_trial)))-u(1:2*length(P_trial))));   
error_p = max(abs(full(sol(2*length(P_trial)+1:end))-u(2*length(P_trial)+1:end)));   

u_comparison_u=[u,full(sol)];


