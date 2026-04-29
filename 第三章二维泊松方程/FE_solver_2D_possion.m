
num_of_element_r=2;%x轴行单元数
num_of_element_c=2;%y轴列单元数

interval_r=[-1,1];interval_c=[-1,1];%左右区间，上下区间

gauss_type=9;%三角形单元n点高斯公式，可供使用3,4,9点（至少满足二次多项式精度要求）

trial_basis_type=202;%试探函数u的基函数类型
test_basis_type=202;%测试函数v的基函数类型
%201表示二维线性triangular，202表示二维二次triangular

%单元数2×2,二维线性maxerror=0.0683
%单元数4×4,二维线性maxerror=0.0249
%单元数8×8,二维线性maxerror=0.0067
%单元数16×16,二维线性maxerror=0.0017
%单元数32×32,二维线性maxerror=4.3521e-4
%误差阶收敛于2次，四倍
%单元数2×2,二维二次maxerror=0.0161
%单元数4×4,二维二次maxerror=0.0019
%单元数8×8,二维二次maxerror=2.0116e-4
%单元数16×16,二维二次maxerror=1.7836e-5
%单元数32×32,二维二次maxerror=1.3528e-6
%误差阶收敛于3次，八倍


%以下是二维满足相应弱形式方程问题的求解器（采用三角形单元）
% -------------------------------------
%区分网格概念和有限元概念，在线性元下网格剖分跟有限元是相同的，故默认P，T是201的
basis_type=201;
[P,T] = triangular_generate_PbTb(num_of_element_r,num_of_element_c,interval_r,interval_c,basis_type);
%P是所有单元网格节点的坐标矩阵，T是所有单元上网格节点的全局编号矩阵

[P_trial,T_trial] = triangular_generate_PbTb(num_of_element_r,num_of_element_c,interval_r,interval_c, trial_basis_type);
[P_test,  T_test] = triangular_generate_PbTb(num_of_element_r,num_of_element_c,interval_r,interval_c,  test_basis_type);
%两种函数根据有限元类型分别设置相应的Pb和Tb

matrix_size = [length(P_test),length(P_trial)];%A刚度矩阵（后期不仅限于刚度矩阵）
vector_size = length(P_trial);%b载荷向量

[boundary_nodes,boundary_edges]=boundary_nodes_edges(trial_basis_type,num_of_element_r,num_of_element_c);
%存储边界边和顶点的代数信息（顶点信息用于处理固定边界条件）

[gauss_weight,gauss_point]=generate_Gauss_reference_triangle(gauss_type);
%根据确定的高斯积分公式提取相应的高斯节点和权重（2D使用三角单元高斯节点和权重）

trial_basis_der_A_x = 1;trial_basis_der_A_y = 0;%矩阵A中试探函数的导数阶数
test_basis_der_A_x  = 1;test_basis_der_A_y  = 0;%矩阵A中测试函数的导数阶数
A1 = assemble_matrix_2D(matrix_size, gauss_point, gauss_weight, 'function_c', P, T, T_trial, T_test, trial_basis_type, trial_basis_der_A_x,trial_basis_der_A_y, test_basis_type, test_basis_der_A_x,test_basis_der_A_y);    


trial_basis_der_A_x = 0;trial_basis_der_A_y = 1;
test_basis_der_A_x  = 0;test_basis_der_A_y  = 1;
A2 = assemble_matrix_2D(matrix_size, gauss_point, gauss_weight, 'function_c', P, T, T_trial, T_test, trial_basis_type, trial_basis_der_A_x,trial_basis_der_A_y, test_basis_type, test_basis_der_A_x,test_basis_der_A_y);    

A=A1+A2;%组装刚度矩阵A，核心算法

test_basis_der_b_x = 0;test_basis_der_b_y=0;%向量b中测试函数的导数阶数
b = assemble_vector_2D(vector_size, gauss_point, gauss_weight, 'function_f', P, T, T_test, test_basis_type, test_basis_der_b_x,test_basis_der_b_y);
%组装载荷向量b，核心算法


[A,b] = dirichlet_boundary_process('function_g',A,b,boundary_nodes,P_trial);
%处理狄利克雷边界条件
%处理诺伊曼和罗宾等混合边界条件，根据变分形式相应的添加到对应的分量（本算例没有）

sol = A\b;%数值解


u=function_u(P_trial(1,:),P_trial(2,:));%离散解析解


%最大模误差
error = max(abs(full(sol)'-u));    

u_comparison=[u',full(sol)];





