
num_of_element=64;%选择单元数（剖分份数）

gauss_type=3;%n点高斯公式

trial_basis_type=103;%试探函数u的基函数类型
test_basis_type=103;%测试函数v的基函数类型
%101表示一维线性元，102表示一维二次元，103表示一维三次元

%单元数4： maxerror_1=0.0023；   maxerror_2=4.6600e-5； maxerror_3=2.3825e-7
%单元数8： maxerror_1=5.8317e-4；maxerror_2=2.9919e-6； maxerror_3=1.0714e-8
%单元数16：maxerror_1=1.4645e-4；maxerror_2=1.8901e-7； maxerror_3=3.8460e-10
%单元数32：maxerror_1=3.6675e-5；maxerror_2=1.1869e-8； maxerror_3=1.2835e-11
%单元数64：maxerror_1=9.1700e-6；maxerror_2=7.4332e-10；maxerror_3=6.9573e-13
%误差阶分别收敛于2次（4倍）、         3次（8倍）、        4次（16倍）

%由于一维分片多项式形式简单，便没有仿射到参考单元而直接写出相应需要的局部基函数

%以下是一维满足相应弱形式方程问题的求解器
% -------------------------------------
%区分网格概念和有限元概念，在线性元下网格剖分跟有限元是相同的，故默认P，T是101的
basis_type=101;
[P,T] = generate_PbTb(num_of_element,basis_type);
%P是所有网格节点的坐标矩阵，T是所有单元上网格节点的全局编号矩阵

[P_trial,T_trial] = generate_PbTb(num_of_element, trial_basis_type);
[P_test,  T_test] = generate_PbTb(num_of_element, test_basis_type);
%两种函数根据有限元类型分别设置相应的Pb和Tb

matrix_size = [length(P_test),length(P_trial)];%A刚度矩阵
vector_size = length(P_trial);%b载荷向量

[gauss_weight,gauss_point] = gaussValues_1d(gauss_type); 
%根据确定的高斯积分公式提取相应的高斯节点和权重

trial_basis_der_A = 1;test_basis_der_A  = 1;
%矩阵A的测试函数试探函数的导数阶数

A = assemble_matrix_1d(matrix_size, gauss_point, gauss_weight, 'function_c', P, T, T_trial, T_test, trial_basis_type, trial_basis_der_A, test_basis_type, test_basis_der_A);    
%组装刚度矩阵A

test_basis_der_b = 0;
%向量b中测试函数的导数阶数

b = assemble_vector_1d(vector_size, gauss_point, gauss_weight, 'function_f', P, T, T_test, test_basis_type, test_basis_der_b);
%组装载荷向量b

[A,b] = boundary_process(A,b);
%处理狄利克雷边界条件（只有两个端点）
%处理诺尔曼和罗宾等混合边界条件，根据变分形式相应的添加到对应的分量

sol = A\b;

%L2误差和最大模误差
%error = norm(sol'-P_trial.*cos(P_trial))/sqrt(num_of_element);  
error = max(sol'-P_trial.*cos(P_trial));                         




