
num_of_element_r=2;%x轴行单元数
num_of_element_c=2;%y轴列单元数

interval_r=[0,2];interval_c=[0,1];%左右区间，上下区间
interval_t=[0,1];%时间长度

gauss_type=9;%三角形单元n点高斯公式，可供使用3,4,9点（至少满足二次多项式精度要求）

h=1/num_of_element_c;

%别忘了改基函数类型
theta=1;dt=4*h^2;%三角形单元线性有限元，向后欧拉方案
%theta=1/2;dt=h;%三角形单元线性有限元，C-N中心方案
%theta=1;dt=8*h^3;%三角形单元二次有限元，向后欧拉方案
%theta=1/2;dt=sqrt(h^3);%三角形单元二次有限元，C-N中心方案

trial_basis_type=201;%试探函数u的基函数类型
test_basis_type=201;%测试函数v的基函数类型
%201表示二维三角形单元线性元，202表示二维三角形单元二次元

%单元数---线性FE向后FD--线性FE中心FD--二次FE向后FD--二次FE中心FD--maxerror
%2×2------0.1818-------0.0394--------0.2930-------0.0093
%4×4------0.0757-------0.0162--------0.0551-------0.0021
%8×8------0.0212-------0.0046--------0.0070------2.0808e-4
%16×16----0.0055-------0.0012-------8.8622e-4----2.3133e-5
%32×32----0.0014------2.9333e-4-----1.1048e-4----2.6825e-6
%误差阶----2次（4倍）---2次（4倍）----3次（8倍）----3次（8倍）

%二次FE向后FD的计算成本15分钟，但同样的有限元，CN格式的计算成本很小。
%A中系数项c与时间t无关
%以下是二维满足相应弱形式方程非稳态问题的求解器（热传导方程）
% -------------------------------------
%区分网格概念和有限元概念，在线性元下网格剖分跟有限元是相同的，故默认P，T是101的
[P,T] = triangular_generate_PbTb(num_of_element_r,num_of_element_c,interval_r,interval_c,201);
%P是所有网格节点的坐标矩阵，T是所有单元上网格节点的全局编号矩阵

[P_trial,T_trial] = triangular_generate_PbTb(num_of_element_r,num_of_element_c,interval_r,interval_c, trial_basis_type);
[P_test,  T_test] = triangular_generate_PbTb(num_of_element_r,num_of_element_c,interval_r,interval_c,  test_basis_type);
%两种函数根据有限元类型分别设置相应的Pb和Tb

matrix_size = [length(P_test),length(P_trial)];%A刚度矩阵
vector_size = length(P_trial);%b载荷向量


[gauss_weight,gauss_point]=generate_Gauss_reference_triangle(gauss_type);
%根据确定的高斯积分公式提取相应的高斯节点和权重

[boundary_nodes,boundary_edges]=boundary_nodes_edges(trial_basis_type,num_of_element_r,num_of_element_c);
%存储边界边和顶点的代数信息（顶点信息用于处理固定边界条件）

trial_basis_der_M_x = 0;trial_basis_der_M_y = 0;%矩阵M的试探函数的导数阶数
test_basis_der_M_x  = 0;test_basis_der_M_y  = 0;%矩阵M的测试函数的导数阶数
M=assemble_matrix_2D(matrix_size, gauss_point, gauss_weight, 'function_one', P, T, T_trial, T_test, trial_basis_type, trial_basis_der_M_x,trial_basis_der_M_y, test_basis_type, test_basis_der_M_x,test_basis_der_M_y);
%组装质量矩阵M

trial_basis_der_A_x = 1;trial_basis_der_A_y = 0;%矩阵A的试探函数的导数阶数
test_basis_der_A_x  = 1;test_basis_der_A_y  = 0;%矩阵A的测试函数的导数阶数
A1 = assemble_matrix_2D(matrix_size, gauss_point, gauss_weight, 'function_c', P, T, T_trial, T_test, trial_basis_type, trial_basis_der_A_x,trial_basis_der_A_y, test_basis_type, test_basis_der_A_x,test_basis_der_A_y);    


trial_basis_der_A_x = 0;trial_basis_der_A_y = 1;
test_basis_der_A_x  = 0;test_basis_der_A_y  = 1;
A2 = assemble_matrix_2D(matrix_size, gauss_point, gauss_weight, 'function_c', P, T, T_trial, T_test, trial_basis_type, trial_basis_der_A_x,trial_basis_der_A_y, test_basis_type, test_basis_der_A_x,test_basis_der_A_y);    


A=A1+A2;%组装刚度矩阵A

Atilde=M/dt+theta*A;
Afixed=M/dt-(1-theta)*A;

Xm=inital_vector(P_trial(1,:)',P_trial(2,:)');
%初值条件（当前时间解向量）

Nt=(interval_t(2)-interval_t(1))/dt;

test_basis_der_b_x = 0;test_basis_der_b_y=0;%向量b中测试函数的导数阶数

%遍历时间迭代解向量，（单步法）且下一个时间解向量只需要当前时间解向量参与计算
for i=1:Nt
    tm=(i-1)*dt;
    tmp1=i*dt;
    
    bm = assemble_vector_2D(vector_size, gauss_point, gauss_weight, 'function_f',tm, P, T, T_test, test_basis_type, test_basis_der_b_x,test_basis_der_b_y);
    bmp1 = assemble_vector_2D(vector_size, gauss_point, gauss_weight, 'function_f',tmp1, P, T, T_test, test_basis_type, test_basis_der_b_x,test_basis_der_b_y);
    
    btilde=theta*bmp1+(1-theta)*bm+Afixed*Xm;
    
    [Atilde,btilde] = dirichlet_boundary_process('function_g',Atilde,btilde ,tmp1,boundary_nodes,P_trial);
    
    Xmp1=Atilde\btilde;%下一时间的解向量

    Xm=Xmp1;

end

u=function_u(P_trial(1,:),P_trial(2,:),tmp1);
%离散解析解

%最大模误差
error = max(abs(full(Xm)-u')); 

u_comparison=[u',full(Xm)];





