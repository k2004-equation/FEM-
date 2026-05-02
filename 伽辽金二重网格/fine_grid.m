s=s^2;r=202;tic
num_of_element_r=s;%x轴行单元数
num_of_element_c=s;%y轴列单元数
h=1/s;
interval_r=[0,1];interval_c=[0,1];%左右区间，上下区间
interval_t=[0,1];
dt=h^2;
%theta=1;
gauss_type=9;%三角形单元n点高斯公式，可供使用3,4,9点（至少满足二次多项式精度要求）

trial_basis_type=r;%试探函数u的基函数类型
test_basis_type=r;%测试函数v的基函数类型

%最后一层时间的解误差估计，粗网格s×s，细网格s^2×s^2
%------------无穷范数------L2范数------H1范数-------时间成本
%2×2 4×4    
%4×4 16×16
%8×8 64×64




%以下是二维非稳态非线性抛物线方程问题的求解器（采用三角形单元）
%fine grid细网格
% -------------------------------------
%区分网格概念和有限元概念，在线性元下网格剖分跟有限元是相同的，故默认P，T是201的
basis_type=201;
[P,T] = triangular_generate_PbTb(num_of_element_r,num_of_element_c,interval_r,interval_c,basis_type);
%P是所有单元网格节点的坐标矩阵，T是所有单元上网格节点的全局编号矩阵

[P_trial,T_trial] = triangular_generate_PbTb(num_of_element_r,num_of_element_c,interval_r,interval_c, trial_basis_type);
[P_test,  T_test] = triangular_generate_PbTb(num_of_element_r,num_of_element_c,interval_r,interval_c,  test_basis_type);
%两种函数根据有限元类型分别设置相应的Pb和Tb

matrix_size_A = [length(P_test),length(P_trial)];

vector_size_b = length(P_test);

[boundary_nodes,boundary_edges]=boundary_nodes_edges(trial_basis_type,num_of_element_r,num_of_element_c);
%存储边界边和顶点的代数信息（顶点信息用于处理固定边界条件）

[gauss_weight,gauss_point]=generate_Gauss_reference_triangle(gauss_type);
%根据确定的高斯积分公式提取相应的高斯节点和权重（2D使用三角单元高斯节点和权重）


Ym(:,1)=inital_vector_u(P_trial(1,:)',P_trial(2,:)');
%初值条件（当前时间解向量）

Nt=(interval_t(2)-interval_t(1))/dt;

Xm=uh;
%遍历时间迭代解向量，（单步法）且下一个时间解向量只需要当前时间解向量参与计算
for n=1:Nt

    tn=n*dt;
    
    test_basis_der_b_x = 0;test_basis_der_b_y=0;%向量b中测试函数的导数阶数
    bhn = assemble_vector_2D_t(vector_size_b, gauss_point, gauss_weight, 'function_fn',tn, P, T, T_test, test_basis_type, test_basis_der_b_x,test_basis_der_b_y);
    

    trial_basis_type_c=r;
    trial_basis_der_x_c=0;trial_basis_der_y_c=0;
    trial_basis_der_A_x = 0;trial_basis_der_A_y = 0;%矩阵A中试探函数的导数阶数
    test_basis_der_A_x  = 0;test_basis_der_A_y  = 0;%矩阵A中测试函数的导数阶数
    Ah1 = assemble_matrix_2D_ns(matrix_size_A, gauss_point, gauss_weight, 'function_feval_c1', 'function_gu',Xm(:,n+1),P, T, T_trial, T_test, trial_basis_type_c, trial_basis_der_x_c,trial_basis_der_y_c,trial_basis_type, trial_basis_der_A_x,trial_basis_der_A_y, test_basis_type, test_basis_der_A_x,test_basis_der_A_y);    
    
    trial_basis_der_x_c=0;trial_basis_der_y_c=0;
    trial_basis_der_A_x = 0;trial_basis_der_A_y = 0;
    test_basis_der_A_x  = 0;test_basis_der_A_y  = 0;
    Ah2 = assemble_matrix_2D_ns_h(matrix_size_A, gauss_point, gauss_weight, 'function_feval_c2_h', 'function_gudu',Xm(:,n+1),Xm(:,n),P, T, T_trial, T_test, trial_basis_type_c, trial_basis_der_x_c,trial_basis_der_y_c,trial_basis_type, trial_basis_der_A_x,trial_basis_der_A_y, test_basis_type, test_basis_der_A_x,test_basis_der_A_y);    
    
    
    trial_basis_der_A_x = 1;trial_basis_der_A_y = 0;
    test_basis_der_A_x  = 1;test_basis_der_A_y  = 0;
    Ah3 = assemble_matrix_2D(matrix_size_A, gauss_point, gauss_weight, 'function_c3', P, T, T_trial, T_test,trial_basis_type, trial_basis_der_A_x,trial_basis_der_A_y, test_basis_type, test_basis_der_A_x,test_basis_der_A_y);    
    
    
    trial_basis_der_A_x = 0;trial_basis_der_A_y = 1;
    test_basis_der_A_x  = 0;test_basis_der_A_y  = 1;
    Ah4 = assemble_matrix_2D(matrix_size_A, gauss_point, gauss_weight, 'function_c4', P, T, T_trial, T_test, trial_basis_type, trial_basis_der_A_x,trial_basis_der_A_y, test_basis_type, test_basis_der_A_x,test_basis_der_A_y);    
    
      
    Ah=Ah1/dt+Ah2/dt+Ah3+Ah4;
    

    trial_basis_der_b_x = 0;trial_basis_der_b_y = 0;
    test_basis_der_b_x = 0;test_basis_der_b_y=0;%向量b中测试函数的导数阶数
    bh1 = assemble_vector_2D_ns_ns(vector_size_b, gauss_point, gauss_weight, 'function_feval_c1','function_gu', Xm(:,n+1),Ym(:,n),trial_basis_type, trial_basis_der_b_x,trial_basis_der_b_y,P, T,P_trial,T_trial, T_test, test_basis_type, test_basis_der_b_x,test_basis_der_b_y);
    
    trial_basis_der_b_x = 0;trial_basis_der_b_y = 0;
    test_basis_der_b_x = 0;test_basis_der_b_y = 0;%向量b中测试函数的导数阶数
    bh2 = assemble_vector_2D_ns_h(vector_size_b, gauss_point, gauss_weight, 'function_feval_c3_h','function_gudu', Xm(:,n+1),Xm(:,n),trial_basis_type, trial_basis_der_b_x,trial_basis_der_b_y,P, T,P_trial,T_trial, T_test, test_basis_type, test_basis_der_b_x,test_basis_der_b_y);
    
      
    bh=bhn+bh1/dt+bh2/dt;

    [Ah,bh] = dirichlet_boundary_process('function_dirichlet',Ah,bh,tn,boundary_nodes,P_trial);
    
    Yh=Ah\bh;%下一时间的解向量
    
        
    

    Ym(:,n+1)=Yh;
end



for n=1:Nt+1

uhh(:,n)=function_u(P_trial(1,:)',P_trial(2,:)',(n-1)*dt);

end

%离散解析解


u_comparison_uh=[uhh,full(Ym)];


%无穷范数误差估计
error_max_u = max(abs(full(Ym)-uhh));   

%L2范数误差估计
error_L2_u=zeros(1,Nt+1);
trial_basis_der_A_x=0;trial_basis_der_A_y=0;
for n=1:Nt+1
error_L2_u(n)=error_L2_assemble_2D(uhh(:,n),Ym(:,n),gauss_point, gauss_weight,P, T, T_trial, T_test, trial_basis_type, trial_basis_der_A_x,trial_basis_der_A_y);
error_L2_u(n)=sqrt(error_L2_u(n));
end

%H1半范数误差估计
error_H1_u=zeros(1,Nt+1);
trial_basis_der_A_x1=1;trial_basis_der_A_y1=0;
trial_basis_der_A_x2=0;trial_basis_der_A_y2=1;
for n=1:Nt+1
error_H1_u(n)=error_L2_assemble_2D(uhh(:,n),Ym(:,n),gauss_point, gauss_weight,P, T, T_trial, T_test, trial_basis_type, trial_basis_der_A_x1,trial_basis_der_A_y1)+error_L2_assemble_2D(uhh(:,n),Ym(:,n),gauss_point, gauss_weight,P, T, T_trial, T_test, trial_basis_type, trial_basis_der_A_x2,trial_basis_der_A_y2);

end

error_H1_uh=sqrt(sum(error_H1_u)*dt);
error=[error_max_u(Nt+1);error_L2_u(Nt+1);error_H1_uh]

toc