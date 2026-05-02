clear;tic
s=4;r=201;
num_of_element_r=s;%x轴行单元数
num_of_element_c=s;%y轴列单元数
h=1/s;
interval_r=[0,1];interval_c=[0,1];%左右区间，上下区间
interval_t=[0,1];
dt=h^(r-200+1);
%theta=1;
gauss_type=9;%三角形单元n点高斯公式，可供使用3,4,9点（至少满足二次多项式精度要求）

trial_basis_type=r;%试探函数u的基函数类型
test_basis_type=r;%测试函数v的基函数类型


%以下是二维非稳态非线性抛物线方程问题的求解器（采用三角形单元）
% -------------------------------------
%区分网格概念和有限元概念，在线性元下网格剖分跟有限元是相同的，故默认P，T是201的
basis_type=201;
[P,T] = triangular_generate_PbTb(num_of_element_r,num_of_element_c,interval_r,interval_c,basis_type);
%P是所有单元网格节点的坐标矩阵，T是所有单元上网格节点的全局编号矩阵

[P_trial,T_trial] = triangular_generate_PbTb(num_of_element_r,num_of_element_c,interval_r,interval_c, trial_basis_type);
[P_test,  T_test] = triangular_generate_PbTb(num_of_element_r,num_of_element_c,interval_r,interval_c,  test_basis_type);
%两种函数根据有限元类型分别设置相应的Pb和Tb（支持201和202）

matrix_size_A = [length(P_test),length(P_trial)];

vector_size_b = length(P_test);

[boundary_nodes,boundary_edges]=boundary_nodes_edges(trial_basis_type,num_of_element_r,num_of_element_c);
%存储边界边和顶点的代数信息（顶点信息用于处理固定边界条件）

[gauss_weight,gauss_point]=generate_Gauss_reference_triangle(gauss_type);
%根据确定的高斯积分公式提取相应的高斯节点和权重（2D使用三角单元高斯节点和权重）


Xm=inital_vector_u(P_trial(1,:)',P_trial(2,:)');
%初值条件（当前时间解向量）

Nt=(interval_t(2)-interval_t(1))/dt;


%时间迭代
for n=1:Nt

    tn=n*dt;
    
    test_basis_der_b_x = 0;test_basis_der_b_y=0;
    bn = assemble_vector_2D_t(vector_size_b, gauss_point, gauss_weight, 'function_fn',tn, P, T, T_test, test_basis_type, test_basis_der_b_x,test_basis_der_b_y);

    XL=Xm(:,n);
    dXL=ones(length(P_trial),1);

    while max(abs(dXL))>1e-9

        trial_basis_type_c = r;%牛顿迭代被积函数系数的基函数类型

        trial_basis_der_x_c = 0;trial_basis_der_y_c = 0;
        trial_basis_der_A_x = 0;trial_basis_der_A_y = 0;
        test_basis_der_A_x  = 0;test_basis_der_A_y  = 0;
        A1 = assemble_matrix_2D_ns(matrix_size_A, gauss_point, gauss_weight, 'function_feval_c1', 'function_gu',XL,P, T, T_trial, T_test, trial_basis_type_c, trial_basis_der_x_c,trial_basis_der_y_c,trial_basis_type, trial_basis_der_A_x,trial_basis_der_A_y, test_basis_type, test_basis_der_A_x,test_basis_der_A_y);    
        
        trial_basis_der_x_c = 0;trial_basis_der_y_c = 0;
        trial_basis_der_A_x = 0;trial_basis_der_A_y = 0;
        test_basis_der_A_x  = 0;test_basis_der_A_y  = 0;
        A2 = assemble_matrix_2D_ns(matrix_size_A, gauss_point, gauss_weight, 'function_feval_c2', 'function_gudu',XL,P, T, T_trial, T_test, trial_basis_type_c, trial_basis_der_x_c,trial_basis_der_y_c,trial_basis_type, trial_basis_der_A_x,trial_basis_der_A_y, test_basis_type, test_basis_der_A_x,test_basis_der_A_y);    
        
        
        trial_basis_der_A_x = 1;trial_basis_der_A_y = 0;
        test_basis_der_A_x  = 1;test_basis_der_A_y  = 0;
        A3 = assemble_matrix_2D(matrix_size_A, gauss_point, gauss_weight, 'function_c3', P, T, T_trial, T_test,trial_basis_type, trial_basis_der_A_x,trial_basis_der_A_y, test_basis_type, test_basis_der_A_x,test_basis_der_A_y);    
        
        
        trial_basis_der_A_x = 0;trial_basis_der_A_y = 1;
        test_basis_der_A_x  = 0;test_basis_der_A_y  = 1;
        A4 = assemble_matrix_2D(matrix_size_A, gauss_point, gauss_weight, 'function_c4', P, T, T_trial, T_test, trial_basis_type, trial_basis_der_A_x,trial_basis_der_A_y, test_basis_type, test_basis_der_A_x,test_basis_der_A_y);    
        

        trial_basis_der_A_x_p = 0;trial_basis_der_A_y_p = 0;
        trial_basis_der_A_x_q = 0;trial_basis_der_A_y_q = 0;
        trial_basis_der_A_x = 0;trial_basis_der_A_y = 0;
        test_basis_der_A_x  = 0;test_basis_der_A_y  = 0;
        At = assemble_matrix_2D_ns_ns(matrix_size_A, gauss_point, gauss_weight, 'function_feval_c1','function_gudu',XL,Xm(:,n), P, T, T_trial, T_test, trial_basis_type_c, trial_basis_der_A_x_p,trial_basis_der_A_y_p,trial_basis_type_c, trial_basis_der_A_x_q,trial_basis_der_A_y_q,trial_basis_type, trial_basis_der_A_x,trial_basis_der_A_y, test_basis_type, test_basis_der_A_x,test_basis_der_A_y);    
        
        A=A1/dt+A2/dt+A3+A4+At/dt;
            

        trial_basis_der_b_x = 0;trial_basis_der_b_y = 0;
        test_basis_der_b_x = 0;test_basis_der_b_y = 0;
        b1 = assemble_vector_2D_ns(vector_size_b, gauss_point, gauss_weight, 'function_feval_c2','function_gu', XL,trial_basis_type, trial_basis_der_b_x,trial_basis_der_b_y,P, T,P_trial,T_trial, T_test, test_basis_type, test_basis_der_b_x,test_basis_der_b_y);
        
        trial_basis_der_b_x = 1;trial_basis_der_b_y = 0;
        test_basis_der_b_x = 1;test_basis_der_b_y = 0;
        b2 = assemble_vector_2D_ns_c(vector_size_b, gauss_point, gauss_weight, 'function_feval_c3','function_c3', XL,trial_basis_type, trial_basis_der_b_x,trial_basis_der_b_y,P, T,P_trial,T_trial, T_test, test_basis_type, test_basis_der_b_x,test_basis_der_b_y);
        
        trial_basis_der_b_x = 0;trial_basis_der_b_y = 1;
        test_basis_der_b_x = 0;test_basis_der_b_y = 1;
        b3 = assemble_vector_2D_ns_c(vector_size_b, gauss_point, gauss_weight, 'function_feval_c3','function_c4', XL,trial_basis_type, trial_basis_der_b_x,trial_basis_der_b_y,P, T, P_trial,T_trial,T_test, test_basis_type, test_basis_der_b_x,test_basis_der_b_y);
        
        trial_basis_der_b_x = 0;trial_basis_der_b_y = 0;
        test_basis_der_b_x = 0;test_basis_der_b_y = 0;
        bt = assemble_vector_2D_ns_ns(vector_size_b, gauss_point, gauss_weight, 'function_feval_c1','function_gu', XL,Xm(:,n),trial_basis_type, trial_basis_der_b_x,trial_basis_der_b_y,P, T,P_trial,T_trial, T_test, test_basis_type, test_basis_der_b_x,test_basis_der_b_y);
        
        b=bn+b1/dt+b2+b3+bt/dt;

        [AL,bL] = dirichlet_boundary_process('function_dirichlet',A,b,tn,boundary_nodes,P_trial);
        
        dXL=AL\bL;
        XL=XL+dXL;
        
    end

    Xm(:,n+1)=XL;
end


ul=Xm;
for n=1:Nt+1

u(:,n)=function_u(P_trial(1,:)',P_trial(2,:)',(n-1)*dt);

end

%离散解析解
u_comparison_u=[u,full(ul)];

%无穷范数误差估计
error_max_u = max(abs(full(ul)-u));   

%L2范数误差估计
error_L2_u=zeros(1,Nt+1);
trial_basis_der_A_x=0;trial_basis_der_A_y=0;
for n=1:Nt+1
error_L2_u(n)=error_L2_assemble_2D(u(:,n),ul(:,n),gauss_point, gauss_weight,P, T, T_trial, T_test, trial_basis_type, trial_basis_der_A_x,trial_basis_der_A_y);
error_L2_u(n)=sqrt(error_L2_u(n));
end

%H1半范数误差估计
error_H1_u=zeros(1,Nt+1);
trial_basis_der_A_x1=1;trial_basis_der_A_y1=0;
trial_basis_der_A_x2=0;trial_basis_der_A_y2=1;
for n=1:Nt+1
error_H1_u(n)=error_L2_assemble_2D(u(:,n),ul(:,n),gauss_point, gauss_weight,P, T, T_trial, T_test, trial_basis_type, trial_basis_der_A_x1,trial_basis_der_A_y1)+error_L2_assemble_2D(u(:,n),ul(:,n),gauss_point, gauss_weight,P, T, T_trial, T_test, trial_basis_type, trial_basis_der_A_x2,trial_basis_der_A_y2);

end

error_H1_uh=sqrt(sum(error_H1_u)*dt);
error=[error_max_u(Nt+1);error_L2_u(Nt+1);error_H1_uh]
toc