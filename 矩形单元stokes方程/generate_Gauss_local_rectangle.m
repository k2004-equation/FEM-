function [Gauss_weight_local_rectangle,Gauss_point_local_rectangle]=generate_Gauss_local_rectangle(Gauss_coefficient_reference_rectangle,Gauss_point_reference_rectangle,vertices_rectangle)

x1=vertices_rectangle(1,1);
x2=vertices_rectangle(1,2);
x3=vertices_rectangle(1,3);
x4=vertices_rectangle(1,4);%新加节点
y1=vertices_rectangle(2,1);
y2=vertices_rectangle(2,2);
y3=vertices_rectangle(2,3);
y4=vertices_rectangle(2,4);

h1=x3-x1;h2=y2-y1;
J=h1*h2/4;

%参考区间变换
Gauss_weight_local_rectangle=Gauss_coefficient_reference_rectangle*J;

%参考区间变换
Gauss_point_local_rectangle(:,1)=h1/2*Gauss_point_reference_rectangle(:,1)+x1+h1/2;
Gauss_point_local_rectangle(:,2)=h2/2*Gauss_point_reference_rectangle(:,2)+y1+h2/2;