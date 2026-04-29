function [gauss_point_local,gauss_weight_local]=Gauss_local_values_1D(vertices, gauss_point, gauss_weight)  ;
xl = vertices(1);
xr = vertices(2);

%参考单元高斯区间[-1,1]，局部单元区间为[xl,xr]，需建立参考单元到局部单元的映射
gauss_point_local  = xl + (xr - xl)*(1+gauss_point)/2;
gauss_weight_local = gauss_weight*(xr - xl)/2;

end

