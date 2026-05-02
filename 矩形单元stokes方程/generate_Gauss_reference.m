function [Gauss_coefficient_reference,Gauss_point_reference]=generate_Gauss_reference(Gauss_point_number)

n=0;
Gauss_coefficient=[5/9;8/9;5/9];
Gauss_point=[-sqrt(3/5);0;sqrt(3/5)];

Gauss_coefficient_reference=zeros(1,9);
Gauss_point_reference=zeros(9,2);
for i=1:3
    for j=1:3
        n=n+1;
    Gauss_coefficient_reference(1,n)=Gauss_coefficient(i).*Gauss_coefficient(j);
    Gauss_point_reference(n,1)=Gauss_point(i);
    Gauss_point_reference(n,2)=Gauss_point(j);
    end
end
end