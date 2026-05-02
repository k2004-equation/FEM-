function r=function_g1(x,y)

if x==0
    r=exp(-y);
elseif x==1
    r=y^2+exp(-y);
elseif y==-0.25
    r=1/16*x^2+exp(0.25);
elseif y==0
    r=1;
end
end