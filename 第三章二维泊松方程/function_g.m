function r=function_g(x,y)

if x==-1
    r=-1.5*y.*(1-y).*exp(-1+y);
elseif x==1
    r=0.5*y.*(1-y).*exp(1+y);
elseif y==-1
    r=-2*x.*(1-x/2).*exp(x-1);
elseif y==1
    r=0;
end
end