function r=function_g(x,y,t)

if x==0
    r=exp(y+t);
elseif x==2
    r=exp(2+y+t);
elseif y==0
    r=exp(x+t);
elseif y==1
    r=exp(x+1+t);
end
end