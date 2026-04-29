function r=function_g1(x,y,t)

if x==0
    r=exp(-y).*cos(2*pi*t);
elseif x==1
    r=(y.^2+exp(-y)).*cos(2*pi*t);
elseif y==-0.25
    r=(1/16*x^2+exp(0.25)).*cos(2*pi*t);
elseif y==0
    r=cos(2*pi*t);
end
end