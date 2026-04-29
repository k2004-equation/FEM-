function r=function_g2(x,y,t)

if x==0
    r=2*cos(2*pi*t);
elseif x==1
    r=(-2/3*y.^3+2).*cos(2*pi*t);
elseif y==-0.25
    r=(1/96*x+2-pi*sin(pi*x))*cos(2*pi*t);
elseif y==0
    r=(2-pi*sin(pi*x))*cos(2*pi*t);
end
end