function z = function_f(x,y)
z = -y.*(1-y).*(1-x-x.^2./2).*exp(x+y)-x.*(1-x./2).*(-3*y-y.^2).*exp(x+y);
end