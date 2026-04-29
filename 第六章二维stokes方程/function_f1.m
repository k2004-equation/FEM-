function r = function_f1(x,y,nu)

r = -2*nu*x.^2-2*nu*y.^2-nu*exp(-y)+pi^2*cos(pi*x).*cos(2*pi*y);

end