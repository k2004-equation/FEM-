function r = function_f1(x,y,nu)

r = -2*nu*x.^2-2*nu*y.^2-nu*exp(-y)+pi^2*cos(pi*x).*cos(2*pi*y)....
    +2*x.*y.^2.*(x.^2.*y.^2+exp(-y))+(-2/3*x.*y.^3+2-pi*sin(pi*x)).*(2*x.^2.*y-exp(-y));

end