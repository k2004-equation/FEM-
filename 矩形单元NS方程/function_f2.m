function r = function_f2(x,y,nu)

r = 4*nu*x.*y-nu*pi^3*sin(pi*x)+2*pi*(2-pi*sin(pi*x)).*sin(2*pi*y)....
    +(x.^2.*y.^2+exp(-y)).*(-2/3*y.^3-pi^2*cos(pi*x))+(-2/3*x.*y.^3+2-pi*sin(pi*x)).*(-2*x.*y.^2);

end

