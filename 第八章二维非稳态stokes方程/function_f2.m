function r = function_f2(x,y,t,nu)

r = -2*pi*(-2/3*x.*y.^3+2-pi*sin(pi*x)).*sin(2*pi*t)+(4*nu*x.*y-nu*pi^3*sin(pi*x)+2*pi*(2-pi*sin(pi*x)).*sin(2*pi*y)).*cos(2*pi*t);

end

