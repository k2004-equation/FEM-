function r = function_f1(x,y,t,nu)

r = -2*pi*(x.^2.*y.^2+exp(-y)).*sin(2*pi*t)+(-2*nu*x.^2-2*nu*y.^2-nu*exp(-y)+pi^2*cos(pi*x).*cos(2*pi*y)).*cos(2*pi*t);

end