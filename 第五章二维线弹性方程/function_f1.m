function r = function_f1(x,y,lambda,mu)
r = -(lambda+2*mu)*(-pi^2*sin(pi*x).*sin(pi*y))....
    -(lambda+mu)*((2*x-1).*(2*y-1))-mu*(-pi^2*sin(pi*x).*sin(pi*y));
end