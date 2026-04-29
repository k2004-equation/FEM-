function r = function_f2(x,y,lambda,mu)
r = -(lambda+2*mu)*(2*x.*(x-1))....
    -(lambda+mu)*(pi^2*cos(pi*x).*cos(pi*y))-mu*(2*y.*(y-1));
end

