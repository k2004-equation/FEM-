function r = function_fn(x,y,t)

r = (t.^2.*x.^2.*(x-1).^2.*y.^2.*(y-1).^2+1).*x.*(x-1).*y.*(y-1)-t*((4*x+1).*y.*(y-1)+x.*(x-1).*(4*y+1));

end