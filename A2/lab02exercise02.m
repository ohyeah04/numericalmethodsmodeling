muller(@(x) (x^3-0.5*x^2+4*x-2), 0, 1, 0.2, 1000);

function muller(func, x0, x1, x2, maxit)

iter = 0;

while(1)
    
    h0 = x1 - x0;
    h1 = x2 - x1;
    
    d0 = (func(x1) - func(x2))/h0;
    d1 = (func(x2) - func(x1))/h1;
    
    a = (d1 - d0)/(h1 + h0);
    b = a*h1 + d1;
    c = func(x2);
    
    rad = sqrt(b*b - 4*a*c);
    iter = iter + 1;
    
    if abs(b + rad) >= abs(b - rad)
        den = b + rad;
    else 
        den = b - rad;
    end
    dxr = -2*c/den;
    xr = x2 + dxr;
    
    if (abs(dxr) <= eps*xr || iter >= maxit), break, end
    
    x0=x1;
    x1=x2;
    x2=xr;
end

fprintf('The real positive root is %.2f\n', x0)

end