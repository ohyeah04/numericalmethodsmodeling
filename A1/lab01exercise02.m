[h, ea] = newtraph(@(h) 30 - pi*h.^2.*(3*3-h)/3, @(h) -2*3*pi*h + pi*h.^2, 3, 0.01, 4);

function [root, ea, iter] = newtraph(func, dfunc, xr, es, maxit)

if nargin<3, error('at least 3 input arguments required'),end
if nargin<4, es = 0.0001; end
if nargin<5, maxit = 50; end

iter = 0;

fprintf('Iter#  |   Value   |   approx. rel. err \n');
fprintf(' %d         %.4f           - \n', iter, xr);

while iter < maxit
    x1 = xr-func(xr)/dfunc(xr);
    ea = 100*abs((x1-xr)/x1);
    iter = iter+1;
    fprintf(' %d         %.4f        %.4f%%\n', iter, x1, ea);
    xr = x1;
    if ea<=es || iter >=maxit, break, end
end

fprintf('\nSo the result is iteration number %d\n', iter);
root = xr;

end