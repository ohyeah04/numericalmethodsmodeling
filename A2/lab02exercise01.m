newtraph(@(x) (x^3-5*x^2+7*x-3), @(x) (3*x^2-10*x+7), @(x) (6*x-10), 0.01, 0.0001, 3);

function[root,ea,iter] = newtraph(func,dfunc,ddfunc,xr,es,maxit)
if nargin < 3, error(''); end
if nargin < 4 || isempty(es), es=0.0001; end
if nargin < 5 || isempty(maxit), maxit=3; end

iter = 3;
ea = 100;

er = xr;
errorpercentage = ea;
nof = iter;

while(1)
    xrprev = xr;
    xr = xr - func(xr)*dfunc(xr)/(dfunc(xr)^2 - func(xr)*ddfunc(xr));
    iter = iter - 1;
    if xr ~= 0, ea = abs((xr-xrprev)/xr) * 100; end
    
    er = [er; xr];
    errorpercentage = [errorpercentage; ea];
    nof = [nof; iter];
    
    if ea <= es || iter == 0, break, end
end

root = xr;

table(er,nof,errorpercentage)

end