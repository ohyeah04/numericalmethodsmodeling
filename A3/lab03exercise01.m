%lab #3; exercise #1
%applying initial conditions from the manual

y = @(x) (2.5/(120*50000*30000*600))*(-x.^5 + 2*600^2*x.^3 - 600^4*x);
xl = 0;
xu = 600;
es = 1;
[x, fx, iter] = goldensection(y, xl, xu, es);

%output
fprintf('Maximum deflection is point x = %f, \nAnd the maximum deflections is %f\n', x, abs(fx));

%creating a function goldensection
function [x,fx,ea,iter]=goldensection(f,xl,xu,es,maxit,varargin)
    %conditions - not enough arguements - end code
    if nargin<3,error('at least 3 input arguments required'),end
    if nargin<4 || isempty(es), es=0.0001;end
    if nargin<5 || isempty(maxit), maxit=50;end
    phi=(1+sqrt(5))/2;
    iter=0;
    
    %applying golden section search method, based on the pseudocode
    %provided
    while(1)
        d = (phi-1)*(xu - xl);
        x1 = xl + d;
        x2 = xu - d;
        if f(x1,varargin{:}) < f(x2,varargin{:})
          xopt = x1;
          xl = x2;
        else
          xopt = x2;
          xu = x1;
        end
        iter=iter+1;
        if xopt~=0, ea = (2 - phi) * abs((xu - xl) / xopt) * 100;end
        if ea <= es ||  iter >= maxit,break,end
    end
    
    %output
    x=xopt;fx=f(xopt,varargin{:});
end

