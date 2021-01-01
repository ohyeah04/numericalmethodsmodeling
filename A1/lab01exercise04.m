function lab01exercise04()

syms x; 
value1 = 'To calculate sinx approximation, \ngive the value for n-th estimate (integer only):';
n = input(value1);
value2 = 'Give the value for x:';
x1 = input(value2);

func = sin(x);

t1 = taylor(func, 'Order', nthprime(n)+1);
t2 = taylor(func, 'Order', nthprime(n-1));
truevalue = sin(x1);

approximation1 = subs(t1,x,x1);
approximation1 = vpa(approximation1);

previousapproximation = subs(t2,x,x1);
previousapproximation = vpa(previousapproximation);

truepercenterror = ((truevalue - approximation1)/truevalue)*100;
estimateerror = ((approximation1 - previousapproximation)/approximation1)*100;

fprintf('\nTrue percent error (et) for the %d-th estimate is %.6f percent\n and approximated error (ea) %.6f percent\n', n, truepercenterror, estimateerror);

end