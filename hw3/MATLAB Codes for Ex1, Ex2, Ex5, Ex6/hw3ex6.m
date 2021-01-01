% hw#3 ex.6
% defining given values:
x = [2.5 3.5 5 6 7.5 10 12.5 15 17.5 20];
y = [13 11 8.5 8.2 7 6.2 5.2 4.8 4.6 4.3];

ex6(x,y)

function [] = ex6(x,y)
%linearizing equation
xi = log(x);
yi = log(y);
multiplication = xi.*yi;
xisquared = xi.^2;

n = numel(yi); %getting number of elements in a array

%calculating a1 and a0
a1 = (n*sum(multiplication)-sum(xi)*sum(yi))/(n*sum(xisquared)-(sum(xi))^2);
a0 = (sum(yi)/n) - (a1*(sum(xi)/n));

a = exp(a0)
b = a1
finalvalue=a*9^b

end