x = [1 2 3 4 5 6 7 8 9];
y = [1 1.5 2 3 4 5 8 10 13];

ex5a(x,y)

function [] = ex5a(x,y)
%we have linear eq. y = a1*x + a0
%so again we need to calcualte a0 and a1 
xi = x;
yi = y;
multiplication = xi.*yi;
xisquared = xi.^2;
yisquared = yi.^2;

n = numel(yi); %getting number of elements in an array

a1 = (n*sum(multiplication)-sum(xi)*sum(yi))/(n*sum(xisquared)-(sum(xi))^2)
a0 = (sum(yi)/n) - (a1*(sum(xi)/n))

%approximating the error&finding correlation coefficient
Sr=sum((yi-a0-a1*xi).^2)
Syx=sqrt(Sr/7)

r=(n*sum(multiplication)-sum(xi)*sum(yi))/(sqrt(n*sum(xisquared)-(sum(xi))^2)*sqrt(n*sum(yisquared)-(sum(yi))^2))

hold on
stem(x, y, 'fill', '-.','LineStyle','none', 'Color', [0.8500, 0.3250, 0.0980])
fplot(@(k) a1*k+a0,[0 10])
grid on
axis tight
hold off

end