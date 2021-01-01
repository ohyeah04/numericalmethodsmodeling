c = [0.5, 0.8, 1.5, 2.5, 4];
k = [1.1, 2.4, 5.3, 7.6, 8.9];
inversek = 1./k;
inversec = 1./(c.^2);

ex2(inversec,inversek, c, k)

function [] = ex2(x,y, c, k)
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

kmax = 1/a0
cs = a1*kmax

hold on
stem(c, k, 'fill', '-.','LineStyle','none', 'Color', [0.8500, 0.3250, 0.0980])
fplot(@(c) (kmax*c^2)/(cs+c^2),[0 5])
grid on
axis tight
hold off

end