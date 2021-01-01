x = [1 2 3 4 5 6 7 8 9];
y = [1 1.5 2 3 4 5 8 10 13];

ex5b(x,y)

function [] = ex5b(x,y)
%we have linear eq. y = a2*x^2+a1*x + a0
%so again we need to calcualte a0 and a1 
xi = x;
yi = y;
multiplication = xi.*yi;
multiplication2 = xi.^2.*yi;
xisquared = xi.^2;
xithree = xi.^3;
xifour = xi.^4;

n = numel(yi); %getting number of elements in an array

%calculating coefficients a0, a1, a2
syms a0 a1 a2
eqn1 = sum(xifour)*a0+sum(xithree)*a1+sum(xisquared)*a2 == sum(multiplication2);
eqn2 = sum(xithree)*a0+sum(xisquared)*a1+sum(xi)*a2 == sum(multiplication);
eqn3 = sum(xisquared)*a0+sum(xi)*a1+n*a2 == sum(yi);

[A,B] = equationsToMatrix([eqn1, eqn2, eqn3], [a0, a1, a2]);

solution = linsolve(A,B);
a2=solution(1)
a1=solution(2)
a0=solution(3)

%approximating the error&finding correlation coefficient
Sr=sum((yi-a0-a1*xi).^2)
Syx=sqrt(Sr/7)

%plotting based on data
hold on
stem(x, y, 'fill', '-.','LineStyle','none', 'Color', [0.8500, 0.3250, 0.0980])
fplot(@(k) a2*(k^2)+a1*k+a0,[0 10])
grid on
axis tight
hold off

end