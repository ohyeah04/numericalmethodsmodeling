%lab #3; exercise #2
%using initial conditions from the given table in the lab, T and mu (u)

T = [0, 5, 10, 20, 30, 40];
u = [1.787, 1.519, 1.307, 1.002, 0.7975, 0.6529];

%this problem can be solved through Least Square Regression (page 336 from
%coursebook) calling a function
leastsquare(T,u)

%creating a function to solve this problem
function [] = leastsquare(T,u)
%Using least-square method
%given is T in celsius, changing it to Ta - in Kelvin: T+273.15
%so finding necessary data
Ta = T+273.15;
xi = 1./Ta;
yi = log(u);
multiplication = xi.*yi;
xisquared = xi.^2;

%numel - counting number of elements in the array
n = numel(yi);

%finding a1 and a0 (which are B and D)
a1 = (n*sum(multiplication)-sum(xi)*sum(yi))/(n*sum(xisquared)-(sum(xi))^2);
a0 = (sum(yi)/n) - (a1*(sum(xi)/n));

D = exp(a0)
B = a1

%final expression
y=D*exp(B./Ta);

%building a graph, plotting both data from table and approximated data hold
%to set up all parameters of graph
hold on
% plot - plotting a graph of arrays Ta and y, setting linewidth of graph
% to 1.25
plot(Ta, y, 'LineWidth', 1.25)

%stem - stemming data from the original table, Ta vs y, fill - filling in
%the datapoints with color, linestyle none - removes line from the axis to
%datapoints, color - setting the color of the datapoints.
stem(Ta, u, 'fill', '-.','LineStyle','none', 'Color', [0.8500, 0.3250, 0.0980])

%turning on grid
grid on

%labelling x and y axis
xlabel('Temperature')
ylabel('Dynamic Viscosity of Water')

%scaling the final graph
axis tight

%switching off hold, showing the graph
hold off

%end, happy with the code
end