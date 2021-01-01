NewInt([0, 1.8, 5, 6, 8.2, 9.2, 12], [26, 16.415, 5.375, 3.5, 2.015, 2.54, 8], 3.5);

function hw4ex3 = NewInt(x,y,p)
%according to the pseudocode on page 498 of the book, 
% writing a matlab code

format long g

n = length(x); %getting a size of array 
a(1) = y(1); %equating two arrays 

for k = 1 : n - 1
   fdd(k, 1) = (y(k+1) - y(k))/(x(k+1) - x(k));
end
for j = 2 : n - 1
   for k = 1 : n - j
      fdd(k, j) = (fdd(k+1, j - 1) - fdd(k, j - 1))/(x(k+j) - x(k)) %evaluating values
   end
end
for j = 2 : n
   a(j) = fdd(1, j-1);
end
Df(1) = 1;
c(1) = a(1);
for j = 2 : n
   Df(j)=(p - x(j-1)) .* Df(j-1);
   c(j) = a(j) .* Df(j);
end
result=sum(c);%calculating the final answer

fprintf("Calculated Value: %f\n",result)
end




