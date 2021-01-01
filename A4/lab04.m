disp("Solution for 18.1(a):")
NewInt([8,12], [0.9030900,1.0791812], 10);
disp("===========================")
disp("Solution for 18.1(b):")
NewInt([9,11], [0.9542425,1.0413927], 10);
disp("===========================")
disp("Solution for 18.2:")
NewInt([8,9,11], [0.9030900,0.9542425,1.0413927], 10);
disp("===========================")
disp("Solution for 18.3:")
NewInt([8,9,11,12], [0.9030900,0.9542425,1.0413927,1.0791812], 10);

function lab4 = NewInt(x,y,p)
%according to the pseudocode on page 498 of the book, 
% writing a matlab code


n = length(x); %getting a size of array 
a(1) = y(1); %equating two arrays 

for k = 1 : n - 1
   fdd(k, 1) = (y(k+1) - y(k))/(x(k+1) - x(k));
end

%evaluating values, if we uncomment this, here we draw a table with
%coefficients
for j = 2 : n - 1
   for k = 1 : n - j
      fdd(k, j) = (fdd(k+1, j - 1) - fdd(k, j - 1))/(x(k+j) - x(k)); 
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

%now calculatng error and final value (Calculated value)
lab4=sum(c);
trval=log10(p);%true value
err=((trval-lab4)/trval)*100;

fprintf("Calculated Value: %f\nTrue Value: %f\nError: %f%%\n",lab4,trval,err)
end




