f =@(x) 6+10*x+9*x^2+16*x^3;
bisection(f,-2,1)
function p = bisection(f,a,b)

if f(a)*f(b)>0 
    disp('Error, change limits!')
else
    p = (a + b)/2;
    xl = p;
    iter = 0;
    val = abs(f(p));
    while val > 0.00293
      if f(a)*f(p)<0 
         b = p;
      else
         a = p;         
      end
       p = (a + b)/2;
       iter=[iter; iter+1];
       xl = [xl; p];
       val = [val; abs(f(p))];
    end
end
table(xl, val)
end
    
    