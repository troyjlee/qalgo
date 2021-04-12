function [a,p,q] = cfr(z, tol)
  a = floor(z);
  y = z - a;
  while abs(y) > tol
    a = [a;floor(1/y)];
    y = 1/y - a(end);
  endwhile 

  % now compute vectors p, q
  % such the p./q gives the value 
  % of convergents
  n = size(a,1);
  x = a(n);
  for i=n-1:-1:1
    x = a(i) + 1/x;  
  endfor

  p = zeros(n,1);
  q = ones(n,1);

  p(1) = a(1);
  p(2) = a(2)*a(1)+1;
  q(1) = 1;
  q(2) = a(2);
  for i=3:n
    p(i) = a(i)*p(i-1)+p(i-2);
    q(i) = a(i)*q(i-1) + q(i-2);
  end 
