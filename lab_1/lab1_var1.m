clear;
clc;

n = 5;
A = zeros(n, n);
f = zeros(n, 1);
E = eye(size(A));

res = zeros(n, n);

for i = 1:n
    f(i) = 1/i;
    for j = 1:n
      if (i == j)
            A(i, j) = 1 + 1 / j;
        else
            A(i, j) = (n - j)^2;
        endif
    endfor
endfor

function [x, d] = Gauss (A,f)
    d = 1;
    [n, m] = size(A);
    x = zeros(n, 1);
    A(:, m + 1) = f;
    for k = 1:n
        d = d * A(k, k);
        A(k, :) = A(k, :) / A(k, k); 
        for j = k + 1:m
            A(j, :) = A(j, :) - A(j, k) * A(k, :);
        endfor
    endfor
    f = A(:, m + 1);
    x(n) = f(n);
    for i = n-1:-1:1
        x(i) = f(i) - A(i,i + 1:end - 1) * x(i+1:end);
    endfor
endfunction

[x, d] = Gauss(A, f);
x
A \ f
d
det(A)  
for i = 1:n
    [res(:, i), _] = Gauss(A, E(:, i));
endfor
res
inv(A)


