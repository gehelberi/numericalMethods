clear;
clc;

n = 10;
v = 1;
eps = 0.5 * 10**(-3);

z = ones(n, 1);
y = ones(n, 1);
A = zeros(n);

for i = 1:n
a = v + 10 / (10 + i);
b = v + 10 / i;
c = v + (1 + i) / 10;
endfor

for i = 1:n
  A(i, i) = a;
  
  if i <= n - 1
    A(i, i + 1) = b;
    A(i + 1, i) = b;
  endif
  
  if i <= n - 2
    A(i, i + 2) = c;
    A(i + 2, i) = c;
  endif
endfor

A

% степенной метод
y = A * z;                                          % 2
lambda = y ./ z;                                    % 3
max_lambda = max(lambda);
min_lambda = min(lambda);
while abs(max_lambda - min_lambda) > eps            % 4
  z = y / norm(y);
  y = A * z;
  lambda = y ./ z;
  max_lambda = max(lambda);
  min_lambda = min(lambda);
endwhile
x1 = z;

% метод скалярных произведений
lambda_s = 0;                                       % 1
y = A * z;                                          % 2
lambda_n = dot(y, z) / dot(z, z);                   % 3
while abs(lambda_s - lambda_n) > eps                % 4
  z = y / norm(y);
  y = A * z;
  lambda_s = lambda_n;
  lambda_n = dot(y, z) / dot(z, z);
endwhile
x2 = z;

[v, D] = eigs(A)

max_lambda
x1
lambda_n
x2

