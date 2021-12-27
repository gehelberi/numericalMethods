clear;
clc;

n = 10;
v = 1;

for i = 1:n
  a = v + 10 / (10 + i);
  b = v + 10 / i;
  c = v + (1 + i) / 10;
endfor

A = zeros(n);                                                % 1
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

copy_A = A;
S = eye(n);

for i = 1:n - 1
  M = eye(n);
  M_inv = eye(n);
  if A(n - i + 1, n - i) != 0
    for j = 1:n
      if i + j == n
        M(n - i, j) = 1 / A(n - i + 1, n - i);
      elseif
        M(n - i, j) = -A(n - i + 1, j) / A(n - i + 1, n - i);
      endif
      M_inv(n - i, j) = A(n - i + 1, j);
    endfor
  endif
  A = M_inv * A * M;
  S = S * M;
endfor

M


eq = (-1)**n * [1 (-1 * A(1, :))];                           % 2

eig_values = roots(eq);                                      % 3
eig_values

y = zeros(n, n);                                             % 4
eig_vectors = zeros(n, n);

for i = 1:n
  for j = 1:n
    y(j, i) = eig_values(i) ** (n - j);                      % 13.10 - 4) y
  endfor
  eig_vectors(:, i) = S * y(:, i);                           % 13.11 - 4) x
endfor

y


eig_vectors

display("------------check-------------");                   % 5
[vectors, D] = eigs(copy_A, n);
values = diag(D);
values
vectors

for i = 1:n                                                  % 6
  display(norm(copy_A * eig_vectors(:, i) - eye(n) * eig_values(i) * eig_vectors(:, i)));
endfor 

