clear;
clc;

n = 10;
v = 1;
eps = 10**(-9);

function[sum] = findNotDiagProd(A, n)
	A = A - diag(diag(A), 0);
	sum = 0;
	for i = 1:n
		for j = 1:n
			sum = sum + A(i, j)**2;
		endfor
	endfor
endfunction

A = zeros(n);
for i = 1:n
  a = v + 10 / (10 + i);
  b = v + 10 / i;
  c = v + (1 + i) / 10;
  
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

copy_A = A;

display(A);

while(findNotDiagProd(A, n) > eps)
	for k = 1:n
		for m = 1:n
			if(k != m)
				T = eye(n);
				w = (2 * A(k, m)) / (A(k, k) - A(m, m));
				s = sign(w) * sqrt((sqrt(1 + w**2) - 1) / (2 * sqrt(1 + w**2)));   %(2)
				c = sqrt((sqrt(1 + w**2) + 1) / (2 * sqrt(1 + w**2)));             %(2)
				T(k, k) = c;
				T(m, m) = c;
				T(m, k) = s;
				T(k, m) = -s;
				A = T' * A * T;
			endif
		endfor
	endfor
endwhile

display(A);
display(diag(A));
display("======================================================check======================================================")
display(eigs(copy_A, n));
