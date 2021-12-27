clear;
clc;

eps = 0.001;

a = randi([-10,10]);
b = a + randi([-10, 10]);

v = 1;
n = 10;
A = zeros(n);
f = zeros(n ,1);

for i = 1:n
    f(i) = a + (b - a) * rand;
    for j = 1:n
        if i == j
            A(i, j) = 100 + v;
        else 
            A(i, j) = 1 / (i + j + v);
        endif
    endfor
endfor


norm_A = norm(A);
tau = [1 / (2 * norm_A), 1 / (4 * norm_A), 1 / (8 * norm_A)];

er = zeros(size(tau)(2), 1);

for i = 1:size(tau)(2)
   
    H = eye(n) - tau(i) * A;
    phi = tau(i) * f;
    q = norm(H);
    xs = phi;
    xn = H * xs + phi;
    iter = 0;
    while q / (1-q) * norm(xs - xn) > eps
        xs = xn;
        xn = H * xs + phi;
        iter = iter + 1;
        er(i, iter) = norm(xs - xn) * q / (1 - q);
    endwhile
    
endfor

x = zeros(size(er)(2),1);
for i = 1:size(er)(2)
    x(i) = i;
endfor
plot(x, er(1, :), x, er(2, :),x , er(3, :));
legend ("tau=1/(2|A|)", "tau=1/(4|A|)", "tau=1/(8|A|)");

