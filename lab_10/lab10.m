clear;
clc;

f=@(x) x.^4 - 8*x.^2 + 9*x-3;     % 1
df=@(x) 4*x.^3-16*x+9;
eps = 10**(-11); 
fi=@(x) x-(f(x)/df(x));           % 2
X = [-0.5:0.1:3.5];

figure
plot(X, f(X))

x(1) = 2;                         % 3
k=1;                              % 4
Err=1;

while Err > eps && k<1000         % 5
  x(k+1) = fi(x(k));
  Err = abs(f(x(k+1)));
  k = k + 1;
endwhile

sprintf('answer x=%0.11f',x(k-1)) % 6
figure
semilogy (abs(f(x)));             % 7
