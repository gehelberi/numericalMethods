clear;
clc;

f=@(x) x.^4 - 8*x.^2 + 9*x - 3; % 1

tau = -1/41; %нахожу из уравнения 1+tau*f'(x)=0, за x беру решение уравнения x^4 - 8*x^2 + 9*x-3 = 0 (2)
eps = 10**(-11); 
fi=@(x) x+tau*f(x); % 3
X = [-0.5:0.1:3.5];
figure
plot(X, X, X, fi(X))

x(1) = 1,5; % 4
k=1; % 5
Err=1; 

while Err > eps % 6
  x(k+1) = fi(x(k));
  Err = abs(f(x(k+1)));
  k = k + 1;
endwhile

sprintf('answer x=%0.11f',x(k-1)) % 7 
figure % 8
semilogy(abs(f(x)));
