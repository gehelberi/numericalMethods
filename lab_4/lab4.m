clc; 
clear; 

n = 7; %�������
N = n^2;
eps = 1.e-5;
A = gallery('poisson', n);

f = rand(N,1);
k_max = 1000;

B = diag(diag(A)); %����� �����
invB = inv(B);

D = invB * A; %��������������� �������
g = invB * f; %������ ����� g 
x = zeros(N, 1); %������� �����������
r = D * x - g;
err = norm (r) / norm(g);
k1 = 0;

while (err > eps && k1 < k_max)
   Ar = A * r;
   tau = (r'*r)/((Ar')*r);
   x = x - tau*r;
   r = D*x - g;
   err = norm (r)/norm(g);
   k1 = k1 + 1;
   Err1(k1) = err;
end

B = triu(A); %����� �������
invB = inv(B);

D = invB * A;
g = invB * f;
x = zeros(N, 1);
k2 = 0;
r = D*x - g;
err = norm (r)/norm(g);

while (err > eps && k2 < k_max)
   Ar = A * r;
   tau = (r' * r)/((Ar') * r);
   x = x - tau * r;
   r = D * x - g;
   k2 = k2 + 1;
   err = norm (r)/norm(g);
   Err2(k2) = err;
end

R1 = zeros(N); %�����������-����������� �����
for k = 1:N
   for m = 1:N
       if (k > m)
          R1(k,m) = A(k, m);
       elseif (k == m)
          R1(k, m) = A (k, m) / 2;
       end
   end
end

R2 = zeros(N);
for k = 1:N
   for m = 1:N
       if (k < m)
          R2(k,m) = A(k, m);
       elseif (k == m)
          R2(k, m) = A (k, m) / 2;
       end
   end
end

w = 0.5; %���. ����. ��������
B = (eye(N) + w * R1) * (eye(N) + w * R2); %������������ �����. ������
invB = inv(B);

D = invB * A;
g = invB * f;
x = zeros(N, 1);
k3 = 0;
r = D * x - g;
err = norm (r)/norm(g);

while (err > eps && k3 < k_max)
   Ar = A*r;
   tau = (r'*r)/((Ar')*r);
   x = x - tau*r;
   r = D*x - g;
   k3 = k3 + 1;
   err = norm (r)/norm(g);
   Err3(k3) = err;
end

semilogy(1:k1, Err1, 1:k2, Err2, 1:k3, Err3)
legend('Jacoby','Zeidel','Triag');
grid
