clc;
clear;
format long;

A = [0 1; 0 -2];
B = [0; 1];
C = [1 0];
D = 0;

G = ss(A, B, C, D);
 K = [25, 5];
 Ke=[198;9604];
ki = K(1);

T = 0.01;
tfinal = 5;
Kmax = floor(tfinal / T) + 1;
k = linspace(0, Kmax, Kmax + 1);

ref = ones(1, length(k));

x1 = 0;
x2 = 0;
x1_ponto = 0;
x2_ponto = 0;
xtil1 = 0;
xtil2 = 0;
xtil1_ponto = 0;
xtil2_ponto = 0;
y = C(1,1) * x1;
ytil = C(1,1) * xtil1;
erro = 0;
u = 0;

for j = 2:Kmax+1
    % sistema original
    % Equações dos integradores
    x1(j) = T * x1_ponto(j-1) + x1(j-1);
    x2(j) = T * x2_ponto(j-1) + x2(j-1);    
    xtil1(j) = T * xtil1_ponto(j-1) + xtil1(j-1); 
    xtil2(j) = T * xtil2_ponto(j-1) + xtil2(j-1); 
    % equação de saída do observador de ordem mínima
    % Equação de Saída: Y=C*X+D*U
    y(j) = C(1,1) * x1(j);
    ytil(j) = C(1,1) * xtil1(j);
    erro(j) = ref(j) - xtil1(j);
    u(j) = -(K(1,2) * xtil2(j)) + (ki * (erro(j)));
    % Equação diferencial de estados:  Xponto=A*X+B*U
    x1_ponto(j) = A(1,1) * x1(j) + A(1,2) * x2(j) + B(1) * u(j);  
    x2_ponto(j) = A(2,1) * x1(j) + A(2,2) * x2(j) + B(2) * u(j); 
    xtil1_ponto(j) = A(1,1) * xtil1(j) + A(1,2) * xtil2(j) + B(1) * u(j) + Ke(1) * (y(j) - ytil(j));  
    xtil2_ponto(j) = A(2,1) * xtil1(j) + A(2,2) * xtil2(j) + B(2) * u(j) + Ke(2) * (y(j) - ytil(j)); 
end
figure(1)
plot(T*k, y);
title('y(t)');
grid on;

figure (2);
plot(T*k, x1);
title('x1(t)');
grid on;

figure (3);
plot(T*k, x2);
title('x2(t)');
grid on;

figure (4);
plot(T*k, u);
title('u(t)');
grid on;
