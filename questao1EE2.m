clc;
clear ;

% Declaração de variáveis
K = 0.2;
C = 1000e-6;
L = 10e-3;
R = 5;

% Matrizes de estados cálculadas 
% X1 = Vc1; X2 = iL
A = [-1/(R*C), 1/C; -1/L, 0];
B = [K/(R*C); 1/L];
C = [1, 0];
D = K;

% Sistema no espaço de estados
sys = ss(A, B, C, D);

% Convertendo para função de transferência diretamente
G = tf(sys);
[num, den] = tfdata(G, 'v');
figure(1)
step(sys)

%recursiva
Ts = 1e-3;
Fs = 1 / Ts;
tfinal = 0.1;

Kmax = round(tfinal / Ts) + 1;
t = linspace(0, tfinal, Kmax);

% declaração dos vetores
u = ones(1, Kmax);
x1 = zeros(1, Kmax);
x2 = zeros(1, Kmax);
x1_ponto = zeros(1, Kmax);
x2_ponto = zeros(1, Kmax);
y = zeros(1, Kmax);

% Condições iniciais
% Equação diferencial de estados: Xponto = A*X + B*U
x1_ponto(1) = A(1,1)*x1(1) + A(1,2)*x2(1) + B(1)*u(1);
x2_ponto(1) = A(2,1)*x1(1) + A(2,2)*x2(1) + B(2)*u(1);

% Equação de Saída: Y = C*X + D*U
y(1) = C(1,1)*x1(1) + C(1,2)*x2(1) + D*u(1);

for j = 2:Kmax
    % Equações dos integradores
    x1(j) = Ts*x1_ponto(j-1) + x1(j-1);
    x2(j) = Ts*x2_ponto(j-1) + x2(j-1);

    % Equação diferencial de estados: Xponto = A*X + B*U
    x1_ponto(j) = A(1,1)*x1(j) + A(1,2)*x2(j) + B(1)*u(j);
    x2_ponto(j) = A(2,1)*x1(j) + A(2,2)*x2(j) + B(2)*u(j);

    % Equação de Saída: Y = C*X + D*U
    y(j) = C(1,1)*x1(j) + C(1,2)*x2(j) + D*u(j);
end

figure(2);
plot(t, x1);
title('Resposta do estado x1 da planta ao degrau unitário');
xlabel('t');
ylabel('x1(t)');
legend('x1(t)');
grid on;

figure(3);
plot(t, x2);
title('Resposta do estado x2 da planta ao degrau unitário');
xlabel('t');
ylabel('x2(t)');
legend('x2(t)');
grid on;