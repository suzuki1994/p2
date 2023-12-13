clc;
clear;
format long;
C1 = 1000e-6;
C2 = C1;
L = 20e-3;
R1 = 5;
R2 = R1;
tf = 100e-3;

A = [-R1/C1, 0, 1/C1; 0, -1/(R2*C2), -1/C2; -1/L, 1/L, 0];
B = [0; 1; 0];
C = [1, 0, 0];
D = 0;

G = ss(A, B, C, D);

p=pole(G);
z=zero(G);
u1 = 4 * p(1);
u2 = 4 * p(2);
u3 = 4 * p(3);

K = acker(A', C', [u1, u2, u3])'
Ke=K';
wd = imag(u2);
wa = 30 * wd;
T = 2 * pi / wa

Kmax = 250;
k = linspace(0, Kmax, Kmax+1);

x = zeros(length(k), 3);
x_p = zeros(length(k), 3);
x_obs = zeros(length(k), 3);
x_obs_p = zeros(length(k), 3);
u = ones(1, length(k));
y = zeros(1, length(k));

for i = 2:length(k)
    x(i, :) = T*x_p(i-1, :) + x(i-1, :);
    x_obs(i, :) = T*x_obs_p(i-1, :) + x_obs(i-1, :);
    
    y(i) = C*x(i, :)' + D*u(i);
    
    x_p(i, :) = A*x(i, :)' + B*u(i);
    x_obs_p(i, :) = A*x_obs(i, :)' + B*u(i) + Ke'*(y(i) - C*x_obs(i, :)');
end

figure('Position', [100, 100, 1200, 800]);

subplot(3, 3, 1);
plot(T*k, x(:, 1), 'LineWidth', 8, 'DisplayName', 'Planta');
hold on;
plot(T*k, x_obs(:, 1), 'w', 'LineWidth', 3, 'DisplayName', 'Estimado');
hold off;
title('x1(t)');
legend;
grid on;

subplot(3, 3, 2);
plot(T*k, x(:, 2), 'LineWidth', 8, 'DisplayName', 'Planta');
hold on;
plot(T*k, x_obs(:, 2), 'w', 'LineWidth', 3, 'DisplayName', 'Estimado');
hold off;
title('x2(t)');
legend;
grid on;

subplot(3, 3, 3);
plot(T*k, x(:, 3), 'LineWidth', 8, 'DisplayName', 'Planta');
hold on;
plot(T*k, x_obs(:, 3), 'w', 'LineWidth', 3, 'DisplayName', 'Estimado');
hold off;
title('x3(t)');
legend;
grid on;

subplot(3, 3, 4);
plot(T*k, x_p(:, 1), 'LineWidth', 8, 'DisplayName', 'Planta');
hold on;
plot(T*k, x_obs_p(:, 1), 'w', 'LineWidth', 3, 'DisplayName', 'Estimado');
hold off;
title('x1_ponto(t)');
legend;
grid on;

subplot(3, 3, 5);
plot(T*k, x_p(:, 2), 'LineWidth', 8, 'DisplayName', 'Planta');
hold on;
plot(T*k, x_obs_p(:, 2), 'w', 'LineWidth', 3, 'DisplayName', 'Estimado');
hold off;
title('x2_ponto(t)');
legend;
grid on;

subplot(3, 3, 6);
plot(T*k, x_p(:, 3), 'LineWidth', 8, 'DisplayName', 'Planta');
hold on;
plot(T*k, x_obs_p(:, 3), 'w', 'LineWidth', 3, 'DisplayName', 'Estimado');
hold off;
title('x3_ponto(t)');
legend;
grid on;

subplot(3, 3, [7, 8]);
plot(T*k, y);
title('y(t)');
grid on;

subplot(3, 3, 9);
plot(T*k, u);
title('u(t)');
grid on;