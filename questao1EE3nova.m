clc;
clear ;

A=[-5 1;-45 0];
B=[0; 45];
C=[1 0];
D=0;

G=ss(A,B,C,D);
G1=tf(G);

zeta= 0.6;
Ts=0.5;
wn=4/(zeta*Ts)

% Verificando polos da planta
autovalores_planta = eig(A);
polos_planta = autovalores_planta

s1=-zeta*wn+1j*wn*sqrt(1-zeta^2);
s2= conj(s1);

wd = imag(s1);
wa = 60*wd;
T =2*pi/wa %s
T1 = (2*pi/wa)*1000; %ms

Mc= [B A*B]

A_c = [A zeros(2,1); -C 0];
B_c = [B; 0];

K_c = acker(A_c, B_c, [s1, s2, -5*zeta*wn])
K = [K_c(1, 1), K_c(1, 2)]
Ki = -K_c(1, 3)

Kmax = 250;
k = linspace(0, Kmax, Kmax+1);

 ref = ones(1, length(k));
 x = zeros(2, length(k));
 x_p = zeros(2, length(k));

 qsi = zeros(1, length(k));
 qsi_p = zeros(1, length(k));
 u = zeros(1, length(k));
 y = zeros(1, length(k));

qsi_p(1) = ref(1) - y(1);

for i = 2:length(k)
    x(:, i) = T*x_p(:, i-1) + x(:, i-1);
    qsi(i) = T*qsi_p(i-1) + qsi(i-1);
    
    u(i) = qsi(i)*Ki - K*x(:, i);
    y(i) = C*x(:, i) + D*u(i);
    
    x_p(:, i) = A*x(:, i) + B*u(i);
    qsi_p(i) = ref(i) - y(i);
end

figure('Position', [100, 100, 1200, 800]);

subplot(3, 2, 1);
plot(T*k, x(1, :));
title('x1(t)');
grid on;

subplot(3, 2, 2);
plot(T*k, x(2, :));
title('x2(t)');
grid on;

subplot(3, 2, 3);
plot(T*k, x_p(1, :));
title('x1_ ponto(t)');
grid on;

subplot(3, 2, 4);
plot(T*k, x_p(2, :));
title('x2_ ponto(t)');
grid on;

subplot(3, 2, 5);
plot(T*k, y);
title('y(t)');
grid on;

subplot(3, 2, 6);
plot(T*k, u);
title('u(t)');
grid on;