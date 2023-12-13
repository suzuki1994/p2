clc;
clear;
format long;

A = [0, 1; -3, -2];
B = [0; 1];
C = [1, 0];
D = 0;
sys=ss(A,B,C,D);
figure(1)
step(sys)
[ysys,tsys]=step (sys);
mp=max(ysys)-ysys(length(ysys-1));

At = conj(A.');
Ct = conj(C.');

N = [Ct, At * Ct];

% Teste de observabilidade
zeta = 0.7;
wn = 4;
Mp = exp((-pi*zeta)/sqrt(1-zeta^2));
ts5 = 3 / (zeta*wn);

wd = wn * sqrt(1 - zeta^2);
tp = pi / wd;
s1 = -zeta*wn + 1i*wd;
s2 = -zeta*wn - 1i*wd;
s3 = 10 * real(s1); % polo não dominante

A_chapeu = [A, zeros(2, 1); -C, 0];
B_chapeu = [B; 0];
%contrabilidade
Mc = [A, B; -C, 0];
K_chapeu = acker(A_chapeu,B_chapeu, [s1, s2, s3]);
[autovalores, autovetores] = eig(A_chapeu - B_chapeu * K_chapeu);

K = [K_chapeu(1,1), K_chapeu(1,2)]; % chapeu(1), chapeu(2)
Ki = -K_chapeu(1,3);

AA = [A - B*K, B*Ki; -C, 0];
BB = [0; 0; 1];
CC = [1, 0, 0];
DD = 0;

ctr = ss(AA, BB, CC, DD); % Sistema Compensado

tfinal = 7;
Kmax = tfinal / 0.01 + 1;
t = linspace(0, tfinal, Kmax);
k = linspace(0, Kmax, Kmax);

% Resposta para entrada degrau
[Yc, Tc] = step(ctr, tfinal);

% Comportamento da saída y
mp = max(Yc) - Yc(end);
fprintf('Mp pratico = %f\n', mp);
fprintf('Mp teorico = %f\n', Mp);

% Tempo de acomodação para 5%
i = length(Tc);
delta = 0;
while delta < 0.05
    delta = abs(Yc(end) - Yc(i)) / Yc(end);
    ts = Tc(i);
    i = i - 1;
end
fprintf('ts5%% pratico = %f\n', ts);
fprintf('ts5%% teorico = %f\n', ts5);

% Tempo de pico
j = 1;
while Yc(j) < max(Yc)
    tp_pratico = Tc(j);
    j = j + 1;
end
fprintf('tp pratico = %f\n', tp_pratico);
fprintf('tp teorico = %f\n', tp);
figure(2)
plot(Tc, Yc); % Sistema original vs compensado

