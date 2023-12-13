% Representação do Sistema no Espaço de Estados
clc;
clear;
close all;
format long;

R1=51000;
R2=18000;
C1=0.0000001;
C2=0.00000068;

A=[0 1/C1;-1/(R1*R2*C2) (-1/(R1*C2))+(-1/(R2*C2))];
B=[0;1/(R1*R2*C2)];
C=[1 0];
D=0;

[NUM,DEN] = ss2tf(A,B,C,D)

% Resposta ao degrau da função original
sistema=ss(A,B,C,D);
figure()
step(sistema)
[Y,t,X] = step(sistema);
title('Resposta ao Degrau da Planta Original')
% Separando as partes
figure()
subplot(2,1,1)
plot(t,X(:,1))
title('x1')
subplot(2,1,2)
plot(t,X(:,2))
title('x2')

% Informacoes de Projeto
Mp = 0.1
ts5 = 0.025                     

% Calculo do zeta
syms zeta
zeta = solve(Mp == exp(-pi*(zeta/sqrt(1-zeta^2))), zeta);
zeta = eval(zeta(1))

% Calculo da frequencia natutral nao amortecida (Wn), dos polos dominantes
% de malha fechada.
syms wn
wn = 3/(zeta*ts5)

% Calculo da frequencia natural amortecida (Wd).
wd = wn*sqrt(1-zeta^2)

s1=-zeta*wn+j*wn*sqrt(1-zeta^2)
s2=-zeta*wn-j*wn*sqrt(1-zeta^2)
s3=-10*abs(s1)   % Autovalor adicional para o projeto, valor mais alto que a parte real de s1 para maior estabilidade.

% Autovalores (polos de malha fechada desejados)
u1=s1;
u2=s2;
u3=s3;

% Matriz de Controlabilidade
M=ctrb(A,B)
% Teste de controlabilidade
rank(M)  

% Projeto do controlador / servosistema
A_chapeu = [A zeros(2,1);-C 0];
B_chapeu = [B;0];

% Matriz de ganhos do controlador K_chapeu = [k1 k2 -ki]
K_chapeu=acker(A_chapeu,B_chapeu,[u1 u2 u3])

% Verificando
eig(A_chapeu-B_chapeu*K_chapeu)

K = [K_chapeu(1) K_chapeu(2)];
Ki = -K_chapeu(3);
AA = [A-B*K B*Ki;-C 0];
BB = [0;0;1];
CC = [1 0 0];
DD = 0;

% Resposta ao degrau do controlador e original
figure()
step(sistema)
hold on
step(AA,BB,CC,DD)
title('Resposta ao Degrau da Planta Original X Controlador')
legend('Planta Original', 'Controlador')
hold off

% Projeto do observador de ordem plena

% Polinomio característico desejado
Pd=conv([1 -u1],[1 -u2]);
alfa1=Pd(2);
alfa2=Pd(3);
% Polinomio característico original
P=poly(A);
a1=P(2);
a2=P(3);
% Matriz de observabilidade
N=[conj(C') conj(A')*conj(C')];
rank(N)
W=[a1 1;1 0];
% Matriz de Transformação
Q=inv(W*conj(N'));
% Matriz de Ganho do Observador
Ke=Q*[alfa2-a2;alfa1-a1]
eig(A-Ke*C)

figure()
t=0:0.00001:0.14;
AAA = [A zeros(length(A)); Ke*C A-Ke*C];
BBB = [B;B];
CCC = eye(2*length(A));
DDD = zeros(2*length(A),1);

step(sistema)
hold on
so = step(AAA,BBB,CCC,DDD,1,t);
plot(t,so(:,1),'y')
title('Resposta ao Degrau da Planta Original X Observador de Ordem Plena')
legend('Planta Original', 'Observador de Ordem Plena')
hold off

figure()
ax = gca
ax.YLim = [0 1.4];
hold on
so = step(AAA,BBB,CCC,DDD,1,t);
plot(t,so(:,1),'ko')
step(AA,BB,CC,DD,1,t);
sori = step(A,B,C,D,1,t);
plot(t,sori,'r','LineWidth',2)
title('Resposta ao Degrau da Planta Original X Observador X Controlador')
legend('Observador','Controlado','Original');
hold off


%Equações recursivas
T = 1e-4;
k = 0:1:(0.15/T); 
r = ones(1,length(k)/2); 
r = [r 1.5.*r];
E(1) = 0;            % Saida do integrador          
x1(1) = 0; 
x2(1) = 0; 
y(1) = 0;            % Sinal de saida
y_til(1) = 0;  
x1_ponto(1) = 0; 
x2_ponto(1) = 0; 
x1_til(1) = 0; 
x2_til(1) = 0; 
x1_til_ponto(1) = 0; 
x2_til_ponto(1) = 0; 
e(1) = 0;          % Erro
v(1) = 0; 
u(1) = 0;          % Sinal de controle          
E_ponto(1) = r(1)-y(1);

for i=2:length(k)
    x1(i)=T*x1_ponto(i-1)+x1(i-1);
    x2(i)=T*x2_ponto(i-1)+x2(i-1);
    x1_til(i) = T*x1_til_ponto(i-1)+x1_til(i-1);
    x2_til(i) = T*x2_til_ponto(i-1)+x2_til(i-1);
    E(i) =T*E_ponto(i-1)+E(i-1);
    y(i) = C(1)*x1(i);
    y_til(i)=C(1)*x1_til(i);
    E_ponto(i) = r(i)-y(i);
    v(i) = K(1)*x1(i)+K(2)*x2(i);
    u(i) = Ki*E(i)-v(i);
    x1_ponto(i) = A(1,1)*x1(i)+A(1,2)*x2(i);
    x2_ponto(i) = A(2,1)*x1(i)+A(2,2)*x2(i)+B(2)*u(i);
    e(i) = y(i)-y_til(i);
    x1_til_ponto(i) = A(1,1)*x1_til(i)+A(1,2)*x2_til(i)+B(1)*u(i)+Ke(1)*e(i);
    x2_til_ponto(i) = A(2,1)*x1_til(i)+A(2,2)*x2_til(i)+B(2)*u(i)+Ke(2)*e(i);
end
figure()
plot(k*T, u)
title('Ação de Controle')

figure()
% step(sistema, 0.075)
% hold on
plot(k*T,y,'*g')
title('Saida do Sistema')
xlim([0.075 0.15]);
%legend('Planta Original', 'Planta Controlada')
hold off

figure()
subplot(2,1,1)
plot(k*T,x1,'*g')
title('x1')
hold on
subplot(2,1,2)
plot(k*T,x2,'*g')
title('x2')
hold on

subplot(2,1,1)
plot(k*T,x1_til,'or')
hold on
subplot(2,1,2)
plot(k*T,x2_til,'or')
hold off

figure();
plot(k*T,E_ponto,'*','DisplayName','Erro');

figure();
subplot(2,2,1); plot(k*T,x1);
title('Comportamento variável x1');grid;
subplot(2,2,3); plot(k*T,x2);
title('Comportamento variável x2');grid;
subplot(2,2,2); plot(k*T,x1_til);
title('Comportamento variável x1t');grid;
subplot(2,2,4); plot(k*T,x2_til);
title('Comportamento variável x2t');grid;

