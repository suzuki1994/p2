clear
clc
close all;
format long g;

% Matrizes do Sistema
R1=43000;
R2=18000;
C1=0.0000001;
C2=0.00000068;

A=[0 1/C1;-1/(R1*R2*C2) (-1/(R1*C2))+(-1/(R2*C2))];
B=[0;1/(R1*R2*C2)];
C=[1 0];
D=0;

figure(1)
sistema=ss(A,B,C,D);
[Y,t,X]=step(sistema);
step(sistema)
title('Resposta ao Degrau da Planta Original')
% Informacoes de Projeto
Mp = 0.14;
ts5 = 0.022;
% Calculo do zeta
syms zeta
zeta = solve(Mp == exp(-pi*(zeta/sqrt(1-zeta^2))), zeta);
zeta = eval(zeta(1));

% Calculo da frequencia natutral nao amortecida (Wn)
syms wn
wn = 3/(zeta*ts5);
% Calculo da frequencia natural amortecida (Wd).
wd = wn*sqrt(1-zeta^2);
% Autovalores (polos de malha fechada desejados)
s1=-zeta*wn+j*wn*sqrt(1-zeta^2);
s2=-zeta*wn-j*wn*sqrt(1-zeta^2);
s3=-10*abs(s1);   % Autovalor adicional para o projeto, valor mais alto que a parte real de s1 para maior estabilidade.

% Projeto do controlador
A_chapeu = [A zeros(2,1);-C 0];
B_chapeu = [B;0];

% Matriz de ganhos do controlador K_chapeu = [k1 k2 -ki]
K_chapeu=acker(A_chapeu,B_chapeu,[s1 s2 s3]);

% Verificando
eig(A_chapeu-B_chapeu*K_chapeu)

K = [K_chapeu(1) K_chapeu(2)];
Ki = -K_chapeu(3);
AA = [A-B*K B*Ki;-C 0];
BB = [0;0;1];
CC = [1 0 0];
DD = 0;

% Resposta ao degrau do controlador e original
figure(2)
step(sistema)
hold on
step(AA,BB,CC,DD)
title('Resposta ao Degrau da Planta Original X Controlador')
legend('Planta Original', 'Controlador')
hold off

% observador de ordem mínima
Aaa = A(1, 1);
Aab = A(1, 2);
Aba = A(2, 1);
Abb = A(2, 2);
Ba = B(1);
Bb = B(2);

% Matriz de Controlabilidade do sistema dual para ordem mínima
Mc = ctrb(Aab, Abb);

% Teste de Controlabilidade
rank_Mc = rank(Mc); % número de postos do sistema

% Polo do observador, 10 vezes mais rápido do que os polos do controlador
L=-2570.441299881;% valor de S3 obtido no arquivo Projeto2_controlador
% Ganho do observador
Ke = acker(Abb, Aab, L);
Ke = Ke(1); % Transformando em escalar
Achapeu = Abb - Ke * Aab;
Bchapeu = Achapeu * Ke + Aba - Ke * Aaa;
Fchapeu = Bb - Ke * Ba;
Cchapeu = [0; 1];
Dchapeu = [1; Ke];

%Equacoes recursivas
clear u
T = 1e-4;
k = 0:1:(0.15/T); 
uu=ones(1,length(k));
r = ones(1,length(k)/2); 
r = [r 1.5.*r];
x1(1)= 0;  
x2(1)= 0;  
E(1)= 0;
u(1)=0;
x1_ponto(1)=0;  
x2_ponto(1)=0;  
y(1)=0;  
E_ponto(1)= r(1)-y(1);

x1_obs(1)=0;
x2_obs(1)=0;
x1_ponto_obs(1)=0;
x2_ponto_obs(1)=0;
eta1(1)=0;
eta1_ponto(1)=0;
xtil1(1)=0;
xtil2(1)=0;
y_obs(1) = 0;
for j=2:length(k)
    x1(j)=T*x1_ponto(j-1)+x1(j-1);
    x2(j)=T*x2_ponto(j-1)+x2(j-1); 
    x1_obs(j) = T * x1_ponto_obs(j - 1) + x1_obs(j - 1);
    x2_obs(j) = T * x2_ponto_obs(j - 1) + x2_obs(j - 1);
    
	x1_ponto_obs(j) = A(1,1)*x1_obs(j) + A(1,2)*x2_obs(j) + B(1)*uu(j);
    x2_ponto_obs(j) = A(2,1)*x1_obs(j) + A(2,2)*x2_obs(j) + B(2)*uu(j);  
    y_obs(j) = C(1) * x1_obs(j);
    eta1(j) = T * eta1_ponto(j - 1) + eta1(j - 1);
    eta1_ponto(j) = Achapeu * eta1(j) + Bchapeu * y_obs(j) + Fchapeu * uu(j);
    xtil1(j) = Cchapeu(1) * eta1(j) + Dchapeu(1) * y_obs(j);
    xtil2(j) = Cchapeu(2) * eta1(j) + Dchapeu(2) * y_obs(j);
       
    E(j)=T*E_ponto(j-1)+E(j-1);
    u(j)=-(K(1)*x1(j)+K(2)*x2(j))+Ki*E(j);
    x1_ponto(j)=A(1,1)*x1(j)+A(1,2)*x2(j)+B(1)*u(j);  
    x2_ponto(j)=A(2,1)*x1(j)+A(2,2)*x2(j)+B(2)*u(j); 
    y(j)=C(1)*x1(j); 
    E_ponto(j)= r(j)-y(j);
    
end
figure(3)
plot(k*T, u)
title('Ação de Controle')

figure(4)
plot(k*T,y,'*g')
title('Saída do Sistema')
xlim([0.075 0.15]);
hold off

figure(5);
plot(k*T,E_ponto);
title('Erro do Sistema')

figure(6);
subplot(2,2,1); plot(k*T,x1);
title('Comportamento variavel x1 controlador');grid;
subplot(2,2,3); plot(k*T,x2);
title('Comportamento variavel x2 controlador');grid;
subplot(2,2,2); plot(k*T,x1_obs);
title('Comportamento variavel x1 observador');grid;
subplot(2,2,4); plot(k*T,x2_obs);
title('Comportamento variavel x2 observador');grid;
