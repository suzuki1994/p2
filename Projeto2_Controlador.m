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
% separando os elementos do vetor X
x1a=X(:,1);
x1b=X(:,2);


% Informacoes de Projeto
Mp = 0.14;
ts5 = 0.022;

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

% Autovalores (polos de malha fechada desejados)
s1=-zeta*wn+j*wn*sqrt(1-zeta^2)
s2=-zeta*wn-j*wn*sqrt(1-zeta^2)
s3=-10*abs(s1)   % Autovalor adicional para o projeto, valor mais alto que a parte real de s1 para maior estabilidade.

% Projeto do controlador
A_chapeu = [A zeros(2,1);-C 0];
B_chapeu = [B;0];

% Teste de controlabilidade
rank([A B;-C 0])

% Matriz de ganhos do controlador K_chapeu = [k1 k2 -ki]
K_chapeu=acker(A_chapeu,B_chapeu,[s1 s2 s3])

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

clear u
T = 1e-4;
k = 0:1:(0.15/T);
r=ones(1,length(k));

% condicoes iniciais
x1(1)= 0; 
x2(1)= 0;  
zeta(1)= 0;
u(1)=Ki*zeta(1);
x1_ponto(1)=B(1)*u(1);  
x2_ponto(1)=B(2)*u(1);  
y(1)=C(1)*x1(1)+C(2)*x2(1); 
zeta_ponto(1)= r(1)-y(1);

for j=2:length(k)
    % Equacoes dos integradores
    x1(j)=T*x1_ponto(j-1)+x1(j-1);
    x2(j)=T*x2_ponto(j-1)+x2(j-1);    

    zeta(j)=T*zeta_ponto(j-1)+zeta(j-1);

    % U=-K*X+Ki*zeta
    u(j)=-(K(1)*x1(j)+K(2)*x2(j))+Ki*zeta(j);
   
    % Equacao diferencial de estados:  Xponto=A*X+B*U
    x1_ponto(j)=A(1,1)*x1(j)+A(1,2)*x2(j)+B(1)*u(j);  
    x2_ponto(j)=A(2,1)*x1(j)+A(2,2)*x2(j)+B(2)*u(j); 
   
    % Equacao de Saaida: Y=C*X+D*U
    y(j)=C(1)*x1(j)+C(2)*x2(j);

    %Erro 
    zeta_ponto(j)= r(j)-y(j);


end
figure (3)
plot(t,Y,'b')
hold on
plot(k*T,y)
ylabel('y')
title('Resposta ao Degrau da Planta Original X Equação recursiva Controlador')
legend('Planta Original', 'Controlador')

hold off