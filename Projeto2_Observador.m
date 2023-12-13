% Projeto de um observador de estados de ordem mínima

clc;
clear;
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
sistema=ss(A,B,C,D);

% resposta ao degrau
[Y,t,X] = step(sistema);
% separando as saídas
x1a=X(:,1);
x1b=X(:,2);

figure(1)
subplot(2,2,1)
plot(t,x1a)
hold on
subplot(2,2,2)
plot(t,x1b)
hold on
subplot(2,2,3)
plot(t,Y)
hold on

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
Ke = Ke(1) % Transformando em escalar

Achapeu = Abb - Ke * Aab
Bchapeu = Achapeu * Ke + Aba - Ke * Aaa
Fchapeu = Bb - Ke * Ba

Cchapeu = [0; 1]
Dchapeu = [1; Ke]



%recursiva 
clear u
T = 1e-4;
k=0:0.15/T;
uu=ones(1,length(k));

% definindo os vetores
x1_obs(1)=0;
x2_obs(1)=0;
x1_ponto_obs(1)=0;
x2_ponto_obs(1)=0;
eta1(1)=0;

eta1_ponto(1)=Fchapeu(1)*uu(1);
xtil1(1)=0;
xtil2(1)=0;
% condições iniciais (podem ser diferentes de zero)
y_obs(1) = C(1,1)*x1_obs(1) + C(1,2)*x2_obs(1) + D*uu(1);  


for j=2:length(k)
    % sistema original
    % Equações dos integradores
    x1_obs(j) = T * x1_ponto_obs(j - 1) + x1_obs(j - 1);
    x2_obs(j) = T * x2_ponto_obs(j - 1) + x2_obs(j - 1);
    
    % Equação diferencial de estados:  Xponto=A*X+B*U
	x1_ponto_obs(j) = A(1,1)*x1_obs(j) + A(1,2)*x2_obs(j) + B(1)*uu(j);
    x2_ponto_obs(j) = A(2,1)*x1_obs(j) + A(2,2)*x2_obs(j) + B(2)*uu(j);  
    
    % Equação de Saída: Y=C*X+D*U
    y_obs(j) = C(1,1) * x1_obs(j) + C(1,2) * x2_obs(j) + D * uu(j);

    % observador de ordem mínima
    % Equações dos integradores do observador
    eta1(j) = T * eta1_ponto(j - 1) + eta1(j - 1);
    
    % Equação diferencial de estados do observador
    eta1_ponto(j) = Achapeu * eta1(j) + Bchapeu * y_obs(j) + Fchapeu * uu(j);
    
    % Transformação
    % equação de saída do observador de ordem mínima
    xtil1(j) = Cchapeu(1) * eta1(j) + Dchapeu(1) * y_obs(j);
    xtil2(j) = Cchapeu(2) * eta1(j) + Dchapeu(2) * y_obs(j);
end


subplot(2,2,1)
plot(k*T,x1_obs,'--')
ylabel('x1')
title('Resposta ao Degrau do x1 Planta Original X Equação recursiva Observador')
subplot(2,2,2)
plot(k*T,x2_obs,'--')
ylabel('x2')
title('Resposta ao Degrau do x2 Planta Original X Equação recursiva Observador')

subplot(2,2,3)
plot(k*T,y_obs,'--')
ylabel('y')
legend('Planta Original', 'Observador')
title('Resposta ao Degrau da Planta Original X Equação recursiva Observador')

hold off



