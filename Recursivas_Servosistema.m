clear
clc
close all;

% Matrizes do Sistema
R1=51000;
R2=18000;
C1=0.0000001;
C2=0.00000068;

A=[0 1/C1;-1/(R1*R2*C2) (-1/(R1*C2))+(-1/(R2*C2))];
B=[0;1/(R1*R2*C2)];
C=[1 0];
D=0;

figure()
sistema=ss(A,B,C,D);
[Y,t,X]=step(sistema);

% separando os elementos do vetor X
x1a=X(:,1);
x1b=X(:,2);

figure(1)
subplot(2,2,1)
plot(t,x1a,'b')
hold on
subplot(2,2,2)
plot(t,x1b,'b')
hold on
subplot(2,2,3)
plot(t,Y,'b')
hold on

% autovalores desejados
s1 = -120.0000000000000 + 163.7251624610209i ;
s2 = -120.0000000000000 - 163.7251624610209i;
s3 = -2029.924353834096;
A_chapeu = [A zeros(2,1);-C 0];
B_chapeu = [B;0];

% teste de controlabilidade   = n+1 ?
rank([A B;-C 0])

% Matriz de ganhos do controlador K_chapeu = [k1 k2 -ki]
K_chapeu=acker(A_chapeu,B_chapeu,[s1 s2 s3])

% Verificando
eig(A_chapeu-B_chapeu*K_chapeu)

K=[K_chapeu(1) K_chapeu(2)];
Ki=-K_chapeu(3);
AA=[A-B*K B*Ki;-C 0];
BB=[0;0;1];
CC=[1 0 0];
DD=0;
sistema= ss(AA,BB,CC,DD);
% [Y,t,X] = step(sistema);
% x1a=X(:,1);
% x1b=X(:,2);
% 
% subplot(2,2,1)
% plot(t,x1a)
% hold on
% subplot(2,2,2)
% plot(t,x1b)
% hold on
% subplot(2,2,3)
% plot(t,Y)
% hold on


clear u
T = 1e-4;
k = 0:1:(0.15/T);
r=ones(1,length(k));

% condições iniciais
x1(1)= 0;  %  para k = 0
x2(1)= 0;  %  para k = 0
zeta(1)= 0;

u(1)=Ki*zeta(1);

x1_ponto(1)=B(1)*u(1);  %  para k = 0
x2_ponto(1)=B(2)*u(1);  %  para k = 0

y(1)=C(1)*x1(1)+C(2)*x2(1);  %  para k = 0

zeta_ponto(1)= r(1)-y(1);

for j=2:length(k)
    % Equações dos integradores
    x1(j)=T*x1_ponto(j-1)+x1(j-1);
    x2(j)=T*x2_ponto(j-1)+x2(j-1);    

    zeta(j)=T*zeta_ponto(j-1)+zeta(j-1);

    % U=-K*X+Ki*zeta
    u(j)=-(K(1)*x1(j)+K(2)*x2(j))+Ki*zeta(j);
   
    % Equação diferencial de estados:  Xponto=A*X+B*U
    x1_ponto(j)=A(1,1)*x1(j)+A(1,2)*x2(j)+B(1)*u(j);  
    x2_ponto(j)=A(2,1)*x1(j)+A(2,2)*x2(j)+B(2)*u(j); 
   
    % Equação de Saída: Y=C*X+D*U
    y(j)=C(1)*x1(j)+C(2)*x2(j);

    %Erro 
    zeta_ponto(j)= r(j)-y(j);


end

subplot(2,2,1)
plot(k*T,x1)
ylabel('x1')
subplot(2,2,2)
plot(k*T,x2)
ylabel('x2')
subplot(2,2,3)
plot(k*T,y)
ylabel('y')
legend('Planta Original', 'Servosistema')

hold off