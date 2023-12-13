clc
clear
close all;

%Equações recursivas para do sistema

% Matrizes do Sistema
R1=51000;
R2=18000;
C1=0.0000001;
C2=0.00000068;

A=[0 1/C1;-1/(R1*R2*C2) (-1/(R1*C2))+(-1/(R2*C2))];
B=[0;1/(R1*R2*C2)];
C=[1 0];
D=0;

sistema=ss(A,B,C,D);

% Matriz do observador
Ke = [129.4655901576317;0.0010876046533];

% Condições Iniciais


figure()
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

% simulação com discretização dos integradores
% I(z)=T/(z-1)=Xi(z)/Xi_ponto(z)
% xi(k)=T*xi_ponto(k-1)+xi(k-1)

clear u
T = 1e-4;
k=0:0.15/T;
u=ones(1,length(k));

% condições iniciais
x1(1)=0;  %  para k = 0
x2(1)=0;  %  para k = 0
x1_ponto(1)=B(1);  %  para k = 0
x2_ponto(1)=B(2);  %  para k = 0

y(1)=C(1)*x1(1)+C(2)*x2(1)+D*u(1);  %  para k = 0

% condições iniciais do observador
x1_obs(1)=0;  %  para k = 0
x2_obs(1)=0;  %  para k = 0

y_obs(1)=C(1)*x1_obs(1)+C(2)*x2_obs(1)+D*u(1)  %  para k = 0

x1_ponto_obs(1)=B(1)*u(1)+Ke(1)*(y(1)-y_obs(1));  %  para k = 0
x2_ponto_obs(1)=B(2)*u(1)+Ke(2)*(y(1)-y_obs(1));  %  para k = 0

for j=2:length(k)
    % Equações dos integradores
    x1(j)=T*x1_ponto(j-1)+x1(j-1);
    x2(j)=T*x2_ponto(j-1)+x2(j-1);    
    
   
    % Equação diferencial de estados:  Xponto=A*X+B*U
    x1_ponto(j)=A(1,1)*x1(j)+A(1,2)*x2(j)+B(1)*u(j);  
    x2_ponto(j)=A(2,1)*x1(j)+A(2,2)*x2(j)+B(2)*u(j); 
    % Equação de Saída: Y=C*X+D*U
    y(j)=C(1)*x1(j)+C(2)*x2(j)+D*u(j);
    
    % Equações dos integradores do observador
    x1_obs(j)=T*x1_ponto_obs(j-1)+x1_obs(j-1);
    x2_obs(j)=T*x2_ponto_obs(j-1)+x2_obs(j-1);    
    
     % Equação de Saída observador: Y=C*X+D*U
     y_obs(j)=C(1)*x1_obs(j)+C(2)*x2_obs(j)+D*u(j);
    
    % Equação diferencial de estados do observador:
    % Xponto=A*X+B*U+Ke(y-C*x)
    x1_ponto_obs(j)=A(1,1)*x1_obs(j)+A(1,2)*x2_obs(j)+B(1)*u(j)+Ke(1)*(y(j)-y_obs(j));  
    x2_ponto_obs(j)=A(2,1)*x1_obs(j)+A(2,2)*x2_obs(j)+B(2)*u(j)+Ke(2)*(y(j)-y_obs(j)); 
   
end

subplot(2,2,1)
plot(k*T,x1_obs,'*')
ylabel('x1')
subplot(2,2,2)
plot(k*T,x2_obs,'*')
ylabel('x2')
subplot(2,2,3)
plot(k*T,y_obs,'*')
ylabel('y')
legend('Planta Original', 'Observador')



