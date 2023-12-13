clc;
clear;
format long;

A = [0 1; -3 -2];
B = [0; 1];
C = [1 0];
D = 0;

s1=-2;
s2=-3;
s3=-15;
figure(1)
step(A,B,C,D)

A_chapeu = [A zeros(2,1);-C 0];
B_chapeu = [B;0];

rank([A B;-C 0]);

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
figure(2)
step(AA,BB,CC,DD)

%B
X0 = [0; 0; 0];
syms s
U=1/s;

I = eye(3);
X=inv(s*eye(3)-AA)*X0+inv(s*eye(3)-AA)*BB*U
Y=CC*X
y=ilaplace(Y)


t=0:.01:3.5;
y_expected=(5*exp(-3*t))/2 - (45*exp(-2*t))/13 - 0.038461538461538*exp(-15*t) + 1;

figure(3)

plot(t,y_expected)
hold on
u = ones(size(t));

sistema = ss(AA,BB,CC,DD); 
y2=lsim(sistema,u,t,X0);
plot(t,y2,'--')
hold off