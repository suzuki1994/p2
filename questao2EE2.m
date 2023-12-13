clc;
clear;

% Matrizes de estados cálculadas 
A = [0, 1; -8, -6];
B = [0; 1];
C = [1, 0];
D = 0;
% Conversão para função de transferência
sys = ss(A, B, C, D);
figure (1)
step(sys)
G=tf(sys);
% Verificando polos da planta
autovalores_planta = eig(A);
polos_planta = autovalores_planta;

% Matrizes reduzidas
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

% Polo do observador, 5 vezes mais rápido do que os polos do controlador
L = -5 * abs(polos_planta(2))

% Ganho do observador
Ke = acker(Abb, Aab, L);
Ke = Ke(1) % Transformando em escalar

Achapeu = Abb - Ke * Aab
Bchapeu = Achapeu * Ke + Aba - Ke * Aaa
Fchapeu = Bb - Ke * Ba

Cchapeu = [0; 1]
Dchapeu = [1; Ke]

%recursiva 
T = 0.02;
k=0:round(7/T);
u=ones(1,length(k));

% definindo os vetores
x1(1)=0;
x2(1)=0;
x1_ponto(1)=0;
x2_ponto(1)=0;
eta1(1)=0;

eta1_ponto(1)=Fchapeu(1)*u(1);
xtil1(1)=0;
xtil2(1)=0;
% condições iniciais (podem ser diferentes de zero)
y(1) = C(1,1)*x1(1) + C(1,2)*x2(1) + D*u(1);  % para k = 0


for j=2:length(k)
    % sistema original
    % Equações dos integradores
    x1(j) = T * x1_ponto(j - 1) + x1(j - 1);
    x2(j) = T * x2_ponto(j - 1) + x2(j - 1);
    
    % Equação diferencial de estados:  Xponto=A*X+B*U
	x1_ponto(j) = A(1,1)*x1(j) + A(1,2)*x2(j) + B(1)*u(j);
    x2_ponto(j) = A(2,1)*x1(j) + A(2,2)*x2(j) + B(2)*u(j);  
    
    % Equação de Saída: Y=C*X+D*U
    y(j) = C(1,1) * x1(j) + C(1,2) * x2(j) + D * u(j);

    % observador de ordem mínima
    % Equações dos integradores do observador
    eta1(j) = T * eta1_ponto(j - 1) + eta1(j - 1);
    
    % Equação diferencial de estados do observador
    eta1_ponto(j) = Achapeu * eta1(j) + Bchapeu * y(j) + Fchapeu * u(j);
    
    % Transformação
    % equação de saída do observador de ordem mínima
    xtil1(j) = Cchapeu(1) * eta1(j) + Dchapeu(1) * y(j);
    xtil2(j) = Cchapeu(2) * eta1(j) + Dchapeu(2) * y(j);
end

figure (2);
plot(k * T, xtil1)
title('xtil da planta ao degrau unitário');
xlabel('t');
ylabel('y(t)');
legend('y(t)');
grid on;


%   % observador de ordem mínima
%     % Equações dos integradores do observador
%     eta1(j) = T * eta1_ponto(j - 1) + eta1(j - 1);
%     
%     % Equação diferencial de estados do observador
%     eta1_ponto(j) = Achapeu * eta1(j) + Bchapeu * y(j) + Fchapeu * u(j);
%     
%     % Transformação
%     % equação de saída do observador de ordem mínima
%     xtil1(j) = Cchapeu(1) * eta1(j) + Dchapeu(1) * y(j);
%     xtil2(j) = Cchapeu(2) * eta1(j) + Dchapeu(2) * y(j);
