function [g]=cacode(svnum)
% g2 contiene las combinaciones de taps para generar los 32 c�odigos de gps
g2=[2 6; 3 7; 4 8; 5 9; 1 9; 2 10; 1 8; 2 9; 3 10; 2 3; 3 4; 5 6; 6 7; ...
7 8; 8 9; 9 10; 1 4; 2 5; 3 6; 4 7; 5 8; 6 9; 1 3; 4 6; 5 7; 6 8; ...
7 9; 8 10; 1 6; 2 7; 3 8; 4 9];
a=g2(svnum,1);
b=g2(svnum,2);
x = -ones(1,10);
y = -ones(1,10);
g=zeros(1,1023);

for j = 1 : 1023
g1 = x(10);
g2 = y(a)*y(b);
g(j) = g1*g2; % Secuencia Gold: XOR de los codigos G1 y G2.
% x1 e y1 son las realimentaciones a los bloques de G1 y G2
x1 = x(3)*x(10);
y1 = y(2)*y(3)*y(6)*y(8)*y(9)*y(10);
% Corrimiento de los bits una posicion en los bloques G1 y G2.
% El bit que estaba en la posicion 10 se pierde.
x(2:10)=x(1:9);
y(2:10)=y(1:9);
x(1) = x1;
y(1) = y1;
end
